#!/usr/bin/env nextflow

/*
 * verifybamid-py — contamination QC across an HPC (SLURM) cluster.
 *
 * Input: a TSV of  sample_id <tab> panel_s3_prefix <tab> cram_s3_path
 * (panel_s3_prefix + ".vcf.gz" is the project's population/genotype VCF, shared
 *  by every sample in the project).
 *
 * Smart bits:
 *   - the marker panel and the CHIP genotype matrix are built ONCE per unique
 *     panel (not per sample) — the central efficiency for 10,000s of samples
 *     that all share one project panel;
 *   - storeDir persists those artifacts so re-runs / future projects reuse them;
 *   - each CRAM is one independent SLURM task that STREAMS from S3 (no full
 *     download); -resume skips completed samples after any failure/preemption;
 *   - fast mode (downsampled panel + targeted .crai fetch) is the default to
 *     cut S3 egress to the cluster. MEASURED: ~13-15 GB per ~19 GB CRAM at the
 *     20k panel (markers are ~155 kb apart genome-wide, so targeted fetch still
 *     touches most CRAM slices) -- a ~25% saving over a full download, not the
 *     order-of-magnitude one might assume. At 100k samples this WAN egress, not
 *     CPU, is the binding constraint; see README "Scaling / egress".
 */

nextflow.enable.dsl = 2

// marker-set tag woven into cached artifact names (fast vs full panel).
// A function, not a top-level variable: strict DSL (NF >=25) forbids
// statements outside process/workflow/function declarations.
def mtag() { params.fast ? "fast${params.fast_n}" : "full" }

workflow {
    rows = Channel.fromPath(params.samples)
        | splitCsv(sep: '\t', strip: true)
        | map { r -> tuple(r[1], r[0], r[2]) }      // (panel_prefix, sample, cram)

    // unique panels -> (prefix, panel_name)
    panel_ch = rows.map { it[0] }.unique()
        .map { pfx -> tuple(pfx, pfx.tokenize('/').last()) }

    PANEL(panel_ch)                                  // (prefix, run_panel.parquet)

    if (params.chipmix) {
        CHIP(PANEL.out)                              // (prefix, chip.parquet)
        artifacts = PANEL.out.join(CHIP.out)         // (prefix, panel, chip)
    } else {
        no_chip = file("${projectDir}/assets/NO_CHIP")
        artifacts = PANEL.out.map { pfx, panel -> tuple(pfx, panel, no_chip) }
    }

    // broadcast each project's panel+chip to all its samples
    jobs = rows.combine(artifacts, by: 0)            // (prefix, sample, cram, panel, chip)

    VERIFYBAMID(jobs)

    // Merge per-sample rows into one table, single header. collectFile streams the
    // concatenation (no giant arg list / 100k files staged into one task like a
    // collect()+process would), and drops each file's header but the first.
    VERIFYBAMID.out.selfsm
        .collectFile(name: 'contamination.selfSM', storeDir: params.outdir,
                     keepHeader: true, skip: 1, sort: true)

    // Failed-sample manifest = requested ids minus ids that produced a .selfSM
    // (terminally-failed samples are 'ignore'd above, so they just never appear).
    all_ids = rows.map { it[1] }
        .collectFile(name: 'all_ids.txt', newLine: true, sort: true)
    ok_ids  = VERIFYBAMID.out.selfsm.map { it.baseName }
        .collectFile(name: 'ok_ids.txt', newLine: true, sort: true)
    REPORT(all_ids, ok_ids)
}

/*
 * Build the marker panel from the project VCF (once per panel, persisted).
 * In fast mode also downsamples to a sparse common-SNP panel; emits whichever
 * panel the per-sample step will actually use.
 */
process PANEL {
    tag "${pname}"
    storeDir params.panel_cache
    cpus 2
    memory '8 GB'
    time '2h'

    input:
    tuple val(prefix), val(pname)

    output:
    tuple val(prefix), path("${pname}.${mtag()}.panel.parquet")

    script:
    if (params.fast)
        """
        aws s3 cp ${prefix}.vcf.gz panel.vcf.gz --region ${params.region}
        ${params.bindir}/build-panel --vcf panel.vcf.gz --out full.parquet
        ${params.bindir}/downsample --panel full.parquet \
            --out ${pname}.${mtag()}.panel.parquet -n ${params.fast_n} --min-maf ${params.min_maf}
        """
    else
        """
        aws s3 cp ${prefix}.vcf.gz panel.vcf.gz --region ${params.region}
        ${params.bindir}/build-panel --vcf panel.vcf.gz --out ${pname}.${mtag()}.panel.parquet
        """
}

/*
 * Precompute the per-sample genotype matrix (once per panel, persisted), aligned
 * to the run panel so CHIPMIX is a single-column lookup per CRAM.
 */
process CHIP {
    tag "${pname}"
    storeDir params.panel_cache
    cpus 2
    memory '8 GB'
    time '2h'

    input:
    tuple val(prefix), path(panel)

    output:
    tuple val(prefix), path("${pname}.${mtag()}.chip.parquet")

    script:
    pname = prefix.tokenize('/').last()
    """
    aws s3 cp ${prefix}.vcf.gz panel.vcf.gz --region ${params.region}
    ${params.bindir}/build-chip --vcf panel.vcf.gz --panel ${panel} \
        --out ${pname}.${mtag()}.chip.parquet
    """
}

/*
 * Per-CRAM contamination estimate: stream from S3, pileup, FREEMIX/FREELK (+CHIPMIX).
 * One SLURM task per sample.
 */
process VERIFYBAMID {
    tag "${sample}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.selfSM'
    cpus params.cpus
    memory params.mem
    time params.time
    maxForks params.max_streams            // S3-stream throttle (see config)
    // Retry transient failures (S3 throttle -> coverage-guard exit, preemption);
    // after that, IGNORE so one bad CRAM can't abort a 100k-sample batch. Ignored
    // samples emit no .selfSM and are listed in failed_samples.txt by REPORT.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(prefix), val(sample), val(cram), path(panel), path(chip)

    output:
    path "${sample}.selfSM", emit: selfsm

    script:
    def chip_arg = chip.name == 'NO_CHIP' ? '' :
                   (params.best ? "--chip-matrix ${chip} --best" : "--chip-matrix ${chip} --chip-id ${sample}")
    def fast_arg = params.fast ? "--max-span ${params.max_span}" : ''
    """
    export AWS_REGION=${params.region}
    ${params.bindir}/verifybamid \
        --cram ${cram} \
        --ref ${params.ref} \
        --panel ${panel} \
        ${chip_arg} ${fast_arg} \
        --seq-id ${sample} \
        --jobs ${task.cpus} \
        --expires ${params.presign_expires} \
        --out ${sample}
    """
}

/*
 * Write failed_samples.txt = requested ids that produced no .selfSM (terminal
 * failures after retries). Empty file when every sample succeeded.
 */
process REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 1
    memory '1 GB'
    time '15m'

    input:
    path all_ids
    path ok_ids

    output:
    path "failed_samples.txt"

    script:
    """
    set -e
    sort -u ${all_ids} > a.txt
    sort -u ${ok_ids}  > o.txt
    comm -23 a.txt o.txt > failed_samples.txt || true
    echo "[report] \$(wc -l < failed_samples.txt) sample(s) failed after retries" >&2
    """
}
