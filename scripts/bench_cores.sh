#!/usr/bin/env bash
# Map the cores-per-sample efficient frontier for verifybamid-py IN YOUR ENVIRONMENT.
#
# The right cores-per-sample depends on the I/O regime: on local/in-region-fast storage
# the work is CPU-bound and scales ~linearly (cores are "free" for latency, neutral for
# throughput); on S3 over a high-latency link it is network-bound and extra cores mostly
# idle. So measure where you'll actually run (the HPC node pulling from S3), not on a dev
# box. Runs N replicates per core count and reports the MEDIAN to cut network variance.
#
# Decision metric: CPU-seconds/sample = cores * wall. Batch wall-time on a fixed slot
# budget is proportional to it, so the core count that MINIMISES cpu_s is the throughput
# sweet spot. wall alone is the latency axis. util% shows how much each core idles on I/O.
#
# Usage:
#   scripts/bench_cores.sh <cram> <ref> <panel> [crai] [reps] [core_list]
#   # local:  scripts/bench_cores.sh sample.cram ref.fa fast20k.parquet
#   # S3:     scripts/bench_cores.sh s3://b/sample.cram ref.fa fast20k.parquet "" 3 "1 2 4 8"
set -euo pipefail

CRAM=${1:?cram path or s3:// uri}
REF=${2:?reference fasta}
PANEL=${3:?panel parquet}
CRAI=${4:-}
REPS=${5:-3}
CORES=${6:-"1 2 4 8"}

BIN="$(cd "$(dirname "$0")/.." && pwd)/.venv/bin"
MAXSPAN=$([[ "$PANEL" == *fast* ]] && echo "--max-span 150000" || echo "")

# presign if streaming from S3 (pileup needs https, not s3://)
if [[ "$CRAM" == s3://* ]]; then
  CRAM=$(aws s3 presign "$CRAM" --expires-in 14400)
  [[ -z "$CRAI" ]] && CRAI=$(aws s3 presign "${4:?pass s3 crai or it is derived}" --expires-in 14400) || true
fi
CRAI_ARG=$([[ -n "$CRAI" ]] && echo "--crai $CRAI" || echo "")

printf "%-6s %-8s %-9s %-9s %-7s\n" cores wall_med cpu_med cpu_s util%
for c in $CORES; do
  walls=(); cpus=()
  for r in $(seq 1 "$REPS"); do
    read e u s < <( { /usr/bin/time -f '%e %U %S' \
        "$BIN/pileup" --cram "$CRAM" $CRAI_ARG --ref "$REF" --panel "$PANEL" \
        $MAXSPAN --jobs "$c" >/dev/null 2>/tmp/_bench.err; } 2>&1; tail -1 /tmp/_bench.err )
    # drop runs that under-read (degraded stream) so they don't pollute timing
    frac=$(grep -oE 'covered=[0-9,]+/[0-9,]+ \(([0-9.]+)\)' /tmp/_bench.err | grep -oE '\([0-9.]+\)$' | tr -d '()')
    awk "BEGIN{exit !(${frac:-0} >= 0.8)}" || { echo "  (skip c=$c rep=$r: covered=$frac)" >&2; continue; }
    walls+=("$e"); cpus+=("$(awk "BEGIN{print $u+$s}")")
  done
  [[ ${#walls[@]} -eq 0 ]] && { printf "%-6s ALL-DEGRADED\n" "$c"; continue; }
  wmed=$(printf '%s\n' "${walls[@]}" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  cmed=$(printf '%s\n' "${cpus[@]}"  | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
  cpu_s=$(awk "BEGIN{printf \"%.0f\", $c*$wmed}")
  util=$(awk "BEGIN{printf \"%.0f\", $cmed/$wmed/$c*100}")
  printf "%-6s %-8s %-9s %-9s %-7s\n" "$c" "$wmed" "$cmed" "$cpu_s" "$util"
done
echo "# pick the core count with the lowest cpu_s (throughput) ; lower wall = latency"
