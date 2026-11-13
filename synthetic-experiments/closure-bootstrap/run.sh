#!/usr/bin/env bash
# Runs all parameter sets of this collection and writes OUT/keys.jsonl.gz and
# the summary.
#
# Environment variables:
#   SAGE      SageMath launcher (default: sage)
#   JOBS      parallel processes (default: number of CPUs)
#   KEYS_TOY  keys per small parameter set (default 100)
#   KEYS_CM   keys per shortened Classic McEliece set (default 20)
#   SKIP_CM   set to 1 to skip the Classic McEliece sets
#   OUT       output directory (default: results)
#
# Each job is one process on one core; jobs run in parallel.
set -euo pipefail
cd "$(dirname "$0")"

SAGE=${SAGE:-sage}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}
KEYS_TOY=${KEYS_TOY:-100}
KEYS_CM=${KEYS_CM:-20}
SKIP_CM=${SKIP_CM:-0}
OUT=${OUT:-results}

rm -rf runs/parts
mkdir -p runs/parts "$OUT"
JOBLIST=runs/jobs.txt
: > "$JOBLIST"

# family m n t mode extra keys keys-per-job
SETS=()
for p in "6 64 6" "7 112 10" "7 122 11" "8 192 16" "8 256 24"; do
  SETS+=("goppa $p saving 2 $KEYS_TOY 10")     # n-5t+1 closures
  SETS+=("goppa $p saving 0 $KEYS_TOY 10")     # n-5t-1 closures
  SETS+=("goppa $p proved 0 $KEYS_TOY 10")     # n-4t-1 closures
done
for p in "6 56 6" "7 99 10" "8 241 24"; do     # too few binary quadrics
  SETS+=("goppa $p saving 0 $KEYS_TOY 10")
done
for p in "6 64 5" "7 120 9" "8 240 20"; do
  SETS+=("alternant $p proved 0 $KEYS_TOY 10")
done
if [ "$SKIP_CM" != 1 ]; then
  for p in "12 922 64" "13 1888 128"; do
    SETS+=("goppa $p saving 2 $KEYS_CM 1")
    SETS+=("goppa $p saving 0 $KEYS_CM 1")
  done
fi

for spec in "${SETS[@]}"; do
  read -r family m n t mode extra keys chunk <<< "$spec"
  lo=1
  while [ "$lo" -le "$keys" ]; do
    hi=$((lo + chunk - 1))
    [ "$hi" -gt "$keys" ] && hi=$keys
    echo "$SAGE -python src/closure_chain.py --family $family --m $m --n $n --t $t --mode $mode --extra $extra --seeds $lo-$hi --output runs/parts/$family-m$m-n$n-t$t-$mode-$extra--$lo.jsonl" >> "$JOBLIST"
    lo=$((hi + 1))
  done
done

echo "running $(wc -l < "$JOBLIST") jobs with $JOBS processes"
xargs -P "$JOBS" -I{} sh -c '{} > /dev/null' < "$JOBLIST"
cat runs/parts/*.jsonl | gzip -9 > "$OUT/keys.jsonl.gz"
python3 summarize.py "$OUT"
