#!/bin/ksh
# Sort all winters (all models) by mean NOV-APR snowfall, lowest first.
# Output: tab-separated RESULTS/snowfall_sorted.txt with columns snowfall, year, model

RESDIR=/highres/svenk/SEKEstorylines_snow/RESULTS
OUT=${RESDIR}/snowfall_sorted.txt

{
    printf 'snowfall\tyear\tmodel\n'
    for f in "${RESDIR}"/result_*_gwl3.0.txt; do
        model=${f##*/result_}
        model=${model%_gwl3.0.txt}
        awk -v m="$model" '!/^#/ && NF >= 2 {printf "%s\t%s\t%s\n", $2, $1, m}' "$f"
    done | sort -t "$(printf '\t')" -k1,1g |
        LC_NUMERIC=C awk -F '\t' '{printf "%.2f\t%s\t%s\n", $1, $2, $3}'
} > "$OUT"
