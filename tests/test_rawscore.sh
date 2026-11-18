#!/usr/bin/env bash
# The raw score sine_scan.sh writes (column 7, computed from the BTOP) must equal ssearch36's own
# "Smith-Waterman score" for every alignment. Usage: tests/test_rawscore.sh <queries.fa> <subject.fa>
# Pairs the -m 8CB rows with the -m 0 alignments of two runs by query id and coordinates.
set -euo pipefail
Q=$1; S=$2
HERE=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
# the rawscore() function exactly as sine_scan.sh defines it
awk '/^      function rawscore/{p=1} p{print} p && /^      }$/{exit}' "$HERE/sine_scan.sh" > "$T/fn.awk"
grep -q 'function rawscore' "$T/fn.awk" || { echo "FAIL: rawscore() not found in sine_scan.sh"; exit 1; }
# two runs: ssearch36 36.3.8e prints no -m 0 alignments when combined with -m 8CB; alignments are paired by coordinates
ssearch36 -T 4 -E "1000000 1" -m 8CB "$Q" "$S" > "$T/out.txt" 2>/dev/null
ssearch36 -T 4 -E "1000000 1" -m 0 "$Q" "$S" > "$T/aln.txt" 2>/dev/null
export LC_ALL=C   # join needs the same collation as sort
# key = query id + qs-qe:ss-se (several bank fragments of one SINE hit the same place with the same coordinates,
# so coordinates alone pair different alignments); value = raw score. Keys seen more than once are dropped.
cp "$T/fn.awk" "$T/prog.awk"
echo 'BEGIN{OFS="\t"} !/^#/ && NF>=13 && $3 ~ /^[0-9.]+$/ {print $1"@"$7"-"$8":"$9"-"$10, rawscore($13)}' >> "$T/prog.awk"
awk -f "$T/prog.awk" "$T/out.txt" | awk -F'\t' '{c[$1]++; v[$1]=$2} END{for(k in c) if(c[k]==1) print k"\t"v[k]}' | sort > "$T/btop.tsv"
# ssearch36 alignments: query block " N>>>QUERYID ...", then "Smith-Waterman score: N; ... overlap (qs-qe:ss-se)"
awk '/^ *[0-9]+>>>/{q=$0; sub(/^ *[0-9]+>>>/,"",q); split(q,a," "); q=a[1]; next}
     /Smith-Waterman score:/ && match($0,/\([0-9]+-[0-9]+:[0-9]+-[0-9]+\)/){
       s=$0; sub(/.*Smith-Waterman score: /,"",s); sub(/;.*/,"",s)
       print q"@"substr($0,RSTART+1,RLENGTH-2)"\t"s }' "$T/aln.txt" \
 | awk -F'\t' '{c[$1]++; v[$1]=$2} END{for(k in c) if(c[k]==1) print k"\t"v[k]}' | sort > "$T/sw.tsv"
n=$(join -t $'\t' "$T/btop.tsv" "$T/sw.tsv" | wc -l)
bad=$(join -t $'\t' "$T/btop.tsv" "$T/sw.tsv" | awk -F'\t' '$2!=$3' | wc -l)
echo "alignments compared: $n; differing raw scores: $bad"
join -t $'\t' "$T/btop.tsv" "$T/sw.tsv" | awk -F'\t' '$2!=$3' | head -5
[ "$n" -gt 0 ] && [ "$bad" -eq 0 ] && echo PASS || { echo FAIL; exit 1; }
