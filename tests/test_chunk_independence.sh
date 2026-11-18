#!/usr/bin/env bash
# The hits sine_scan.sh reports must not depend on CHUNK_BP, and must be decided by the raw-score cutoff (MIN_RAW).
# Synthetic 3 Mb random sequence with 120 planted copies of a 50 bp fragment, six per 150 kb piece, none within 2 kb of
# a chunk border: 60 strong copies at 88 % identity (6 substitutions in positions 6-45, raw score 196) and 60 weak ones
# at 74 % (13 substitutions, raw score 133, below the default MIN_RAW 170). Scanned with CHUNK_BP=250000 (12 chunks)
# and CHUNK_BP=3000000 (1 chunk). PASS = both runs report exactly the 60 strong copies: no weak copy, no chance hit
# elsewhere, identical hit sets. With ssearch36's default E-value gate (before 2026-11-18) a hit's fate depended on the
# chunk size (8001 vs 4044 hits on a real 20 Mb slice); with the gate fully open and no raw cutoff, random sequence
# gave chance hits (106 of 226 hits in an earlier version of this test).
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
SCAN="${SCAN:-$HERE/sine_scan.sh}"
T=$(mktemp -d "${TMPDIR:-$HOME/tmp}/chunktest.XXXX"); trap 'rm -rf "$T"' EXIT
cd "$T"
awk 'BEGIN{
  srand(7); split("A C G T",B," ");
  q="GTCAGCTAGGCTTACGATCGGATCCATGCAAGTCTAGCTAGGCTAACGTA";            # 50 bp query
  print ">frag|5p:1-50\n" q > "q.fa";
  L=3000000; s=""; for(i=1;i<=L;i++) g[i]=B[int(rand()*4)+1];
  n=0
  for(p=0;p<20;p++) for(k=0;k<6;k++){                                 # 6 copies in each 150 kb piece
    pos=p*150000+5000+k*24000                                         # never within 2 kb of a 250 kb border
    ns=(k%2==0)?6:13                                                  # even k strong (88 %), odd k weak (74 %)
    m=q; delete hit; for(j=1;j<=ns;){ x=int(rand()*40)+6; if(x in hit) continue; hit[x]=1; j++   # ns distinct substitutions in 6..45; the alignment still spans 1..50
      c=substr(m,x,1); do{r=B[int(rand()*4)+1]}while(r==c); m=substr(m,1,x-1) r substr(m,x+1) }
    for(j=1;j<=50;j++) g[pos+j]=substr(m,j,1); n++
    print "chr1\t" pos "\t" pos+50 > (ns==6 ? "strong.bed" : "weak.bed")
  }
  printf ">chr1\n" > "g.fa"; line=""
  for(i=1;i<=L;i++){ line=line g[i]; if(length(line)==80){ print line > "g.fa"; line="" } }
  if(line!="") print line > "g.fa"
}'
run(){ mkdir -p "$1"; ( cd "$1" && cp ../g.fa ../q.fa . && CHUNK_BP=$2 MAX_CONCURRENT=4 SSEARCH_THREADS=2 bash "$SCAN" g.fa q.fa > scan.log 2>&1 ) || true
       [ -f "$1"/sine_search_out/g/all_hits.to_genome.bed ] || { echo "FAIL: scan with CHUNK_BP=$2 wrote no hits file"; tail -5 "$1/scan.log"; exit 1; }
       cut -f1-4,6 "$1"/sine_search_out/g/all_hits.to_genome.bed | sort -u > "$1.hits"; }
run small 250000; run big 3000000
found(){ bedtools intersect -u -f 0.9 -a "$2" -b <(cut -f1-3 "$1.hits") | wc -l; }
other(){ bedtools intersect -v -a <(cut -f1-3 "$1.hits") -b strong.bed | wc -l; }
fs=$(found small strong.bed); fb=$(found big strong.bed); ws=$(found small weak.bed); wb=$(found big weak.bed)
os=$(other small); ob=$(other big); ns=$(wc -l < small.hits); nb=$(wc -l < big.hits)
echo "60 strong + 60 weak planted | CHUNK_BP 250000: $ns hits, strong $fs, weak $ws, other $os | CHUNK_BP 3000000: $nb hits, strong $fb, weak $wb, other $ob"
diff <(cat small.hits) <(cat big.hits) > d.txt && same=yes || same=no
echo "identical hit sets: $same"; head -5 d.txt
[ "$fs" -eq 60 ] && [ "$fb" -eq 60 ] && [ "$ws" -eq 0 ] && [ "$wb" -eq 0 ] && [ "$os" -eq 0 ] && [ "$ob" -eq 0 ] && [ "$same" = yes ] && echo PASS || { echo FAIL; exit 1; }
