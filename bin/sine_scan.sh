#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SINE FRAGMENT SEARCH (chunked genome, multi-query per chunk)
#
# Main improvement:
#   - ONE ssearch36 call per chunk using the whole multi-FASTA queries file
#     (no per-query splitting) -> huge reduction in process overhead.
#
# Fixes kept:
#   1) filter out non-hit/meta lines from ssearch36 (-m 8C) (e.g. "Fields:")
#   2) normalize subject name to scaffold and convert chunk-relative coords
#      to true scaffold coords using region start parsed from subject header
#
# Name sanitization:
#   - sanitize query IDs in BED/FASTA for downstream tools (EMBOSS cons):
#     replace '|' and ':' with '_'
###############################################################################

########################
# CONFIG
########################
CHUNK_BP=100000000
FLANK=50
MIN_ID=65
MIN_COV=0.90

########################
# INPUT
########################
[[ $# -eq 2 ]] || {
  echo "Usage: $0 <genome.fa> <queries.fa>" >&2
  exit 1
}

GENOME="$1"
QFA="$2"

for f in "$GENOME" "$QFA"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

########################
# TOOLS
########################
for t in samtools bedtools ssearch36 awk sort wc tee date basename rm mkdir cat; do
  command -v "$t" >/dev/null || { echo "ERROR: missing $t" >&2; exit 1; }
done

########################
# OUTPUT
########################
OUTDIR="sine_search_out/$(basename "${GENOME%.*}")"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

########################
# INDEX + GENOME SIZE
########################
[[ -f "$GENOME.fai" ]] || samtools faidx "$GENOME"

GENOME_SIZE=$(awk '{s+=$2} END{print s+0}' "$GENOME.fai")
[[ "$GENOME_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: failed to parse genome.fai" >&2; exit 1; }

TOTAL_CHUNKS=$(( (GENOME_SIZE + CHUNK_BP - 1) / CHUNK_BP ))

log "Genome size : $GENOME_SIZE bp"
log "Chunk size  : $CHUNK_BP bp"
log "Total chunks: $TOTAL_CHUNKS"
log "MIN_ID=$MIN_ID  MIN_COV=$MIN_COV  FLANK=$FLANK"

########################
# QUERY COUNT + QUERY LENGTHS (for coverage filter)
########################
NQUERIES=$(awk '/^>/{n++} END{print n+0}' "$QFA")
log "Queries in FASTA: $NQUERIES"

QLEN_TSV="$OUTDIR/query_lengths.tsv"
awk '
  BEGIN{ id=""; len=0 }
  /^>/{
    if(id!=""){ print id "\t" len }
    id=substr($0,2)
    # keep only first token as "id" (matches how most tools name queries)
    split(id,a,/[ \t\r\n]+/); id=a[1]
    len=0
    next
  }
  {
    gsub(/[ \t\r\n]/,"")
    len += length($0)
  }
  END{
    if(id!=""){ print id "\t" len }
  }
' "$QFA" > "$QLEN_TSV"

########################
# GLOBAL HIT COLLECTOR
########################
ALL_HITS="$OUTDIR/all_hits.bed"
: > "$ALL_HITS"

########################
# CHUNK LOOP
########################
chunk_i=1
offset=1

while (( offset <= GENOME_SIZE )); do
  log "=== Processing chunk $chunk_i / $TOTAL_CHUNKS ==="

  CHUNK_FASTA="$OUTDIR/chunk_${chunk_i}.fa"
  rm -f "$CHUNK_FASTA"

  remaining=$CHUNK_BP
  awk -v off="$offset" -v rem="$remaining" '
    BEGIN{FS=OFS="\t"}
    {
      if(off>$2){off-=$2; next}
      take=$2-off+1; if(take>rem) take=rem
      print $1 ":" off "-" (off+take-1)
      rem-=take; off=1
      if(rem<=0) exit
    }
  ' "$GENOME.fai" | while read -r r; do
      samtools faidx "$GENOME" "$r" >> "$CHUNK_FASTA"
  done

  TMP="$OUTDIR/tmp.$$"
  mkdir -p "$TMP"

  # ONE ssearch36 call per chunk (multi-query)
  ssearch36 -m 8C "$QFA" "$CHUNK_FASTA" \
  | awk -v QLENFILE="$QLEN_TSV" -v MINID="$MIN_ID" -v MINCOV="$MIN_COV" '
      BEGIN{
        OFS="\t"
        # load query lengths
        while((getline < QLENFILE) > 0){
          if(NF>=2) qlen[$1]=$2+0
        }
        close(QLENFILE)
      }

      # skip meta / non-hit lines
      /^#/ {next}
      $1=="Fields:" || $2=="Fields:" {next}
      NF<12 {next}

      # require numeric fields used (pident, alen, sstart, send)
      ($3 !~ /^([0-9]+(\.[0-9]+)?)$/) {next}
      ($4 !~ /^[0-9]+$/) {next}
      ($9 !~ /^[0-9]+$/ || $10 !~ /^[0-9]+$/) {next}

      # thresholds
      ($3+0 < MINID) {next}

      {
        q=$1
        # match qlen map keying (first token)
        split(q,qa,/[ \t\r\n]+/); q=qa[1]
        ql = (q in qlen ? qlen[q] : 0)
        alen=$4+0
        if(ql>0 && (alen/ql) < MINCOV) next

        subj=$2

        # subject is chunk header scf:start-end -> convert to true scaffold coords
        region_start=1
        scf=subj
        if (match(subj, /^([^:]+):([0-9]+)-([0-9]+)$/, m)) {
          scf=m[1]
          region_start=m[2]+0
        }

        a=$9+0; b=$10+0

        # convert to scaffold coordinates (1-based)
        a2 = region_start + a - 1
        b2 = region_start + b - 1

        if(a2<b2){st=a2-1; en=b2; str="+"}
        else     {st=b2-1; en=a2; str="-"}

        if(st<0) st=0
        if(en<=st) next

        # sanitize query ID for downstream tools (EMBOSS cons hates |)
        qsafe=q
        gsub(/\|/,"_",qsafe)
        gsub(/:/,"_",qsafe)

        # BED6-ish: chrom, start, end, name, pident, strand
        print scf, st, en, qsafe, $3, str
      }
    ' > "$TMP/hits.bed"

  nhits=$(wc -l < "$TMP/hits.bed")
  if (( nhits > 0 )); then
    ex=$(awk 'NR==1{printf "%s:%d-%d,%s",$1,$2,$3,$6}' "$TMP/hits.bed")
    log "chunk $chunk_i/$TOTAL_CHUNKS : $nhits hits ($ex)"
    cat "$TMP/hits.bed" >> "$ALL_HITS"
  else
    log "chunk $chunk_i/$TOTAL_CHUNKS : no hits"
  fi

  rm -rf "$TMP"

  log "=== Finished chunk $chunk_i / $TOTAL_CHUNKS ==="

  offset=$(( offset + CHUNK_BP ))
  ((chunk_i++))
done

########################
# FINAL MERGE
########################
[[ -s "$ALL_HITS" ]] || {
  log "No hits found in genome"
  exit 0
}

log "Merging loci"

CLEAN_HITS="$OUTDIR/all_hits.clean.bed"
awk 'BEGIN{OFS="\t"}
     NF>=6 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>$2 {print}
' "$ALL_HITS" > "$CLEAN_HITS"

[[ -s "$CLEAN_HITS" ]] || {
  log "No valid hits after cleaning"
  exit 0
}

bedtools slop -b "$FLANK" -g "$GENOME.fai" -i "$CLEAN_HITS" \
| sort -k1,1 -k2,2n \
| bedtools merge -s -c 4,5 -o distinct,max \
> "$OUTDIR/merged_loci.bed"

########################
# FASTA EXTRACTION
########################
log "Extracting FASTA"

bedtools getfasta \
  -fi "$GENOME" \
  -bed "$OUTDIR/merged_loci.bed" \
  -s -name \
> "$OUTDIR/merged_loci.fa"

# Extra: sanitize FASTA headers again (belt-and-suspenders for EMBOSS cons/SubFam)
awk '/^>/{h=substr($0,2); gsub(/\|/,"_",h); gsub(/:/,"_",h); print ">"h; next} {print}' \
  "$OUTDIR/merged_loci.fa" > "$OUTDIR/merged_loci.safe.fa"

########################
# DONE
########################
log "==============================================="
log "DONE"
log "BED   : $OUTDIR/merged_loci.bed"
log "FASTA : $OUTDIR/merged_loci.fa"
log "FASTA (safe for cons): $OUTDIR/merged_loci.safe.fa"
log "==============================================="
