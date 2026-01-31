#!/usr/bin/env bash
set -euo pipefail
set +m

###############################################################################
# SINE FRAGMENT SEARCH (chunked genome, multi-query per chunk, PARALLELIZED)
#
# Requires: samtools bedtools ssearch36 awk sort wc tee parallel
#
# Outputs:
#   merged_loci.bed
#   merged_loci.fa
#   merged_loci.safe.fa
#
# Notes:
# - Uses GNU parallel to run chunks concurrently.
# - Avoids exported bash functions by re-invoking this script in --chunk mode.
###############################################################################

########################
# CONFIG (edit if needed)
########################
CHUNK_BP=100000000
FLANK=50
MIN_ID=65
MIN_COV=0.90
MAX_CONCURRENT="${MAX_CONCURRENT:-$(nproc 2>/dev/null || echo 8)}"

########################
# LOG helper
########################
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

########################
# --chunk mode: run exactly one chunk (called by GNU parallel)
########################
if [[ "${1:-}" == "--chunk" ]]; then
  chunk_i="${2:-}"
  [[ "$chunk_i" =~ ^[0-9]+$ ]] || { echo "ERROR: --chunk requires numeric chunk index" >&2; exit 2; }

  # Expect required env vars
  : "${GENOME:?missing env GENOME}"
  : "${GENOME_SIZE:?missing env GENOME_SIZE}"
  : "${QFA:?missing env QFA}"
  : "${QLEN_TSV:?missing env QLEN_TSV}"
  : "${OUTDIR:?missing env OUTDIR}"
  : "${TOTAL_CHUNKS:?missing env TOTAL_CHUNKS}"
  : "${CHUNK_BP:?missing env CHUNK_BP}"
  : "${MIN_ID:?missing env MIN_ID}"
  : "${MIN_COV:?missing env MIN_COV}"

  offset=$(( (chunk_i-1)*CHUNK_BP + 1 ))
  if (( offset > GENOME_SIZE )); then
    log "Chunk $chunk_i/$TOTAL_CHUNKS skipped (offset $offset > genome size $GENOME_SIZE)"
    exit 0
  fi

  remaining=$CHUNK_BP
  (( offset + remaining - 1 > GENOME_SIZE )) && remaining=$(( GENOME_SIZE - offset + 1 ))
  endpos=$(( offset + remaining - 1 ))

  log "Starting chunk $chunk_i/$TOTAL_CHUNKS (positions $offset-$endpos)"

  CHUNK_FASTA="$OUTDIR/chunk_${chunk_i}.fa"
  rm -f "$CHUNK_FASTA"

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

  CHUNK_HITS="$OUTDIR/chunk_${chunk_i}.hits.bed"
  : > "$CHUNK_HITS"

  ssearch36 -m 8C "$QFA" "$CHUNK_FASTA" \
  | awk -v QLENFILE="$QLEN_TSV" -v MINID="$MIN_ID" -v MINCOV="$MIN_COV" -v HITSFILE="$CHUNK_HITS" '
      BEGIN{
        OFS="\t"
        while((getline < QLENFILE) > 0){
          if(NF>=2) qlen[$1]=$2+0
        }
        close(QLENFILE)
      }
      /^#/ {next}
      $1=="Fields:" || $2=="Fields:" {next}
      NF<12 {next}

      ($3 !~ /^([0-9]+(\.[0-9]+)?)$/) {next}
      ($4 !~ /^[0-9]+$/) {next}
      ($9 !~ /^[0-9]+$/ || $10 !~ /^[0-9]+$/) {next}

      ($3+0 < MINID) {next}

      {
        q=$1
        split(q,qa,/[ \t\r\n]+/); q=qa[1]
        ql = (q in qlen ? qlen[q] : 0)
        alen=$4+0
        if(ql>0 && (alen/ql) < MINCOV) next

        subj=$2

        region_start=1
        scf=subj
        if (match(subj, /^([^:]+):([0-9]+)-([0-9]+)$/, m)) {
          scf=m[1]
          region_start=m[2]+0
        }

        # Use QUERY coords for strand (correct way)
        qs=$7+0; qe=$8+0
        # Subject coords always ascending
        a=$9+0; b=$10+0
        a2 = region_start + a - 1
        b2 = region_start + b - 1
        st=a2-1; en=b2

        if(st<0) st=0
        if(en<=st) next

        # Strand from query direction
        strand = (qe - qs > 0) ? "+" : "-"

        qsafe=q
        gsub(/\|/,"_",qsafe)
        gsub(/:/,"_",qsafe)

        print scf, st, en, qsafe, $3, strand >> HITSFILE
      }
    '

  nhits=$(wc -l < "$CHUNK_HITS" 2>/dev/null || echo 0)
  if (( nhits > 0 )); then
    ex=$(awk 'NR==1{printf "%s:%d-%d,%s",$1,$2+1,$3,$6}' "$CHUNK_HITS")
    log "Chunk $chunk_i/$TOTAL_CHUNKS done: $nhits hits (e.g. $ex)"
  else
    log "Chunk $chunk_i/$TOTAL_CHUNKS done: no hits"
    rm -f "$CHUNK_HITS"
  fi

  rm -f "$CHUNK_FASTA"
  exit 0
fi

########################
# MAIN MODE
########################
[[ $# -eq 2 ]] || { echo "Usage: $0 <genome.fa> <queries.fa>" >&2; exit 1; }

GENOME="$1"
QFA="$2"

for f in "$GENOME" "$QFA"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

for t in samtools bedtools ssearch36 awk sort wc tee date basename rm mkdir cat parallel; do
  command -v "$t" >/dev/null || { echo "ERROR: missing $t" >&2; exit 1; }
done

OUTDIR="sine_search_out/$(basename "${GENOME%.*}")"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1

[[ -f "$GENOME.fai" ]] || samtools faidx "$GENOME"

GENOME_SIZE=$(awk '{s+=$2} END{print s+0}' "$GENOME.fai")
[[ "$GENOME_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: failed to parse genome.fai" >&2; exit 1; }

TOTAL_CHUNKS=$(( (GENOME_SIZE + CHUNK_BP - 1) / CHUNK_BP ))

log "Genome size : $GENOME_SIZE bp"
log "Chunk size  : $CHUNK_BP bp"
log "Total chunks: $TOTAL_CHUNKS"
log "Parallel jobs: $MAX_CONCURRENT"
log "MIN_ID=$MIN_ID  MIN_COV=$MIN_COV  FLANK=$FLANK"

NQUERIES=$(awk '/^>/{n++} END{print n+0}' "$QFA")
log "Queries in FASTA: $NQUERIES"

QLEN_TSV="$OUTDIR/query_lengths.tsv"
awk '
  BEGIN{ id=""; len=0 }
  /^>/{
    if(id!=""){ print id "\t" len }
    id=substr($0,2)
    split(id,a,/[ \t\r\n]+/); id=a[1]
    len=0
    next
  }
  { gsub(/[ \t\r\n]/,""); len += length($0) }
  END{ if(id!=""){ print id "\t" len } }
' "$QFA" > "$QLEN_TSV"

ALL_HITS="$OUTDIR/all_hits.bed"
: > "$ALL_HITS"

# export env vars for --chunk mode
export GENOME GENOME_SIZE QFA QLEN_TSV OUTDIR TOTAL_CHUNKS CHUNK_BP MIN_ID MIN_COV

# get absolute script path (best effort)
SELF="$0"
if command -v readlink >/dev/null 2>&1; then
  SELF_ABS="$(readlink -f "$SELF" 2>/dev/null || true)"
else
  SELF_ABS=""
fi
[[ -n "${SELF_ABS:-}" ]] && SELF="$SELF_ABS"

log "Launching parallel processing of $TOTAL_CHUNKS chunks (-j $MAX_CONCURRENT)"
seq 1 "$TOTAL_CHUNKS" | parallel --eta --progress -j "$MAX_CONCURRENT" "$SELF" --chunk {}

log "All chunks completed. Collecting hits..."
cat "$OUTDIR"/chunk_*.hits.bed > "$ALL_HITS" 2>/dev/null || true
rm -f "$OUTDIR"/chunk_*.hits.bed "$OUTDIR"/chunk_*.fa

[[ -s "$ALL_HITS" ]] || { log "No hits found in genome"; exit 0; }

log "Merging loci"

CLEAN_HITS="$OUTDIR/all_hits.clean.bed"
awk 'BEGIN{OFS="\t"}
     NF>=6 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>$2 {print}
' "$ALL_HITS" > "$CLEAN_HITS"

[[ -s "$CLEAN_HITS" ]] || { log "No valid hits after cleaning"; exit 0; }

HITS_TOTAL=$(wc -l < "$CLEAN_HITS")
Q_HIT=$(awk 'BEGIN{FS="\t"} {q[$4]=1} END{print length(q)+0}' "$CLEAN_HITS")
log "Summary: queries_with_hits=$Q_HIT/$NQUERIES  total_hits=$HITS_TOTAL"

# ONLY FIX: Add strand (column 6) to sort for proper strand-aware merging
bedtools slop -b "$FLANK" -g "$GENOME.fai" -i "$CLEAN_HITS" \
| sort -k1,1 -k6,6 -k2,2n | bedtools sort \
| bedtools merge -s -c 4,5,6 -o distinct,max,distinct \
> "$OUTDIR/merged_loci.bed"

LOCI=$(wc -l < "$OUTDIR/merged_loci.bed")
log "Summary: merged_loci=$LOCI"

log "Extracting FASTA"

bedtools getfasta -fi "$GENOME" -bed "$OUTDIR/merged_loci.bed" -s -name \
> "$OUTDIR/merged_loci.fa"

awk '/^>/{h=substr($0,2); gsub(/\|/,"_",h); gsub(/:/,"_",h); print ">"h; next} {print}' \
  "$OUTDIR/merged_loci.fa" > "$OUTDIR/merged_loci.safe.fa"

log "==============================================="
log "DONE"
log "BED   : $OUTDIR/merged_loci.bed"
log "FASTA : $OUTDIR/merged_loci.fa"
log "FASTA (safe for cons): $OUTDIR/merged_loci.safe.fa"
log "==============================================="
