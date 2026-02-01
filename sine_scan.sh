#!/usr/bin/env bash
set -euo pipefail
set +m

###############################################################################
# SINE FRAGMENT SEARCH (chunked DB, multi-query per chunk, PARALLELIZED)
#
# Usage:
#   TARGET_GENOME=/path/to/original.genome.fa ./sine_scan.sh <search_db.fa> <queries.fa>
#
# QUICK (default):
#   - scans only chunk 1
#   - uses only first query
#   - limits chunk FASTA building to first MAX_DB_CONTIGS contigs from DB.fai
#
# FULL:
#   FULL=1 TARGET_GENOME=... ./sine_scan.sh <search_db.fa> <queries.fa>
#
# Important for minus-bank:
#   DB headers like >Scaffold_1:0-78068()
#   The script maps hits back to original scaffold coords (Scaffold_1)
#   so bedtools MUST use TARGET_GENOME.fai (original genome).
#
# Requires:
#   samtools bedtools ssearch36 awk sort wc tee parallel
###############################################################################

########################
# CONFIG (edit if needed)
########################
CHUNK_BP="${CHUNK_BP:-100000000}"
FLANK="${FLANK:-50}"
MIN_ID="${MIN_ID:-65}"
MIN_COV="${MIN_COV:-0.90}"
MAX_CONCURRENT="${MAX_CONCURRENT:-$(nproc 2>/dev/null || echo 8)}"

# QUICK defaults
FULL="${FULL:-0}"                    # set FULL=1 for full run
MAX_DB_CONTIGS="${MAX_DB_CONTIGS:-50}"  # only used in QUICK mode when building chunk FASTA
STOP_AFTER_CHUNK1="${STOP_AFTER_CHUNK1:-1}"  # only used in QUICK mode

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

########################
# --chunk mode: run exactly one chunk (called by GNU parallel)
########################
if [[ "${1:-}" == "--chunk" ]]; then
  chunk_i="${2:-}"
  [[ "$chunk_i" =~ ^[0-9]+$ ]] || { echo "ERROR: --chunk requires numeric chunk index" >&2; exit 2; }

  : "${SEARCH_DB:?missing env SEARCH_DB}"
  : "${DB_SIZE:?missing env DB_SIZE}"
  : "${QFA:?missing env QFA}"
  : "${QLEN_TSV:?missing env QLEN_TSV}"
  : "${OUTDIR:?missing env OUTDIR}"
  : "${TOTAL_CHUNKS:?missing env TOTAL_CHUNKS}"
  : "${CHUNK_BP:?missing env CHUNK_BP}"
  : "${MIN_ID:?missing env MIN_ID}"
  : "${MIN_COV:?missing env MIN_COV}"
  : "${MAX_DB_CONTIGS:?missing env MAX_DB_CONTIGS}"

  offset=$(( (chunk_i-1)*CHUNK_BP + 1 ))
  if (( offset > DB_SIZE )); then
    log "Chunk $chunk_i/$TOTAL_CHUNKS skipped (offset $offset > db size $DB_SIZE)"
    exit 0
  fi

  remaining=$CHUNK_BP
  if (( offset + remaining - 1 > DB_SIZE )); then
    remaining=$(( DB_SIZE - offset + 1 ))
  fi
  endpos=$(( offset + remaining - 1 ))

  log "Starting chunk $chunk_i/$TOTAL_CHUNKS (positions $offset-$endpos)"

  CHUNK_FASTA="$OUTDIR/chunk_${chunk_i}.fa"
  rm -f "$CHUNK_FASTA"

  # Build chunk FASTA by walking the DB .fai
  # QUICK-mode can cap how many DB contigs we pull into this chunk via MAX_DB_CONTIGS.
  awk -v off="$offset" -v rem="$remaining" -v maxc="$MAX_DB_CONTIGS" -v do_cap="${DO_CAP_DB_CONTIGS:-0}" '
    BEGIN{FS=OFS="\t"; c=0}
    {
      if (do_cap==1 && c>=maxc) exit
      if(off>$2){off-=$2; next}
      take=$2-off+1; if(take>rem) take=rem
      print $1 ":" off "-" (off+take-1)
      c++
      rem-=take; off=1
      if(rem<=0) exit
    }
  ' "$SEARCH_DB.fai" | while read -r r; do
      samtools faidx "$SEARCH_DB" "$r" >> "$CHUNK_FASTA"
  done

  CHUNK_HITS="$OUTDIR/chunk_${chunk_i}.hits.bed"
  : > "$CHUNK_HITS"

  # ssearch36 m8C:
  # $1 qseqid, $2 sseqid, $3 pident, $4 length, $5 mismatch, $6 gapopen,
  # $7 qstart, $8 qend, $9 sstart, $10 send, $11 evalue, $12 bitscore
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
      ($7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/) {next}
      ($9 !~ /^[0-9]+$/ || $10 !~ /^[0-9]+$/) {next}
      ($3+0 < MINID) {next}

      {
        q=$1
        split(q,qa,/[ \t\r\n]+/); q=qa[1]
        ql = (q in qlen ? qlen[q] : 0)
        alen=$4+0
        if(ql>0 && (alen/ql) < MINCOV) next

        subj=$2
        split(subj,sa,/[ \t\r\n]+/); subj=sa[1]

        # SUBJECT HEADER PARSE (minus-bank + chunk headers)
        # Works for:
        #   Scaffold_1:0-78068()
        #   Scaffold_1:0-78068():18400-18449
        #   Scaffold_8:0-77793632():1-74852821
        #
        # Split-based (NOT regex PCRE) so it actually works in awk.
        n = split(subj, p, ":")
        scf = p[1]

        interval_start0 = 0
        if (n >= 2 && p[2] ~ /^[0-9]+-/) {
          interval_start0 = p[2]
          sub(/-.*/, "", interval_start0)
          interval_start0 += 0
        }

        chunk_off1 = 1
        for (i=3; i<=n; i++) {
          if (p[i] ~ /^[0-9]+-[0-9]+$/) {
            tmp = p[i]
            sub(/-.*/, "", tmp)
            chunk_off1 = tmp + 0
          }
        }

        # Subject coords within the CHUNK (1-based, ascending)
        a=$9+0; b=$10+0

        # Convert to ORIGINAL scaffold BED coords (0-based half-open)
        bed_start = interval_start0 + chunk_off1 + a - 2
        bed_end   = interval_start0 + chunk_off1 + b - 1

        if(bed_start < 0) bed_start = 0
        if(bed_end <= bed_start) next

        # Strand from query direction
        qs=$7+0; qe=$8+0
        strand = (qe > qs) ? "+" : "-"

        qsafe=q
        gsub(/\|/,"_",qsafe)
        gsub(/:/,"_",qsafe)

        print scf, bed_start, bed_end, qsafe, $3, strand >> HITSFILE
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
[[ $# -eq 2 ]] || { echo "Usage: TARGET_GENOME=/path/to/original.genome.fa $0 <search_db.fa> <queries.fa>" >&2; exit 1; }

SEARCH_DB="$1"
QFA_IN="$2"

for f in "$SEARCH_DB" "$QFA_IN"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

for t in samtools bedtools ssearch36 awk sort wc tee date basename rm mkdir cat parallel seq; do
  command -v "$t" >/dev/null || { echo "ERROR: missing $t" >&2; exit 1; }
done

TARGET_GENOME="${TARGET_GENOME:-$SEARCH_DB}"
[[ -f "$TARGET_GENOME" ]] || { echo "ERROR: TARGET_GENOME not found: $TARGET_GENOME" >&2; exit 1; }

BASE="sine_search_out/$(basename "${SEARCH_DB%.*}")"
if (( FULL == 1 )); then
  OUTDIR="$BASE"
else
  OUTDIR="${BASE}.quick"
fi

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1

[[ -f "$SEARCH_DB.fai" ]] || samtools faidx "$SEARCH_DB"
[[ -f "$TARGET_GENOME.fai" ]] || samtools faidx "$TARGET_GENOME"

DB_SIZE=$(awk '{s+=$2} END{print s+0}' "$SEARCH_DB.fai")
[[ "$DB_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: failed to parse $SEARCH_DB.fai" >&2; exit 1; }

TOTAL_CHUNKS=$(( (DB_SIZE + CHUNK_BP - 1) / CHUNK_BP ))

log "Search DB   : $SEARCH_DB"
log "Queries FASTA: $QFA_IN"
log "Target genome for bedtools/getfasta: $TARGET_GENOME"
log "Genome size : $DB_SIZE bp"
log "Chunk size  : $CHUNK_BP bp"
log "Total chunks: $TOTAL_CHUNKS"
log "Parallel jobs: $MAX_CONCURRENT"
log "MIN_ID=$MIN_ID  MIN_COV=$MIN_COV  FLANK=$FLANK"
log "OUTDIR      : $OUTDIR"
if (( FULL == 1 )); then
  log "MODE        : FULL"
else
  log "MODE        : QUICK (chunk1 + first query only; cap chunk build to MAX_DB_CONTIGS=$MAX_DB_CONTIGS)"
fi

# QUICK: use first query only (create a tiny FASTA)
QFA="$QFA_IN"
if (( FULL != 1 )); then
  QFA="$OUTDIR/queries.first1.fa"
  awk '
    /^>/{
      if(seen){ exit }
      seen=1
    }
    { if(seen) print }
  ' "$QFA_IN" > "$QFA"
  [[ -s "$QFA" ]] || { echo "ERROR: failed to extract first query from $QFA_IN" >&2; exit 2; }
fi

NQUERIES=$(awk '/^>/{n++} END{print n+0}' "$QFA")
log "Queries used: $NQUERIES"

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

export SEARCH_DB DB_SIZE QFA QLEN_TSV OUTDIR TOTAL_CHUNKS CHUNK_BP MIN_ID MIN_COV MAX_DB_CONTIGS

SELF="$0"
if command -v readlink >/dev/null 2>&1; then
  SELF_ABS="$(readlink -f "$SELF" 2>/dev/null || true)"
else
  SELF_ABS=""
fi
[[ -n "${SELF_ABS:-}" ]] && SELF="$SELF_ABS"

if (( FULL == 1 )); then
  export DO_CAP_DB_CONTIGS=0
  log "Launching parallel processing of $TOTAL_CHUNKS chunks (-j $MAX_CONCURRENT)"
  seq 1 "$TOTAL_CHUNKS" | parallel --eta --progress -j "$MAX_CONCURRENT" "$SELF" --chunk {}
else
  export DO_CAP_DB_CONTIGS=1
  log "QUICK: running chunk 1 only (no GNU parallel fanout)"
  "$SELF" --chunk 1
  if (( STOP_AFTER_CHUNK1 == 1 )); then
    : # continue to downstream steps using whatever hits we got
  fi
fi

log "Collecting hits..."
cat "$OUTDIR"/chunk_*.hits.bed > "$ALL_HITS" 2>/dev/null || true
rm -f "$OUTDIR"/chunk_*.hits.bed "$OUTDIR"/chunk_*.fa

if [[ ! -s "$ALL_HITS" ]]; then
  log "No hits found. (Pipeline OK.)"
  exit 0
fi

log "Cleaning hits..."

CLEAN_HITS="$OUTDIR/all_hits.clean.bed"
awk 'BEGIN{OFS="\t"}
     NF>=6 && $1!="" && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>$2 {print}
' "$ALL_HITS" > "$CLEAN_HITS"

if [[ ! -s "$CLEAN_HITS" ]]; then
  log "No valid hits after cleaning. (Pipeline OK.)"
  exit 0
fi

HITS_TOTAL=$(wc -l < "$CLEAN_HITS" 2>/dev/null || echo 0)
Q_HIT=$(awk 'BEGIN{FS="\t"} {q[$4]=1} END{n=0; for(k in q)n++; print n+0}' "$CLEAN_HITS")
log "Summary: queries_with_hits=$Q_HIT/$NQUERIES  total_hits=$HITS_TOTAL"

# Safety net: normalize any composite chroms (should be 0 after fixed parser, but keep robust)
TO_GENOME_HITS="$OUTDIR/all_hits.to_genome.bed"
awk 'BEGIN{OFS="\t"}
{
  chrom=$1
  if (chrom !~ /:/) { print; next }

  n = split(chrom, p, ":")
  scf = p[1]

  interval_start0 = 0
  if (n >= 2 && p[2] ~ /^[0-9]+-/) {
    interval_start0 = p[2]
    sub(/-.*/, "", interval_start0)
    interval_start0 += 0
  }

  chunk_off1 = 1
  for (i=3; i<=n; i++) {
    if (p[i] ~ /^[0-9]+-[0-9]+$/) {
      tmp = p[i]
      sub(/-.*/, "", tmp)
      chunk_off1 = tmp + 0
    }
  }

  shift = interval_start0 + chunk_off1 - 1
  $1 = scf
  $2 = $2 + shift
  $3 = $3 + shift
  print
}' "$CLEAN_HITS" > "$TO_GENOME_HITS"

bad_left=$(awk -F'\t' '$1 ~ /:.*\(\)/ {n++} END{print n+0}' "$TO_GENOME_HITS")
log "Composite chroms remaining after normalization: $bad_left"

log "Merging loci (bedtools)..."

# Strand-aware merging on TARGET_GENOME coordinates
bedtools slop -b "$FLANK" -g "$TARGET_GENOME.fai" -i "$TO_GENOME_HITS" \
| LC_ALL=C sort -t $'\t' -k1,1 -k6,6 -k2,2n \
| bedtools sort \
| bedtools merge -s -c 4,5,6 -o distinct,max,distinct \
> "$OUTDIR/merged_loci.bed"

LOCI=$(wc -l < "$OUTDIR/merged_loci.bed" 2>/dev/null || echo 0)
log "Summary: merged_loci=$LOCI"

log "Extracting FASTA..."

bedtools getfasta -fi "$TARGET_GENOME" -bed "$OUTDIR/merged_loci.bed" -s -name \
> "$OUTDIR/merged_loci.fa"

awk '/^>/{h=substr($0,2); gsub(/\|/,"_",h); gsub(/:/,"_",h); print ">"h; next} {print}' \
  "$OUTDIR/merged_loci.fa" > "$OUTDIR/merged_loci.safe.fa"

log "==============================================="
log "DONE"
log "BED   : $OUTDIR/merged_loci.bed"
log "FASTA : $OUTDIR/merged_loci.fa"
log "FASTA (safe): $OUTDIR/merged_loci.safe.fa"
if (( FULL != 1 )); then
  log "QUICK OK -> run FULL with: FULL=1 TARGET_GENOME=... $0 $SEARCH_DB $QFA_IN"
fi
log "==============================================="
