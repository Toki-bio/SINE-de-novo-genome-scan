#!/usr/bin/env bash
set -euo pipefail
set +m

###############################################################################
# SINE FRAGMENT SEARCH (chunked DB, multi-query per chunk, PARALLELIZED)
#
# Usage:
#   ./sine_scan.sh <search_db.fa> <queries.fa>
#
# For minus-bank DB headers like:
#   >Scaffold_8:0-77793632()
# the script maps hits back to ORIGINAL scaffold coordinates (Scaffold_8),
# so bedtools must use the ORIGINAL genome .fai:
#   TARGET_GENOME=/path/to/original.genome.fa ./sine_scan.sh minus_bank.fa queries.fa
#
# Requires:
#   samtools (>= 1.10, for faidx -r) bedtools ssearch36 awk sort wc tee parallel seq
#
# Env: MAX_CONCURRENT (parallel chunk jobs, default <= 8), SSEARCH_THREADS (threads per
#   ssearch36 job, default 4), TMPDIR (GNU parallel buffers; default ~/tmp, never /tmp)
# Optional orient QC (after merge, before getfasta): ORIENT_FILTER=1
#   ORIENT_MIN_TRUE_ID=65 ORIENT_DELTA_REJECT=25 ORIENT_DELTA_PASS=30 ORIENT_KEEP_WEAK=1
#
# Outputs (in OUTDIR):
#   all_hits.bed
#   all_hits.clean.bed
#   all_hits.to_genome.bed
#   merged_loci.bed
#   merged_loci.fa
#   merged_loci.safe.fa
###############################################################################

########################
# CONFIG (edit if needed)
########################
CHUNK_BP="${CHUNK_BP:-100000000}"
FLANK="${FLANK:-50}"
MIN_ID="${MIN_ID:-65}"
MIN_COV="${MIN_COV:-0.90}"

# Do NOT default to nproc (steals the node). Cap by default.
MAX_CONCURRENT_DEFAULT="$(nproc 2>/dev/null || echo 8)"
if [[ "$MAX_CONCURRENT_DEFAULT" =~ ^[0-9]+$ ]] && (( MAX_CONCURRENT_DEFAULT > 8 )); then
  MAX_CONCURRENT_DEFAULT=8
fi
MAX_CONCURRENT="${MAX_CONCURRENT:-$MAX_CONCURRENT_DEFAULT}"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

########################
# --chunk mode (called by GNU parallel)
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

  # Build chunk FASTA by walking DB .fai (global offset across concatenated contigs)
  awk -v off="$offset" -v rem="$remaining" '
    BEGIN{FS=OFS="\t"}
    {
      if(off>$2){off-=$2; next}
      take=$2-off+1; if(take>rem) take=rem
      print $1 ":" off "-" (off+take-1)
      rem-=take; off=1
      if(rem<=0) exit
    }
  ' "$SEARCH_DB.fai" > "$CHUNK_FASTA.regions"
  # one samtools call for the whole region list: one call PER scaffold took >1 h per chunk
  # on fragmented assemblies (saq: hundreds of thousands of scaffolds of a few hundred bp)
  samtools faidx "$SEARCH_DB" -r "$CHUNK_FASTA.regions" > "$CHUNK_FASTA"
  rm -f "$CHUNK_FASTA.regions"

  CHUNK_HITS="$OUTDIR/chunk_${chunk_i}.hits.bed"
  : > "$CHUNK_HITS"

  # ssearch36 -m 8C columns:
  # $1 qseqid, $2 sseqid, $3 pident, $4 length,
  # $7 qstart, $8 qend, $9 sstart, $10 send, ...
  # -T is required: ssearch36 otherwise starts one thread per core in EVERY
  # parallel chunk job (13 jobs x 128 threads drove a shared node to load 699).
  ssearch36 -T "${SSEARCH_THREADS:-4}" -m 8C "$QFA" "$CHUNK_FASTA" \
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

        # Robust split-based subject parsing (works for BOTH):
        #   regular chunk headers:   Scaffold_1:100-200
        #   minus-bank + chunk:      Scaffold_8:0-77793632():1-74852821
        #
        # We always map back to ORIGINAL scaffold coords:
        #   scf = Scaffold_8
        #   shift = interval_start0 + chunk_off1 - 1
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

        shift = interval_start0 + chunk_off1 - 1

        # Regular FASTA: the subject name is the samtools region "scf:start-end",
        # whose start is 1-BASED (the minus-bank "()" form above is 0-based from
        # bedtools getfasta). Treating it as 0-based put every hit 1 bp to the right.
        if (n == 2 && p[2] ~ /^[0-9]+-[0-9]+$/) {
          start1 = p[2]; sub(/-.*/, "", start1)
          shift = start1 - 1
        }

        # Strand: ssearch36 (36.3.8) reports a reverse-complement match by REVERSING
        # THE QUERY coordinates ($7 > $8) with the subject ascending, so testing the
        # subject alone labelled every minus-strand hit "+". Minus = exactly one of the
        # two coordinate pairs is descending.
        s1=$9+0; s2=$10+0
        if (s1 <= s2) { ss=s1; ee=s2 } else { ss=s2; ee=s1 }
        strand = (((s1 > s2) + ($7+0 > $8+0)) == 1) ? "-" : "+"

        # Convert to ORIGINAL scaffold BED coords (0-based, half-open)
        bed_start = shift + (ss - 1)
        bed_end   = shift + ee
        if(bed_start < 0) bed_start = 0
        if(bed_end <= bed_start) next

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
# MAIN
########################
[[ $# -eq 2 ]] || { echo "Usage: $0 <search_db.fa> <queries.fa>" >&2; exit 1; }

SEARCH_DB="$1"
QFA="$2"

for f in "$SEARCH_DB" "$QFA"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

for t in samtools bedtools ssearch36 awk sort wc tee date basename rm mkdir cat parallel seq head; do
  command -v "$t" >/dev/null || { echo "ERROR: missing $t" >&2; exit 1; }
done

# Detect minus-bank by first header
first_hdr="$(awk 'BEGIN{h=""} /^>/{print $0; exit}' "$SEARCH_DB" | sed 's/\r$//')"
DB_IS_MINUS=0
if [[ "$first_hdr" =~ ^\>[^:]+:[0-9]+-[0-9]+[[:space:]]*\(\)\ *$ ]] || [[ "$first_hdr" =~ ^\>[^:]+:[0-9]+-[0-9]+\(\) ]]; then
  DB_IS_MINUS=1
fi

# Decide TARGET_GENOME
TARGET_GENOME="${TARGET_GENOME:-$SEARCH_DB}"

# If minus-bank detected but TARGET_GENOME not set (still equals DB), fail EARLY with instruction.
if (( DB_IS_MINUS == 1 )) && [[ "$TARGET_GENOME" == "$SEARCH_DB" ]]; then
  echo "ERROR: minus-bank style DB detected, but TARGET_GENOME is not set." >&2
  echo "Your DB headers contain coordinate wrappers (e.g. >Scaffold_8:0-...())." >&2
  echo "Hits are mapped to ORIGINAL scaffold names (Scaffold_8), so bedtools must use ORIGINAL genome .fai." >&2
  echo "" >&2
  echo "Run like this:" >&2
  echo "  TARGET_GENOME=/path/to/original.genome.fa MAX_CONCURRENT=8 $0 $SEARCH_DB $QFA" >&2
  exit 2
fi

[[ -f "$TARGET_GENOME" ]] || { echo "ERROR: TARGET_GENOME not found: $TARGET_GENOME" >&2; exit 1; }

OUTDIR="sine_search_out/$(basename "${SEARCH_DB%.*}")"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1

[[ -f "$SEARCH_DB.fai" ]] || samtools faidx "$SEARCH_DB"
[[ -f "$TARGET_GENOME.fai" ]] || samtools faidx "$TARGET_GENOME"

DB_SIZE=$(awk '{s+=$2} END{print s+0}' "$SEARCH_DB.fai")
[[ "$DB_SIZE" =~ ^[0-9]+$ ]] || { echo "ERROR: failed to parse $SEARCH_DB.fai" >&2; exit 1; }

TOTAL_CHUNKS=$(( (DB_SIZE + CHUNK_BP - 1) / CHUNK_BP ))

log "Search DB   : $SEARCH_DB"
log "Queries FASTA: $QFA"
log "Genome size : $DB_SIZE bp"
log "Chunk size  : $CHUNK_BP bp"
log "Total chunks: $TOTAL_CHUNKS"
log "Parallel jobs: $MAX_CONCURRENT"
log "MIN_ID=$MIN_ID  MIN_COV=$MIN_COV  FLANK=$FLANK"
if (( DB_IS_MINUS == 1 )); then
  log "DB detected : minus-bank style headers"
else
  log "DB detected : regular FASTA headers"
fi
log "MAP_MODE    : ORIG"
log "bedtools -g : $TARGET_GENOME.fai"
log "getfasta -fi: $TARGET_GENOME"
log "OUTDIR      : $OUTDIR"

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

SSEARCH_THREADS="${SSEARCH_THREADS:-4}"
log "ssearch36 threads per chunk: $SSEARCH_THREADS (x $MAX_CONCURRENT jobs)"
export SEARCH_DB DB_SIZE QFA QLEN_TSV OUTDIR TOTAL_CHUNKS CHUNK_BP MIN_ID MIN_COV SSEARCH_THREADS

SELF="$0"
if command -v readlink >/dev/null 2>&1; then
  SELF_ABS="$(readlink -f "$SELF" 2>/dev/null || true)"
else
  SELF_ABS=""
fi
[[ -n "${SELF_ABS:-}" ]] && SELF="$SELF_ABS"

log "Launching parallel processing of $TOTAL_CHUNKS chunks (-j $MAX_CONCURRENT)"
# GNU parallel buffers job output in its tmpdir; the default /tmp is often a small RAM
# tmpfs or root partition and filled up on real runs. Use $TMPDIR, else ~/tmp.
PAR_TMP="${TMPDIR:-$HOME/tmp}"
mkdir -p "$PAR_TMP"
seq 1 "$TOTAL_CHUNKS" | parallel --tmpdir "$PAR_TMP" --eta --progress -j "$MAX_CONCURRENT" "$SELF" --chunk {}

log "All chunks completed. Collecting hits..."
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

# Safety: if any composite chrom accidentally slipped in, normalize it here.
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

bedtools slop -b "$FLANK" -g "$TARGET_GENOME.fai" -i "$TO_GENOME_HITS" \
| LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k6,6 \
| bedtools merge -s -c 4,5,6 -o distinct,max,distinct \
> "$OUTDIR/merged_loci.bed"

LOCI=$(wc -l < "$OUTDIR/merged_loci.bed" 2>/dev/null || echo 0)
log "Summary: merged_loci=$LOCI"

LOCI_BED="$OUTDIR/merged_loci.bed"
if [[ "${ORIENT_FILTER:-0}" == "1" ]]; then
  ORIENT_DIR="$(cd "$(dirname "$SELF")" && pwd)"
  ORIENT_PY="$ORIENT_DIR/orient_filter.py"
  if [[ ! -f "$ORIENT_PY" ]]; then
    echo "ERROR: ORIENT_FILTER=1 but missing $ORIENT_PY" >&2
    exit 2
  fi
  for t in python3 mafft esl-alistat; do
    command -v "$t" >/dev/null || { echo "ERROR: ORIENT_FILTER needs $t" >&2; exit 1; }
  done
  log "Orientation filter (MAFFT + esl-alistat)..."
  python3 "$ORIENT_PY" \
    --queries "$QFA" \
    --genome "$TARGET_GENOME" \
    --hits "$TO_GENOME_HITS" \
    --loci "$LOCI_BED" \
    --out-ok "$OUTDIR/merged_loci.orient_ok.bed" \
    --out-fail "$OUTDIR/merged_loci.orient_fail.bed" \
    --tsv "$OUTDIR/orient_filter.tsv" \
    --min-true-id "${ORIENT_MIN_TRUE_ID:-65}" \
    --delta-reject "${ORIENT_DELTA_REJECT:-25}" \
    --delta-pass "${ORIENT_DELTA_PASS:-30}" \
    ${ORIENT_KEEP_WEAK:+--keep-weak} \
    --tmpdir "${TMPDIR:-$HOME/tmp}"
  LOCI_OK=$(wc -l < "$OUTDIR/merged_loci.orient_ok.bed" 2>/dev/null || echo 0)
  LOCI_FAIL=$(wc -l < "$OUTDIR/merged_loci.orient_fail.bed" 2>/dev/null || echo 0)
  log "Orient filter: kept=$LOCI_OK rejected=$LOCI_FAIL (see orient_filter.tsv)"
  LOCI_BED="$OUTDIR/merged_loci.orient_ok.bed"
fi

log "Extracting FASTA..."
bedtools getfasta -fi "$TARGET_GENOME" -bed "$LOCI_BED" -s -name \
> "$OUTDIR/merged_loci.fa"

awk '/^>/{h=substr($0,2); gsub(/\|/,"_",h); gsub(/:/,"_",h); print ">"h; next} {print}' \
  "$OUTDIR/merged_loci.fa" > "$OUTDIR/merged_loci.safe.fa"

log "==============================================="
log "DONE"
log "BED   : $OUTDIR/merged_loci.bed"
log "FASTA : $OUTDIR/merged_loci.fa"
log "FASTA (safe): $OUTDIR/merged_loci.safe.fa"
log "==============================================="
