#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash bin/sine_scan.sh genome.fa sine_db.fa outdir
#   bash bin/sine_scan.sh genome.fa fragments.fa outdir --fragments
#
# Outputs:
#   outdir/run.log
#   outdir/fragments.fa               (if input was DB)
#   outdir/all_hits.bed
#   outdir/all_hits.clean.bed
#   outdir/merged_loci.bed            (anchors; small FLANK)
#   outdir/merged_loci.fa
#   outdir/merged_loci.big.bed        (candidates; BIGFLANK)
#   outdir/merged_loci.big.fa
#   outdir/merged_loci.big.safe.fa    (header-sanitized for cons/subfam)

GENOME="${1:-}"
DB_OR_FRAGS="${2:-}"
OUTDIR="${3:-}"

MODE="${4:-}"  # optional: --fragments

[[ -n "$GENOME" && -n "$DB_OR_FRAGS" && -n "$OUTDIR" ]] || {
  echo "Usage: $0 <genome.fa> <sine_db.fa|fragments.fa> <outdir> [--fragments]" >&2
  exit 1
}

for t in samtools bedtools ssearch36 awk sort wc tee; do
  command -v "$t" >/dev/null || { echo "ERROR: missing tool: $t" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
exec > >(tee -a "$LOG") 2>&1
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# ----------------------------
# Tunables (env overrides)
# ----------------------------
CHUNK_BP="${CHUNK_BP:-100000000}"
MIN_ID="${MIN_ID:-65}"
MIN_COV="${MIN_COV:-0.90}"

FLANK="${FLANK:-50}"         # anchor context
BIGFLANK="${BIGFLANK:-500}"  # candidate context (full SINE-ish)

FRAG_LEN="${FRAG_LEN:-50}"
FRAG_STEP="${FRAG_STEP:-25}"

# ----------------------------
# Index genome
# ----------------------------
[[ -f "$GENOME" ]] || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
[[ -f "$GENOME.fai" ]] || samtools faidx "$GENOME"

GENOME_SIZE=$(awk '{s+=$2} END{print s+0}' "$GENOME.fai")
TOTAL_CHUNKS=$(( (GENOME_SIZE + CHUNK_BP - 1) / CHUNK_BP ))

log "Genome: $GENOME"
log "Genome size : $GENOME_SIZE bp"
log "Chunk size  : $CHUNK_BP bp"
log "Total chunks: $TOTAL_CHUNKS"
log "MIN_ID=$MIN_ID  MIN_COV=$MIN_COV  FLANK=$FLANK  BIGFLANK=$BIGFLANK"

# ----------------------------
# Prepare fragments
# ----------------------------
FRAGS="$OUTDIR/fragments.fa"

if [[ "$MODE" == "--fragments" ]]; then
  [[ -f "$DB_OR_FRAGS" ]] || { echo "ERROR: fragments not found: $DB_OR_FRAGS" >&2; exit 1; }
  cp -f "$DB_OR_FRAGS" "$FRAGS"
  log "Using provided fragments: $DB_OR_FRAGS"
else
  [[ -f "$DB_OR_FRAGS" ]] || { echo "ERROR: SINE DB not found: $DB_OR_FRAGS" >&2; exit 1; }
  log "Fragmenting DB (LEN=$FRAG_LEN STEP=$FRAG_STEP): $DB_OR_FRAGS"
  awk -v LEN="$FRAG_LEN" -v STEP="$FRAG_STEP" -v PREFIX="FRAG" \
    -f "$(dirname "$0")/make_fragments.awk" \
    "$DB_OR_FRAGS" > "$FRAGS"
fi

# Split fragments into per-query FASTA files (simple + reliable)
QDIR="$OUTDIR/queries"
mkdir -p "$QDIR"
rm -f "$QDIR"/q*.fa

awk '
  /^>/ { if(n){close(fn)}; n++; fn=sprintf("'"$QDIR"'/q%06d.fa",n) }
  { print > fn }
' "$FRAGS"

NQUERIES=$(ls "$QDIR"/q*.fa | wc -l)
log "Fragments: $NQUERIES"

# ----------------------------
# Collect hits
# ----------------------------
ALL_HITS="$OUTDIR/all_hits.bed"
: > "$ALL_HITS"

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

  qi=0
  for qfa in "$QDIR"/q*.fa; do
    ((++qi))

    qhdr=$(awk 'NR==1{gsub(/^>/,"");print}' "$qfa")
    qlen=$(awk 'NR>1{gsub(/[ \t\r\n]/,""); n+=length($0)} END{print n+0}' "$qfa")

    TMP="$OUTDIR/tmp.$$"
    mkdir -p "$TMP"

    # ssearch36 -m 8C output includes meta-lines; filter hard.
    # Also convert subject "scf:start-end" coords to true scf coords.
    ssearch36 -m 8C "$qfa" "$CHUNK_FASTA" \
    | awk -v Q="$qhdr" -v QLEN="$qlen" -v MINID="$MIN_ID" -v MINCOV="$MIN_COV" '
        BEGIN{OFS="\t"}
        /^#/ {next}
        $1=="Fields:" || $2=="Fields:" {next}
        NF<12 {next}

        ($3 !~ /^([0-9]+(\.[0-9]+)?)$/) {next}   # pident
        ($4 !~ /^[0-9]+$/) {next}                # alen
        ($9 !~ /^[0-9]+$/ || $10 !~ /^[0-9]+$/) {next} # sstart/send

        ($3+0 < MINID) {next}
        (QLEN>0 && ($4/QLEN) < MINCOV) {next}

        {
          subj=$2
          region_start=1
          scf=subj
          if (match(subj, /^([^:]+):([0-9]+)-([0-9]+)$/, m)) {
            scf=m[1]
            region_start=m[2]+0
          }

          a=$9+0; b=$10+0
          a2 = region_start + a - 1
          b2 = region_start + b - 1

          if(a2<b2){st=a2-1; en=b2; str="+"}
          else     {st=b2-1; en=a2; str="-"}

          if(st<0) st=0
          if(en<=st) next

          print scf,st,en,Q,$3,str
        }' > "$TMP/hits.bed"

    nhits=$(wc -l < "$TMP/hits.bed")
    if (( nhits > 0 )); then
      ex=$(awk 'NR==1{printf "%s:%d-%d,%s",$1,$2,$3,$6}' "$TMP/hits.bed")
      log "chunk $chunk_i/$TOTAL_CHUNKS : query $qi/$NQUERIES : $nhits hits ($ex)"
      cat "$TMP/hits.bed" >> "$ALL_HITS"
    fi

    rm -rf "$TMP"
  done

  log "=== Finished chunk $chunk_i / $TOTAL_CHUNKS ==="
  offset=$(( offset + CHUNK_BP ))
  ((chunk_i++))
done

[[ -s "$ALL_HITS" ]] || { log "No hits found"; exit 0; }

# Clean (paranoia guard)
CLEAN_HITS="$OUTDIR/all_hits.clean.bed"
awk 'BEGIN{OFS="\t"} NF>=6 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>$2 {print}' \
  "$ALL_HITS" > "$CLEAN_HITS"

[[ -s "$CLEAN_HITS" ]] || { log "No valid hits after cleaning"; exit 0; }

# ----------------------------
# Merge anchors
# ----------------------------
log "Merging anchors (FLANK=$FLANK)"

bedtools slop -b "$FLANK" -g "$GENOME.fai" -i "$CLEAN_HITS" \
| sort -k1,1 -k2,2n \
| bedtools merge -s -c 4,5 -o distinct,max \
> "$OUTDIR/merged_loci.bed"

bedtools getfasta -fi "$GENOME" -bed "$OUTDIR/merged_loci.bed" -s -name \
> "$OUTDIR/merged_loci.fa"

# ----------------------------
# Expand to candidates (full SINE context)
# ----------------------------
log "Expanding to candidates (BIGFLANK=$BIGFLANK)"

bedtools slop -b "$BIGFLANK" -g "$GENOME.fai" -i "$OUTDIR/merged_loci.bed" \
| sort -k1,1 -k2,2n \
| bedtools merge -s -c 4,5 -o distinct,max \
> "$OUTDIR/merged_loci.big.bed"

bedtools getfasta -fi "$GENOME" -bed "$OUTDIR/merged_loci.big.bed" -s -name \
> "$OUTDIR/merged_loci.big.fa"

# ----------------------------
# Sanitize FASTA headers for cons/subfam (removes pipes etc.)
# ----------------------------
awk -f "$(dirname "$0")/sanitize_fasta_headers.awk" \
  "$OUTDIR/merged_loci.big.fa" > "$OUTDIR/merged_loci.big.safe.fa"

log "DONE"
log "Anchors BED : $OUTDIR/merged_loci.bed"
log "Anchors FA  : $OUTDIR/merged_loci.fa"
log "Cand BED    : $OUTDIR/merged_loci.big.bed"
log "Cand FA     : $OUTDIR/merged_loci.big.fa"
log "Cand SAFE   : $OUTDIR/merged_loci.big.safe.fa"
