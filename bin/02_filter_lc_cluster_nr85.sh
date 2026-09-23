#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 02_filter_lc_cluster_nr85.sh
#
# Low-complexity filter fragments + cluster to nr85 (centroids).
#
# Usage:
#   bash 02_filter_lc_cluster_nr85.sh  <prefix.fragments.norm.fa>  [threads]
#
# Outputs:
#   <prefix>.fragments.lc.fa
#   <prefix>.fragments.nr85.fa
#   <prefix>.lowcomplexity.log
# ==============================================================================

IN="${1:-SINEBase.fragments.norm.fa}"
THREADS="${2:-8}"

# derive prefix robustly from input filename
base="$(basename "$IN")"
PREFIX="${base%.fragments.norm.fa}"
[[ "$PREFIX" != "$base" ]] || PREFIX="${base%.*}"

OUT_LC="${PREFIX}.fragments.lc.fa"
OUT_NR85="${PREFIX}.fragments.nr85.fa"
LOG="${PREFIX}.lowcomplexity.log"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "[INFO] Low-complexity filtering started" > "$LOG"
date >> "$LOG"

# =========================
# 1. FASTA -> TSV (temp)
# =========================
seqkit fx2tab -w 0 "$IN" > "$TMP/fragments.tsv"

# =========================
# 2. Repeat-based filtering (your logic unchanged)
# =========================
awk '
BEGIN{
  FS=OFS="\t"
}

function is_homopolymer_trimer(t) {
  return (t=="AAA" || t=="TTT" || t=="GGG" || t=="CCC")
}

function trinuc_repeat(s,   i,t,reps) {
  s = toupper(s)
  for (i = 1; i <= length(s) - 17; i++) {  # >=18 bp
    t = substr(s, i, 3)
    if (t !~ /^[ACGT]{3}$/) continue
    if (is_homopolymer_trimer(t)) continue
    reps = 1
    while (substr(s, i + reps * 3, 3) == t) {
      reps++
      if (reps >= 6) return 1
    }
  }
  return 0
}

{
  id  = $1
  seq = $2

  # Homopolymer >15 bp
  if (match(seq, /[Aa]{16,}|[Tt]{16,}|[Gg]{16,}|[Cc]{16,}/)) {
    print id, "homopolymer>15bp" >> "'"$TMP"'/filtered.log"
  }
  # Dinucleotide repeats >15 bp (>=8 repeats)
  else if (match(tolower(seq),
      /(ac){8,}|(ca){8,}|(ag){8,}|(ga){8,}|(at){8,}|(ta){8,}|(cg){8,}|(gc){8,}|(ct){8,}|(tc){8,}|(gt){8,}|(tg){8,}/)) {
    print id, "dinucleotide>15bp" >> "'"$TMP"'/filtered.log"
  }
  # Trinucleotide repeats >15 bp (>=6 repeats, non-homopolymer)
  else if (trinuc_repeat(seq)) {
    print id, "trinucleotide>15bp" >> "'"$TMP"'/filtered.log"
  }
  else {
    print id >> "'"$TMP"'/pass.ids"
  }
}
' "$TMP/fragments.tsv"

# =========================
# 3. Extract passing fragments
# =========================
seqkit grep -w 0 -f "$TMP/pass.ids" "$IN" > "$OUT_LC"

# =========================
# 4. Logging (your style)
# =========================
in_n=$(seqkit stat -w 0 "$IN" | awk 'END{gsub(",","",$4);print $4}')
out_n=$(seqkit stat -w 0 "$OUT_LC" | awk 'END{gsub(",","",$4);print $4}')
removed=$((in_n - out_n))

{
  echo
  echo "[SUMMARY]"
  echo "Input fragments : $in_n"
  echo "Output fragments: $out_n"
  echo "Filtered        : $removed"
  echo
  echo "[FILTERED FRAGMENTS]"
  echo -e "fragment_id\treason"
  if [[ -s "$TMP/filtered.log" ]]; then
    sort "$TMP/filtered.log"
  else
    echo "(none)"
  fi
  echo
  echo "[FILTER COUNTS]"
  if [[ -s "$TMP/filtered.log" ]]; then
    awk '{c[$2]++} END{for(k in c) printf "%-20s %d\n",k,c[k]}' "$TMP/filtered.log"
  else
    echo "(none)"
  fi
  echo
  echo "[THRESHOLDS]"
  echo "  homopolymer   > 15 bp"
  echo "  dinucleotide  > 15 bp (>=8 repeats)"
  echo "  trinucleotide > 15 bp (>=6 repeats, non-homopolymer)"
} >> "$LOG"

echo "[INFO] Done. Output FASTA: $OUT_LC"
echo "[INFO] Log written to:    $LOG"

# =========================
# 5. Cluster to nr85 (your script logic, with PREFIX outputs)
# =========================
echo ""
echo "=== Fragment Clustering (85% identity) ==="
echo "Input:   $OUT_LC"
echo "Output:  $OUT_NR85"
echo "Threads: $THREADS"
echo ""

TMP_CENTROIDS="$(mktemp)"
trap 'rm -f "$TMP_CENTROIDS"' RETURN

vsearch \
  --cluster_fast "$OUT_LC" \
  --id 0.85 \
  --strand both \
  --centroids "$TMP_CENTROIDS" \
  --threads "$THREADS"

seqkit seq -w 0 -u "$TMP_CENTROIDS" > "$OUT_NR85"

input_count=$(seqkit stat -w 0 "$OUT_LC" | awk 'END{gsub(",","",$4);print $4}')
output_count=$(seqkit stat -w 0 "$OUT_NR85" | awk 'END{gsub(",","",$4);print $4}')

clustered=$((input_count - output_count))
clustered_pct=$(awk "BEGIN {printf \"%.1f\", 100*$clustered/$input_count}")

echo ""
echo "=== Clustering Results ==="
echo "Input fragments:      $input_count"
echo "Centroid fragments:   $output_count"
echo "Clustered/redundant:  $clustered ($clustered_pct%)"
echo ""
echo "Output file:"
echo "  nr85 FASTA:         $OUT_NR85"
