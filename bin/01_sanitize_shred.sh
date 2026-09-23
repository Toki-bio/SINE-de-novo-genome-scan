#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# 01_sanitize_shred.sh
#
# Sanitize + unique + normalize SINEBase FASTA, then shred into fragments with
# coordinate-annotated headers.
#
# Usage:
#   bash 01_sanitize_shred.sh  <input_sine_db.fa>  <prefix>
#
# Output:
#   <prefix>.fragments.norm.fa
# ==============================================================================

IN="${1:-SINEBase.nr95.fa}"
PREFIX="${2:-$(basename "$IN" | sed 's/\.[^.]*$//')}"

FRAGLEN=50
STEP_5P=10
STEP_OTHER=25

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# =========================
# 0. Sanitize FASTA IDs, ensure uniqueness, uppercase + one-line
# =========================
# Keep only the first whitespace-delimited token as ID (same as your original)
# Sanitize only | and : (to avoid EMBOSS cons + coordinate parsing conflicts)
# Append _<count> to enforce uniqueness
IN_U="$TMP_DIR/${PREFIX}.uniq.fa"

seqkit seq -u -w 0 "$IN" \
| awk '
  BEGIN { count=0 }
  /^>/{
    count++
    name=substr($0,2)
    split(name, parts, /[[:space:]]+/)
    id=parts[1]
    gsub(/\|/,"_",id)
    gsub(/:/,"_",id)
    printf(">%s_%d\n", id, count)
    next
  }
  { print }
' > "$IN_U"

# =========================
# 1. Index + gzip (BED order safety)
# =========================
seqkit faidx "$IN_U"
gzip -c "$IN_U" > "$IN_U.gz"

# =========================
# 2. 5′ region: 1–150 bp (offset=0)
# =========================
seqkit subseq -w 0 -r 1:150 "$IN_U" \
| sed 's/^>/\0offset=0:/' > "$TMP_DIR/${PREFIX}.5p.fa"

# =========================
# 3. Middle region: 151 .. (L-100), only if L > 250 (offset=150)
# =========================
seqkit fx2tab -w 0 -n -l "$IN_U" \
| awk 'BEGIN{OFS="\t"}
       $2>250{
         id=$1; L=$2;
         s=150;      # BED: 0-based
         e=L-100;
         if(e>s) print id,s,e
       }' > "$TMP_DIR/mid.bed"

seqkit subseq -w 0 --bed "$TMP_DIR/mid.bed" "$IN_U.gz" \
| sed 's/^>/\0offset=150:/' > "$TMP_DIR/${PREFIX}.mid.fa"

# =========================
# 4. 3′ region: last 99 bp, only if L > 150 (offset=s)
# =========================
seqkit fx2tab -w 0 -n -l "$IN_U" \
| awk 'BEGIN{OFS="\t"}
       $2>150 {
         id=$1; L=$2;
         s=L-99;
         e=L;
         if(s<0) s=0;
         if(e>s) print id,s,e,"offset="s":"id
       }' > "$TMP_DIR/tail.bed"

bedtools getfasta -fi "$IN_U" -bed "$TMP_DIR/tail.bed" -nameOnly \
| seqkit seq -w 0 -u \
> "$TMP_DIR/${PREFIX}.3p.fa"

# =========================
# 5. Sliding-window fragmentation
# =========================
seqkit sliding -w 0 -W "$FRAGLEN" -s "$STEP_5P"    "$TMP_DIR/${PREFIX}.5p.fa"  > "$TMP_DIR/frags.5p.fa"
seqkit sliding -w 0 -W "$FRAGLEN" -s "$STEP_OTHER" "$TMP_DIR/${PREFIX}.mid.fa" > "$TMP_DIR/frags.mid.fa"
seqkit sliding -w 0 -W "$FRAGLEN" -s "$STEP_OTHER" "$TMP_DIR/${PREFIX}.3p.fa"  > "$TMP_DIR/frags.3p.fa"

# =========================
# 6. Normalize fragment headers with correct coordinates
# =========================
for r in 5p mid 3p; do
  awk -v R="$r" '
    /^>/{
      h=substr($0,2)

      # Parse: offset=N:ID:start-end  (or ID:start-end if offset not present)
      offset=0
      if (index(h, "offset=") == 1) {
        colon1 = index(h, ":")
        offset = substr(h, 8, colon1-8)
        h = substr(h, colon1+1)
      }

      # Now h = ID:start-end
      last_colon = 0
      for (i = length(h); i > 0; i--) {
        if (substr(h, i, 1) == ":") { last_colon = i; break }
      }
      if (last_colon == 0) next

      id = substr(h, 1, last_colon-1)
      coords = substr(h, last_colon+1)

      gsub(/_sliding$/, "", id)
      # mid-region fragments come from `seqkit subseq --bed`, which appends
      # "_<start>-<end>:<strand>" to the ID; the true coordinates are rebuilt below
      sub(/_[0-9]+-[0-9]+:[.+-]$/, "", id)

      dash_pos = index(coords, "-")
      if (dash_pos == 0) next

      frag_start = substr(coords, 1, dash_pos-1) + 0
      frag_end   = substr(coords, dash_pos+1) + 0

      actual_start = frag_start + offset
      actual_end   = frag_end + offset

      # Header format matches your established style:
      # >ID|REGION:START-END
      printf(">%s|%s:%d-%d\n", id, R, actual_start, actual_end)
      next
    }
    {print}
  ' "$TMP_DIR/frags.$r.fa" > "$TMP_DIR/frags.$r.norm.fa"
done

cat "$TMP_DIR/frags.5p.norm.fa" "$TMP_DIR/frags.mid.norm.fa" "$TMP_DIR/frags.3p.norm.fa" \
  > "${PREFIX}.fragments.norm.fa"

# =========================
# 7. Final sanity stats
# =========================
seqkit stat -w 0 "${PREFIX}.fragments.norm.fa"
echo "[OK] Output: ${PREFIX}.fragments.norm.fa"
