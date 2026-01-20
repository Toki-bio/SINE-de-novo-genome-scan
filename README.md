# SINE-de-novo-genome-scan
Perform a de novo genome search for SINEs using a database of known SINE sequences.

Find SINE candidates in any genome by:
1) fragmenting a SINE database into short windows (default 50 bp, 25 bp step),
2) searching each fragment against the genome with `ssearch36`,
3) merging hits into anchor loci,
4) expanding anchors to full candidate sequences (BIGFLANK),
5) sanitizing FASTA headers for EMBOSS `cons` / SubFam (no pipes, etc.).

## Requirements
- samtools
- bedtools
- FASTA36 (ssearch36)
- awk, sort, wc, bash

Optional downstream:
- EMBOSS `cons`
- SubFam

## Quick start

### A) Provide a full SINE DB (it will be fragmented)
```bash
bash bin/sine_scan.sh genome.fa sine_db.fa out
```bash

### B) Provide pre-fragmented queries
```bash
bash bin/sine_scan.sh genome.fa fragments.fa out --fragments
```bash
---

### Minimal “workflow” summary (what it does, in one paragraph)
- Build `.fai`, fragment DB (or use provided fragments), search each fragment against genome chunks with `ssearch36`, filter by identity+coverage, convert chunk-relative coords to genome coords, merge on strand into anchor loci, extract anchors, then expand and re-extract longer candidates, and finally sanitize headers for `cons/SubFam`.

Key outputs

out/merged_loci.bed and out/merged_loci.fa
Anchor loci (small context; default FLANK=50)

out/merged_loci.big.bed and out/merged_loci.big.fa
Candidate loci (expanded context; default BIGFLANK=500)

out/merged_loci.big.safe.fa
Same candidates with safe headers for cons/SubFam

Tuning (env vars)
CHUNK_BP=100000000   # genome chunk size
MIN_ID=65            # percent identity threshold
MIN_COV=0.90         # aligned_len / query_len threshold

FLANK=50             # context around raw hits for anchor merge
BIGFLANK=500         # context around anchors for candidate sequences

FRAG_LEN=50          # fragment length when fragmenting DB
FRAG_STEP=25         # step size when fragmenting DB


Example:

BIGFLANK=1000 MIN_COV=0.85 bash bin/sine_scan.sh genome.fa sine_db.fa out
