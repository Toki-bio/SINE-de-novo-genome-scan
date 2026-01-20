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

# SINEBase fragment query builder (sanitize → shred → LC filter → nr85)

This repo builds a robust, search-ready fragment query set from a SINE database FASTA.

It produces:
- `*.fragments.norm.fa`  (shredded fragments with normalized coordinate headers)
- `*.fragments.lc.fa`    (low-complexity filtered fragments)
- `*.fragments.nr85.fa`  (nonredundant fragments clustered at 85% identity)

Designed to be “safe by default” for downstream tools (including EMBOSS `cons`), while keeping
the workflow simple and reproducible.

---

## Requirements

**Core**
- bash
- awk
- seqkit
- bedtools
- vsearch

(plus standard utils: `gzip`, `mktemp`)

---

## Pipeline overview

### Script 1 — sanitize + shred
`01_sanitize_shred.sh`

Input:
- SINE DB FASTA (e.g. `SINEBase.nr95.fa`)

Does:
1) Normalize sequences: uppercase + one-line (`seqkit seq -u -w 0`)
2) Sanitize FASTA **ID token only**: replaces `|` and `:` with `_`
3) Enforce unique IDs by appending `_<record_number>`
4) Shred each record into regions:
   - **5′**: 1–150 bp (sliding step 10)
   - **mid**: 151..(L-100) if L>250 (sliding step 25)
   - **3′**: last 99 bp if L>150 (sliding step 25)
5) Sliding window fragments (default `FRAGLEN=50`)
6) Normalize fragment headers to embed true coordinates within the original SINE record:

Header format:
```
>SEQID|REGION:START-END
```

Output:
- `<prefix>.fragments.norm.fa`

Run:
```bash
bash 01_sanitize_shred.sh SINEBase.nr95.fa SINEBase.nr95
```

---

### Script 2 — low-complexity filter + nr85 clustering
`02_filter_lc_cluster_nr85.sh`

Input:
- `<prefix>.fragments.norm.fa`

Does:
1) Filters obvious low-complexity fragments:
   - homopolymers >15 bp
   - dinucleotide repeats >15 bp (>=8 repeats)
   - trinucleotide repeats >15 bp (>=6 repeats; excluding AAA/TTT/GGG/CCC)
2) Writes a log with per-fragment reasons and counts:
   - `<prefix>.lowcomplexity.log`
3) Clusters passing fragments at 85% identity (centroids only):
   - `vsearch --cluster_fast ... --id 0.85 --strand both`

Outputs:
- `<prefix>.fragments.lc.fa`
- `<prefix>.fragments.nr85.fa`

Run:
```bash
bash 02_filter_lc_cluster_nr85.sh SINEBase.nr95.fragments.norm.fa 8
```

---

## Typical usage (end-to-end)

```bash
bash 01_sanitize_shred.sh SINEBase.nr95.fa SINEBase.nr95
bash 02_filter_lc_cluster_nr85.sh SINEBase.nr95.fragments.norm.fa 8
```

Your final query set for genome searching:
- `SINEBase.nr95.fragments.nr85.fa`

---

## Notes / rationale

- **Uniqueness**: IDs are forced unique early (`_<record_number>`) to avoid collisions during downstream parsing.
- **Sanitization**: only the ID token is sanitized (replacing `|` and `:`), because those are common troublemakers
  (e.g., EMBOSS `cons` failing on `|`, and `:` colliding with coordinate parsing). Descriptions after whitespace
  are dropped (current behavior), consistent with many FASTA-processing tools.
- **Fragment header coordinates** are designed to support later clustering/diagnostics and to make it obvious which
  region of the original SINE produced the fragment.

---

## Output files (summary)

For input prefix `SINEBase.nr95`:

- `SINEBase.nr95.fragments.norm.fa`  
  All fragments (normalized coordinates in headers)

- `SINEBase.nr95.fragments.lc.fa`  
  Low-complexity filtered fragments

- `SINEBase.nr95.fragments.nr85.fa`  
  Nonredundant fragments (centroids at 85% identity)

- `SINEBase.nr95.lowcomplexity.log`  
  Filtering report (IDs + reasons + counts)

---

## Next step: genome search

These scripts only build the query set. Use your genome scanning workflow
(e.g., `ssearch36` chunk scanning + merge + BIGFLANK extraction) with:

`<prefix>.fragments.nr85.fa`
