# SINE-de-novo-genome-scan
Perform a de novo genome search for SINEs using a database of known SINE sequences.
Finding a SINE in a genome can be done with searching for even a faint similarity to already known SINE sequnces, because many of them share at least a limited similarity.

This workflow does the following:
1) fragmenting a SINE database of consensus sequences from various genomes into short windows (default 50 bp, 25 bp step),
2) searching each fragment against the studied genome with `ssearch36`,
3) merging hits into anchor loci,
4) expanding anchors to full candidate sequences,
5) prepares sequences for downstream analysis using SubFam or targeted subfamily analysis.

## Requirements
- samtools
- bedtools
- FASTA36 (ssearch36)
- awk, sort, wc, bash

Initial SINE database can be obtained from various sources
SINEbase https://sines.eimb.ru/
Manuscripts  on SINE analysis with consensus sequences
Repeats annotation databases filtered for consistency and verified 

### DATABASE PREPARATION
Input SINEs.fa file is deduplicated using vsearch with 95% similarity cutoff 

### Script 1: 01_sanitize_shred.sh

Purpose:
Normalize SINE sequences, enforce safe/unique FASTA IDs, and shred each SINE into biologically meaningful fragments with explicit coordinates.

Input:

A SINE database FASTA (e.g. SINEBase.nr95.fa)

Command:

bash 01_sanitize_shred.sh SINEBase.nr95.fa SINEBase.nr95


What the script actually does (from code):

Sanitizes FASTA IDs, keeps only the first whitespace-delimited token, replaces | and : with _ (required for EMBOSS / downstream parsers), appends _<record_number> to guarantee uniqueness, uppercase, single-line FASTA.

Splits each SINE of length L into regions:

 5′ region: bases 1–150 with 10bp step
 middle region: bases 151–(L−100) if L>250 with 25bp step
 3′ region: last 99 bp if L>150  with 25bp step

Header format:
>SINE_ID|REGION:START-END
source SINE ID
region (5p, mid, 3p)
true coordinates within the original SINE

Output:
SINEBase.nr95.fragments.norm.fa

This file contains all fragments, including low-complexity and redundant ones. It is not intended for genome searching yet (needs deduplication).

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

Stage 2 — Genome scanning
### Script 3: sine_scan.sh

Purpose:
Search a genome for SINE-like loci using the fragment query set.

Input:
Genome FASTA

Fragment query FASTA (*.fragments.nr85.fa)

Command (example):

bash sine_scan.sh \
  -q SINEBase.nr95.fragments.nr85.fa \
  -g genome.fa \
  -o sine_search_out/genome_name


What the script does:

Optional genome subsampling

By fraction (e.g. 5% of genome)

Or by fixed bp count

Uses random scaffold selection (not positional bias)

Per-fragment Smith–Waterman search

Tool: ssearch36

Parallelized at the query level

Each fragment is searched independently

Live progress reporting:

[QUERY 412/926] src=5S-Sauria region=5p 11–60 raw=42 filtered=3


Filtering

Minimum identity

Minimum query coverage

Hit merging

All filtered hits merged into a single file

Converted to BED

Candidate locus clustering

Genomic hits merged within 500 bp

Produces candidate SINE loci

Sequence extraction

Flanks added

FASTA extracted for downstream validation

Key outputs:

query_summary.tsv          # per-fragment statistics
all_hits.filtered.m8       # merged filtered hits
candidate_loci.bed         # clustered loci
candidates.fa              # extracted sequences

File provenance summary (important)
File	Produced by	Command
*.fragments.norm.fa	01_sanitize_shred.sh	bash 01_sanitize_shred.sh …
*.fragments.lc.fa	02_filter_lc_cluster_nr85.sh	same
*.fragments.nr85.fa	02_filter_lc_cluster_nr85.sh	same
query_summary.tsv	sine_scan.sh	genome scan
candidate_loci.bed	sine_scan.sh	genome scan
candidates.fa	sine_scan.sh	genome scan
