SINE Fragment Scanner

Fragment-based genome scanning pipeline for discovery of SINE elements

This repository contains a two-stage, fully reproducible pipeline:

Build a high-quality fragment query set from a SINE database

Scan genomes with those fragments using sensitive Smith–Waterman searches

All files produced by the pipeline can be traced back to a specific script and command.

Overview of the workflow
SINE database FASTA
        │
        ▼
01_sanitize_shred.sh
        │
        └──► *.fragments.norm.fa
                │
                ▼
02_filter_lc_cluster_nr85.sh
                │
                ├──► *.fragments.lc.fa
                └──► *.fragments.nr85.fa   ← FINAL QUERY SET
                                │
                                ▼
sine_scan.sh
                                │
                                ├──► per-query search logs
                                ├──► merged hits
                                ├──► clustered loci (BED)
                                └──► extracted candidate SINE sequences (FASTA)

Stage 1 — Build fragment query set

This stage turns a heterogeneous SINE database into fixed-length, non-redundant, low-complexity-filtered fragments suitable for large-scale genome searches.

Script 1: 01_sanitize_shred.sh

Purpose:
Normalize SINE sequences, enforce safe/unique FASTA IDs, and shred each SINE into biologically meaningful fragments with explicit coordinates.

Input:

A SINE database FASTA (e.g. SINEBase.nr95.fa)

Command:

bash 01_sanitize_shred.sh SINEBase.nr95.fa SINEBase.nr95


What the script actually does (from code):

Sanitizes FASTA IDs

Keeps only the first whitespace-delimited token

Replaces | and : with _ (required for EMBOSS / downstream parsers)

Appends _<record_number> to guarantee uniqueness

Normalizes sequences

Uppercase

Single-line FASTA (seqkit seq -u -w 0)

Splits each SINE into regions

5′ region: bases 1–150

Middle region: bases 151–(L−100) if L > 250

3′ region: last 99 bp if L > 150

Sliding-window fragmentation

Fragment length: 50 bp

Step size:

5′ region: 10 bp

middle / 3′: 25 bp

Normalizes fragment headers
Each fragment header encodes:

source SINE ID

region (5p, mid, 3p)

true coordinates within the original SINE

Header format:

>SINE_ID|REGION:START-END


Output (important):

SINEBase.nr95.fragments.norm.fa


This file contains all fragments, including low-complexity and redundant ones.
It is not intended for genome searching yet.

Script 2: 02_filter_lc_cluster_nr85.sh

Purpose:
Remove obvious low-complexity fragments and collapse redundancy across SINE families.

Input:

<prefix>.fragments.norm.fa

Command:

bash 02_filter_lc_cluster_nr85.sh SINEBase.nr95.fragments.norm.fa 8


(8 = number of threads for vsearch)

What the script actually does:

1. Low-complexity filtering (explicit rules)

Fragments are removed if they contain:

Homopolymers >15 bp

Dinucleotide repeats >15 bp (≥8 repeats)

Trinucleotide repeats >15 bp (≥6 repeats), excluding AAA/TTT/CCC/GGG

A full log is written with fragment IDs and reasons.

2. Non-redundant clustering

Remaining fragments are clustered using vsearch:

Identity threshold: 85%

Strand: both

Output: centroids only

This collapses:

overlapping fragments

near-identical fragments from different SINE families

redundant windows within the same SINE

Outputs (this answers your key question):

SINEBase.nr95.fragments.lc.fa     # low-complexity filtered
SINEBase.nr95.fragments.nr85.fa   # FINAL QUERY SET
SINEBase.nr95.lowcomplexity.log


👉 SINEBase.fragments.nr85.fa is produced by this command and only this command.

This is the file you use for genome scanning.

Stage 2 — Genome scanning
Script 3: sine_scan.sh

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
Design principles (as implemented)

Fragment-level sensitivity: detect highly diverged SINEs

Safe FASTA headers: compatible with EMBOSS, BEDTools, samtools

Reproducible outputs: every file has a single, traceable origin

No hidden heuristics: all thresholds are explicit in scripts

Separation of concerns: query construction ≠ genome scanning

Typical end-to-end run
# Build fragment queries
bash 01_sanitize_shred.sh SINEBase.nr95.fa SINEBase.nr95
bash 02_filter_lc_cluster_nr85.sh SINEBase.nr95.fragments.norm.fa 8

# Genome scan
bash sine_scan.sh \
  -q SINEBase.nr95.fragments.nr85.fa \
  -g genome.fa \
  -o sine_search_out/genome
