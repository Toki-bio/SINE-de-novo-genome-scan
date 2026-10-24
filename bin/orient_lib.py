"""MAFFT pair alignment + esl-alistat orient metrics. See ORIENT_NOISE_FLOOR.md."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")

MAFFT = [
    "mafft",
    "--localpair",
    "--maxiterate",
    "1000",
    "--ep",
    "0.123",
    "--nuc",
    "--quiet",
]


def read_fa(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name, buf = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def sanitize_query_id(name: str) -> str:
    return name.replace("|", "_").replace(":", "_")


def build_query_alias(qfa: dict[str, str]) -> dict[str, str]:
    return {sanitize_query_id(k): k for k in qfa}


def genome_slice(genome: dict[str, str], chrom: str, start: int, end: int, strand: str) -> str:
    if chrom not in genome:
        for k in genome:
            if k.split()[0] == chrom:
                chrom = k
                break
    seq = genome[chrom][start:end]
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
    return seq


def parse_aligned_fasta(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    name, parts = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(parts)
            name = line[1:].split()[0]
            parts = []
        elif line.strip():
            parts.append(line.strip())
    if name:
        out[name] = "".join(parts).upper()
    return out


def mafft_pair(qseq: str, sseq: str, work: Path) -> tuple[str, str] | None:
    inf = work / "pair.fa"
    inf.write_text(f">q\n{qseq}\n>s\n{sseq}\n")
    proc = subprocess.run(
        MAFFT + [str(inf)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return None
    aln = parse_aligned_fasta(proc.stdout + proc.stderr)
    if "q" not in aln or "s" not in aln:
        return None
    q_aln, s_aln = aln["q"], aln["s"]
    if len(q_aln) != len(s_aln):
        return None
    return q_aln, s_aln


def esl_avg_identity(q_aln: str, s_aln: str, work: Path) -> float | None:
    af = work / "aln.fa"
    af.write_text(f">q\n{q_aln}\n>s\n{s_aln}\n")
    proc = subprocess.run(
        ["esl-alistat", "-1", "--dna", str(af)],
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if parts[0].isdigit():
            return float(parts[-1])
    return None


def ungapped_bases(gapped: str) -> str:
    return "".join(c for c in gapped if c in "ACGT")


def refill_gapped_row(row: str, new_bases: str) -> str:
    out, i = [], 0
    for c in row:
        if c in "ACGT":
            out.append(new_bases[i])
            i += 1
        else:
            out.append(c)
    return "".join(out)


def orient_null_rows(q_aln: str, s_aln: str) -> tuple[str, str]:
    sb = ungapped_bases(s_aln)
    rev_s = refill_gapped_row(s_aln, sb[::-1])
    rvc_s = refill_gapped_row(s_aln, sb.translate(COMP)[::-1])
    return rev_s, rvc_s


def orient_metrics(
    qseq: str,
    sseq: str,
    work: Path,
) -> dict[str, float] | None:
    pair = mafft_pair(qseq, sseq, work)
    if pair is None:
        return None
    q_aln, s_aln = pair
    true_id = esl_avg_identity(q_aln, s_aln, work)
    if true_id is None:
        return None
    rev_s, rvc_s = orient_null_rows(q_aln, s_aln)
    rev_id = esl_avg_identity(q_aln, rev_s, work)
    rvc_id = esl_avg_identity(q_aln, rvc_s, work)
    if rev_id is None or rvc_id is None:
        return None
    null = max(rev_id, rvc_id)
    return {
        "esl_true": true_id,
        "esl_rev": rev_id,
        "esl_revcomp": rvc_id,
        "esl_null": null,
        "delta": true_id - null,
    }


def classify_orient(
    metrics: dict[str, float],
    *,
    min_true_id: float,
    delta_reject: float,
    delta_pass: float,
) -> str:
    if metrics["esl_true"] < min_true_id:
        return "REJECT_LOW_ID"
    if metrics["delta"] < delta_reject:
        return "REJECT_ORIENT"
    if metrics["delta"] < delta_pass:
        return "WEAK"
    return "PASS"
