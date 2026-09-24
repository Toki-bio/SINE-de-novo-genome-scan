#!/usr/bin/env python3
"""Filter merged scan loci by MAFFT + esl-alistat orientation signal."""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

from orient_lib import (
    build_query_alias,
    classify_orient,
    genome_slice,
    orient_metrics,
    read_fa,
    sanitize_query_id,
)


def load_bed6(path: Path) -> list[list[str]]:
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        rows.append(f)
    return rows


def overlap(hit: list[str], locus: list[str]) -> bool:
    if hit[0] != locus[0]:
        return False
    return int(hit[1]) < int(locus[2]) and int(locus[1]) < int(hit[2])


def pick_best_hit(locus: list[str], hits: list[list[str]]) -> list[str] | None:
    cand = [h for h in hits if overlap(h, locus)]
    if not cand:
        return None
    return max(cand, key=lambda h: float(h[4]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--hits", type=Path, help="all_hits.to_genome.bed (6 col)")
    ap.add_argument("--loci", type=Path, required=True, help="merged_loci.bed")
    ap.add_argument("--out-ok", type=Path, required=True)
    ap.add_argument("--out-fail", type=Path, required=True)
    ap.add_argument("--tsv", type=Path, required=True)
    ap.add_argument("--min-true-id", type=float, default=65.0)
    ap.add_argument("--delta-reject", type=float, default=25.0)
    ap.add_argument("--delta-pass", type=float, default=30.0)
    ap.add_argument("--keep-weak", action="store_true", help="write WEAK loci to --out-ok")
    ap.add_argument("--tmpdir", type=Path, default=None)
    args = ap.parse_args()

    if not args.hits:
        print("ERROR: --hits required", file=sys.stderr)
        return 2

    qfa = read_fa(args.queries)
    alias = build_query_alias(qfa)
    genome = read_fa(args.genome)
    hits = load_bed6(args.hits)
    loci = load_bed6(args.loci)

    ok_lines: list[str] = []
    fail_lines: list[str] = []
    tsv_rows: list[str] = [
        "\t".join(
            [
                "chrom",
                "start",
                "end",
                "query",
                "strand",
                "scan_pident",
                "esl_true",
                "esl_rev",
                "esl_revcomp",
                "delta",
                "verdict",
            ]
        )
    ]

    tmp_parent = args.tmpdir or args.out_ok.parent
    tmp_parent.mkdir(parents=True, exist_ok=True)

    n_err = 0
    with tempfile.TemporaryDirectory(dir=tmp_parent) as td:
        work = Path(td)
        for loc in loci:
            best = pick_best_hit(loc, hits)
            if best is None:
                n_err += 1
                fail_lines.append("\t".join(loc))
                tsv_rows.append(
                    "\t".join(
                        [
                            loc[0],
                            loc[1],
                            loc[2],
                            ".",
                            loc[5] if len(loc) > 5 else ".",
                            loc[4] if len(loc) > 4 else ".",
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "REJECT_NO_HIT",
                        ]
                    )
                )
                continue

            qkey = sanitize_query_id(best[3])
            qname = alias.get(qkey)
            if not qname or qname not in qfa:
                n_err += 1
                fail_lines.append("\t".join(loc))
                tsv_rows.append(
                    "\t".join(
                        [
                            loc[0],
                            loc[1],
                            loc[2],
                            best[3],
                            best[5],
                            best[4],
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "REJECT_QUERY",
                        ]
                    )
                )
                continue

            qseq = qfa[qname]
            sseq = genome_slice(
                genome,
                best[0],
                int(best[1]),
                int(best[2]),
                best[5],
            )
            if not qseq or not sseq:
                n_err += 1
                fail_lines.append("\t".join(loc))
                tsv_rows.append(
                    "\t".join(
                        [
                            loc[0],
                            loc[1],
                            loc[2],
                            qname,
                            best[5],
                            best[4],
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "REJECT_SEQ",
                        ]
                    )
                )
                continue

            m = orient_metrics(qseq, sseq, work)
            if m is None:
                n_err += 1
                fail_lines.append("\t".join(loc))
                tsv_rows.append(
                    "\t".join(
                        [
                            loc[0],
                            loc[1],
                            loc[2],
                            qname,
                            best[5],
                            best[4],
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "REJECT_ALIGN",
                        ]
                    )
                )
                continue

            verdict = classify_orient(
                m,
                min_true_id=args.min_true_id,
                delta_reject=args.delta_reject,
                delta_pass=args.delta_pass,
            )
            tsv_rows.append(
                "\t".join(
                    [
                        loc[0],
                        loc[1],
                        loc[2],
                        qname,
                        best[5],
                        best[4],
                        f"{m['esl_true']:.2f}",
                        f"{m['esl_rev']:.2f}",
                        f"{m['esl_revcomp']:.2f}",
                        f"{m['delta']:.2f}",
                        verdict,
                    ]
                )
            )
            if verdict == "PASS" or (verdict == "WEAK" and args.keep_weak):
                ok_lines.append("\t".join(loc))
            else:
                fail_lines.append("\t".join(loc))

    args.out_ok.write_text("\n".join(ok_lines) + ("\n" if ok_lines else ""))
    args.out_fail.write_text("\n".join(fail_lines) + ("\n" if fail_lines else ""))
    args.tsv.write_text("\n".join(tsv_rows) + "\n")

    print(
        f"orient_filter: loci={len(loci)} pass={len(ok_lines)} fail={len(fail_lines)} errors={n_err}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
