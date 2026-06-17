#!/usr/bin/env python3
"""
Parse Google Benchmark JSON from bench_bls12_381 and bench_blst_complete.

- Single file: amortized ns for batch mod-mul / reduce-expand; grouped table.
- With --blst: align RNS vs BLST via RNS_BLST_MAP and report speedup ratio (blst / rns).
  For batch benchmarks, uses amortized ns per operation on the RNS side.

Unmapped entries are listed at the end for manual pairing.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

BATCH_MODMUL_RE = re.compile(r"^BM_BatchModMul_(\d+)$")
BATCH_REDUCE_RE = re.compile(r"^BM_BatchReduceExpand_(\d+)$")
BENCH_BATCH_MODMUL_RE = re.compile(r"^BM_BatchModMul_(?:Matrix|MatrixNoK)_(\d+)$")
BENCH_BATCH_REDUCE_RE = re.compile(r"^BM_BatchReduceExpand_(?:Matrix|MatrixNoK)_(\d+)$")

# RNS (bench_bls12_381 Matrix rows) name -> BLST (bench_blst_complete) name
# bench_bls12_381 registers BM_<Op>_Matrix; canonical keys omit the suffix.
RNS_BLST_MAP: dict[str, str] = {
    "BM_ModMul": "BM_Fp_Mul",
    "BM_FP2_Mul": "BM_Fp2_Mul",
    "BM_FP2_Sqr": "BM_Fp2_Sqr",
    "BM_G1_Add": "BM_G1_Add",
    "BM_G1_Double": "BM_G1_Double",
    "BM_G1_MixedAdd": "BM_G1_MixedAdd",
    "BM_G1_ScalarMult_GLV": "BM_G1_ScalarMul",
    "BM_G2_Add": "BM_G2_Add",
    "BM_G2_Double": "BM_G2_Double",
    "BM_G2_MixedAdd": "BM_G2_MixedAdd",
    "BM_G2_ScalarMult_GLS": "BM_G2_ScalarMul",
    "BM_FP12_Mul": "BM_Fp12_Mul",
    "BM_FP12_Sqr": "BM_Fp12_Sqr",
    "BM_FP12_CyclotomicSqr": "BM_Fp12_CyclotomicSqr",
    "BM_FP12_Conjugate": "BM_Fp12_Conjugate",
    "BM_FP12_Inverse_BLST": "BM_Fp12_Inverse",
    "BM_FP12_Inverse_RNS_BLST": "BM_Fp12_Inverse",
    "BM_MillerLoop": "BM_MillerLoop",
    # Computed in enrich_rns_rows: merged RNS step minus fp12 sparse mult
    "BM_MillerLoop_LineDbl": "BM_LineDbl",
    "BM_MillerLoop_LineAdd": "BM_LineAdd",
    "BM_FP12_FinalExp_BLST_Inverter": "BM_FinalExp",
    "BM_Pairing_RNS_BLST_Inverter": "BM_FullPairing",
}

# BatchModMul_N and BatchModMul_Matrix_N compare to field mul baseline
BATCH_MODMUL_BLST = "BM_Fp_Mul"

# Omit from all parser output
EXCLUDE_FROM_OUTPUT: frozenset[str] = frozenset({
    "BM_FP12_FrobeniusMap1",
    "BM_FP12_FrobeniusMap2",
    "BM_FP12_FrobeniusMap3",
})

# Pairing verify e()=e(): not benchmarked directly in RNS; sum component costs.
# Two Miller loops, one FP12 mul to combine partial products, one final exponentiation.
PAIRING_VERIFY_COMPUTED = "BM_Pairing_ComputedVerify"
PAIRING_VERIFY_TERMS: tuple[tuple[str, int], ...] = (
    ("BM_MillerLoop", 2),
    ("BM_FP12_Mul", 1),
    ("BM_FP12_FinalExp_BLST_Inverter", 1),
)
PAIRING_VERIFY_BLST_TERMS: tuple[tuple[str, int], ...] = (
    ("BM_MillerLoop", 2),
    ("BM_Fp12_Mul", 1),
    ("BM_FinalExp", 1),
)

# Raw RNS benchmarks superseded by computed rows or not used for BLST compare
EXCLUDE_FROM_COMPARE: frozenset[str] = frozenset({
    "BM_MillerLoop_DoubleStep",
    "BM_MillerLoop_AddStep",
    "BM_FP12_MulByXy0z00",
    "BM_FP12_FinalExpAfterInverse",
    "BM_FP12_FinalExp_RNS_BLST_Inverter",
    "BM_Pairing_BLST_Inverter",
    PAIRING_VERIFY_COMPUTED,
}) | EXCLUDE_FROM_OUTPUT

# Hide merged raw steps from tables (computed rows replace them for reporting)
EXCLUDE_FROM_TABLES: frozenset[str] = EXCLUDE_FROM_COMPARE - {
    "BM_FP12_MulByXy0z00",  # show sparse mult in table
}

# Miller line: RNS merges line eval + sparse multiply; subtract sparse mult for BLST alignment
MILLER_LINE_COMPUTED: tuple[tuple[str, str, str], ...] = (
    ("BM_MillerLoop_LineDbl", "BM_MillerLoop_DoubleStep", "BM_FP12_MulByXy0z00"),
    ("BM_MillerLoop_LineAdd", "BM_MillerLoop_AddStep", "BM_FP12_MulByXy0z00"),
)

DISPLAY_LABELS: dict[str, str] = {
    "BM_FP12_MulByXy0z00": "fp12 sparse mult",
    "BM_Pairing_ComputedVerify": "pairing verify e()=e() (computed)",
}


@dataclass
class Row:
    category: str
    name: str
    label: str
    cpu_ns: float
    batch_size: int | None
    amortized_ns: float | None
    iterations: int

    def comparable_ns(self) -> float:
        """Per-operation time used for RNS vs BLST comparison."""
        if self.amortized_ns is not None:
            return self.amortized_ns
        return self.cpu_ns


def bench_name(entry: dict) -> str:
    return entry.get("name", "")


def bench_label(entry: dict) -> str:
    return entry.get("label", "") or bench_name(entry)


def cpu_time_ns(entry: dict) -> float:
    unit = entry.get("time_unit", "ns")
    cpu = float(entry["cpu_time"])
    if unit == "ns":
        return cpu
    if unit == "us":
        return cpu * 1_000
    if unit == "ms":
        return cpu * 1_000_000
    if unit == "s":
        return cpu * 1_000_000_000
    raise ValueError(f"unsupported time_unit: {unit}")


def canonical_rns_name(name: str) -> str:
    """bench_bls12_381 uses BM_<Op>_Matrix; map keys omit the suffix."""
    if name.endswith("_Matrix"):
        return name[: -len("_Matrix")]
    return name


def row_lookup(by_name: dict[str, Row], canonical: str) -> Row | None:
    if canonical in by_name:
        return by_name[canonical]
    matrix = f"{canonical}_Matrix"
    return by_name.get(matrix)


def excluded_name(name: str, excluded: frozenset[str]) -> bool:
    return name in excluded or canonical_rns_name(name) in excluded


def is_matrix_nok(name: str, label: str) -> bool:
    text = f"{name} {label}"
    return "MatrixNoK" in text or "matrix_nok" in text.lower()


def batch_size(name: str) -> int | None:
    for pat in (BATCH_MODMUL_RE, BATCH_REDUCE_RE, BENCH_BATCH_MODMUL_RE, BENCH_BATCH_REDUCE_RE):
        m = pat.match(name)
        if m:
            return int(m.group(1))
    return None


def is_batch_modmul(name: str) -> bool:
    return (
        BATCH_MODMUL_RE.match(name) is not None
        or BENCH_BATCH_MODMUL_RE.match(name) is not None
    )


def categorize(name: str, label: str) -> str:
    n = name
    if "ModMul" in n or "BatchModMul" in n:
        return "modmul"
    if "BatchReduceExpand" in n:
        return "reduce_expand"
    if n.startswith("BM_FP2_") or n.startswith("BM_Fp2_") or "FP2_" in label:
        return "fp2"
    if n.startswith("BM_Fp_") and "Fp2" not in n and "Fp12" not in n:
        return "fp"
    if any(x in n for x in ("G1_", "G2_")) or label.startswith(("G1_", "G2_")):
        return "ec"
    if "MillerLoop" in n or "Pairing" in n or "FullPairing" in n:
        return "pairing"
    if "Line" in n:
        return "miller_step"
    if "FP12" in n or "Fp12" in n or "Inverse" in n or "FinalExp" in n:
        return "fp12_pairing"
    return "other"


def parse_entries(entries: Iterable[dict], include_nok: bool, *, source: str) -> list[Row]:
    rows: list[Row] = []
    for entry in entries:
        name = bench_name(entry)
        label = bench_label(entry)
        if source == "rns" and not include_nok and is_matrix_nok(name, label):
            continue
        cpu = cpu_time_ns(entry)
        bs = batch_size(name)
        amort = (cpu / bs) if bs else None
        rows.append(
            Row(
                category=categorize(name, label),
                name=name,
                label=label,
                cpu_ns=cpu,
                batch_size=bs,
                amortized_ns=amort,
                iterations=int(entry.get("iterations", 0)),
            )
        )
    return rows


def load_json(path: Path) -> list[dict]:
    data = json.loads(path.read_text())
    return data.get("benchmarks", data)


def display_label(row: Row) -> str:
    return DISPLAY_LABELS.get(row.name, row.label)


def apply_display_labels(rows: list[Row]) -> list[Row]:
    out: list[Row] = []
    for r in rows:
        label = DISPLAY_LABELS.get(r.name, r.label)
        if label != r.label:
            out.append(
                Row(
                    category=r.category,
                    name=r.name,
                    label=label,
                    cpu_ns=r.cpu_ns,
                    batch_size=r.batch_size,
                    amortized_ns=r.amortized_ns,
                    iterations=r.iterations,
                )
            )
        else:
            out.append(r)
    return out


def compute_miller_line_rows(by_name: dict[str, Row]) -> list[Row]:
    """RNS Miller steps include sparse multiply; subtract for BLST line alignment."""
    sparse = row_lookup(by_name, "BM_FP12_MulByXy0z00")
    if sparse is None:
        return []
    sparse_ns = sparse.comparable_ns()
    computed: list[Row] = []
    for computed_name, step_name, sparse_name in MILLER_LINE_COMPUTED:
        step = row_lookup(by_name, step_name)
        if step is None or row_lookup(by_name, sparse_name) is None:
            continue
        line_ns = step.comparable_ns() - sparse_ns
        computed.append(
            Row(
                category="miller_step",
                name=computed_name,
                label=f"{computed_name} ({step_name} - fp12 sparse mult)",
                cpu_ns=line_ns,
                batch_size=None,
                amortized_ns=None,
                iterations=step.iterations,
            )
        )
    return computed


def pairing_verify_formula(terms: tuple[tuple[str, int], ...]) -> str:
    return " + ".join(
        f"{count}×{name}" if count > 1 else name for name, count in terms
    )


def compute_pairing_verify_ns(
    by_name: dict[str, Row],
    terms: tuple[tuple[str, int], ...],
) -> float | None:
    total_ns = 0.0
    for name, count in terms:
        row = row_lookup(by_name, name)
        if row is None:
            return None
        total_ns += count * row.comparable_ns()
    return total_ns


def compute_pairing_verify_row(by_name: dict[str, Row]) -> Row | None:
    """Pairing verify e()=e(): 2×MillerLoop + FP12_Mul + FinalExp (not measured in RNS)."""
    total_ns = compute_pairing_verify_ns(by_name, PAIRING_VERIFY_TERMS)
    if total_ns is None:
        return None
    return Row(
        category="pairing",
        name=PAIRING_VERIFY_COMPUTED,
        label=f"pairing verify e()=e(): {pairing_verify_formula(PAIRING_VERIFY_TERMS)}",
        cpu_ns=total_ns,
        batch_size=None,
        amortized_ns=None,
        iterations=0,
    )


def enrich_rns_rows(rows: list[Row]) -> list[Row]:
    rows = [r for r in rows if not excluded_name(r.name, EXCLUDE_FROM_OUTPUT)]
    rows = apply_display_labels(rows)
    by_name = rows_by_name(rows)
    rows.extend(compute_miller_line_rows(by_name))
    verify = compute_pairing_verify_row(by_name)
    if verify is not None:
        rows.append(verify)
    return rows


def rows_by_name(rows: list[Row]) -> dict[str, Row]:
    return {r.name: r for r in rows}


def blst_target_for_rns(rns_name: str) -> str | None:
    key = canonical_rns_name(rns_name)
    if key in RNS_BLST_MAP:
        return RNS_BLST_MAP[key]
    if is_batch_modmul(rns_name):
        return BATCH_MODMUL_BLST
    return None


def format_rns_table(rows: list[Row]) -> str:
    order = ("modmul", "reduce_expand", "fp", "fp2", "ec", "fp12_pairing", "miller_step", "pairing", "other")
    by_cat: dict[str, list[Row]] = {c: [] for c in order}
    for r in rows:
        if excluded_name(r.name, EXCLUDE_FROM_TABLES):
            continue
        by_cat.setdefault(r.category, []).append(r)

    lines: list[str] = []
    headers = f"{'benchmark':<42} {'cpu_ns':>12} {'batch':>6} {'amort_ns':>12}  label"
    sep = "-" * len(headers)

    for cat in order:
        cat_rows = by_cat.get(cat, [])
        if not cat_rows:
            continue
        lines.append(f"\n## {cat}")
        lines.append(headers)
        lines.append(sep)
        for r in sorted(cat_rows, key=lambda x: x.name):
            batch_s = str(r.batch_size) if r.batch_size else "-"
            amort_s = f"{r.amortized_ns:12.1f}" if r.amortized_ns is not None else f"{'':>12}"
            lines.append(
                f"{r.name:<42} {r.cpu_ns:12.1f} {batch_s:>6} {amort_s}  {r.label}"
            )
    return "\n".join(lines).lstrip() + "\n"


def format_blst_table(rows: list[Row]) -> str:
    lines = ["## blst", f"{'benchmark':<32} {'cpu_ns':>12}  label", "-" * 48]
    for r in sorted(rows, key=lambda x: x.name):
        lines.append(f"{r.name:<32} {r.cpu_ns:12.1f}  {r.label}")
    return "\n".join(lines) + "\n"


def format_pairing_verify(
    rns_rows: list[Row],
    blst_rows: list[Row] | None = None,
) -> str:
    by_name = rows_by_name(rns_rows)
    computed = by_name.get(PAIRING_VERIFY_COMPUTED)
    if computed is None:
        return ""

    computed_ns = computed.comparable_ns()
    lines = [
        "## Pairing verify e()=e() (not measured in RNS; summed from components)",
        "Two Miller loops, one FP12 multiply to combine, one final exponentiation.",
        f"{'row':<42} {'cpu_ns':>12}",
        "-" * 56,
        f"{PAIRING_VERIFY_COMPUTED:<42} {computed_ns:12.1f}",
        f"  formula: {pairing_verify_formula(PAIRING_VERIFY_TERMS)}",
    ]

    if blst_rows is not None:
        blst_ns = compute_pairing_verify_ns(rows_by_name(blst_rows), PAIRING_VERIFY_BLST_TERMS)
        if blst_ns is not None:
            ratio = blst_ns / computed_ns if computed_ns else float("nan")
            lines.extend([
                f"{'BLST (same formula)':<42} {blst_ns:12.1f}",
                f"{'ratio blst/rns':<42} {ratio:12.3f}  (>1 means RNS is faster)",
                f"  blst formula: {pairing_verify_formula(PAIRING_VERIFY_BLST_TERMS)}",
            ])

    return "\n".join(lines) + "\n"


@dataclass
class CompareRow:
    rns_name: str
    blst_name: str
    rns_ns: float
    blst_ns: float
    ratio: float
    note: str


def build_comparisons(rns_rows: list[Row], blst_rows: list[Row]) -> tuple[list[CompareRow], list[Row], list[Row]]:
    rns_by = rows_by_name(rns_rows)
    blst_by = rows_by_name(blst_rows)
    compared: list[CompareRow] = []
    mapped_rns: set[str] = set()
    mapped_blst: set[str] = set()

    for rns_name in sorted(rns_by):
        if excluded_name(rns_name, EXCLUDE_FROM_COMPARE):
            continue
        blst_name = blst_target_for_rns(rns_name)
        if blst_name is None:
            continue
        if blst_name not in blst_by:
            continue
        rns = rns_by[rns_name]
        blst = blst_by[blst_name]
        rns_ns = rns.comparable_ns()
        blst_ns = blst.comparable_ns()
        note = ""
        if rns.amortized_ns is not None:
            note = f"rns amortized (batch={rns.batch_size})"
        if rns_name in ("BM_MillerLoop_LineDbl", "BM_MillerLoop_LineAdd"):
            note = "computed (merged step - fp12 sparse mult)"
        compared.append(
            CompareRow(
                rns_name=rns_name,
                blst_name=blst_name,
                rns_ns=rns_ns,
                blst_ns=blst_ns,
                ratio=blst_ns / rns_ns if rns_ns else float("nan"),
                note=note,
            )
        )
        mapped_rns.add(rns_name)
        mapped_blst.add(blst_name)

    unmapped_rns = [
        r for r in rns_rows
        if r.name not in mapped_rns and not excluded_name(r.name, EXCLUDE_FROM_COMPARE)
    ]
    unmapped_blst = [r for r in blst_rows if r.name not in mapped_blst]
    return compared, unmapped_rns, unmapped_blst


def format_comparison(compared: list[CompareRow], unmapped_rns: list[Row], unmapped_blst: list[Row]) -> str:
    lines = [
        "## RNS vs BLST (ratio = blst_ns / rns_ns; >1 means RNS is faster)",
        f"{'rns_benchmark':<36} {'blst_benchmark':<24} {'rns_ns':>10} {'blst_ns':>10} {'ratio':>8}  note",
        "-" * 100,
    ]
    for c in compared:
        lines.append(
            f"{c.rns_name:<36} {c.blst_name:<24} {c.rns_ns:10.1f} {c.blst_ns:10.1f} {c.ratio:8.3f}  {c.note}"
        )

    if unmapped_rns:
        lines.append("\n## Unmapped RNS (no BLST pair in RNS_BLST_MAP)")
        for r in sorted(unmapped_rns, key=lambda x: x.name):
            extra = f" amort={r.amortized_ns:.1f}" if r.amortized_ns else ""
            lines.append(f"  {r.name:<42} {r.cpu_ns:10.1f} ns{extra}  {r.label}")

    if unmapped_blst:
        lines.append("\n## Unmapped BLST (not referenced by any RNS mapping)")
        for r in sorted(unmapped_blst, key=lambda x: x.name):
            lines.append(f"  {r.name:<32} {r.cpu_ns:10.1f} ns  {r.label}")

    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("rns_json", type=Path, help="RNS benchmark JSON (bench_bls12_381)")
    ap.add_argument("--blst", type=Path, help="BLST benchmark JSON (bench_blst_complete)")
    ap.add_argument("-o", "--output", type=Path, help="Write RNS table here")
    ap.add_argument(
        "--compare-output",
        type=Path,
        help="Write RNS vs BLST comparison here (requires --blst)",
    )
    ap.add_argument(
        "--include-nok",
        action="store_true",
        help="Include MatrixNoK RNS rows (excluded by default)",
    )
    args = ap.parse_args()

    rns_rows = enrich_rns_rows(
        parse_entries(load_json(args.rns_json), args.include_nok, source="rns")
    )
    if not rns_rows:
        print("No RNS benchmark rows matched.", file=sys.stderr)
        return 1

    parts = [format_rns_table(rns_rows)]

    if args.blst:
        blst_rows = parse_entries(load_json(args.blst), include_nok=True, source="blst")
        if not blst_rows:
            print("No BLST benchmark rows matched.", file=sys.stderr)
            return 1
        verify = format_pairing_verify(rns_rows, blst_rows)
        if verify:
            parts.append(verify)
        parts.append(format_blst_table(blst_rows))
        compared, unmapped_rns, unmapped_blst = build_comparisons(rns_rows, blst_rows)
        comparison = format_comparison(compared, unmapped_rns, unmapped_blst)
        parts.append(comparison)
        if args.compare_output:
            args.compare_output.write_text(comparison)
    else:
        verify = format_pairing_verify(rns_rows)
        if verify:
            parts.append(verify)

    report = "\n".join(parts)
    if args.output:
        args.output.write_text(report)
    print(report, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
