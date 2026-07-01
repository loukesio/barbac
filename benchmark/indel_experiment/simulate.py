"""Simulate barcode reads with substitution + insertion + deletion errors.

Ground-truth design: a fixed library of N distinct random barcodes of length L,
each assigned a log-normal abundance that sums to the total read count. Each
read is one mutated copy of its parent barcode under a per-base error model.

Outputs three files in the target directory:
    true_counts.csv  - barcode, true_count    (the planned, error-free counts)
    input.csv        - barcode, counts        (collapsed observed reads)
    shepherd_input.txt - same as input.csv but TSV no-header, for shepherd_t0.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np

ALPHABET = np.array(list("ACGT"))


def random_barcodes(n: int, length: int, rng: np.random.Generator) -> list[str]:
    idx = rng.integers(0, 4, size=(n, length))
    return ["".join(ALPHABET[row]) for row in idx]


def lognormal_abundances(n: int, total_reads: int, sigma: float,
                         rng: np.random.Generator) -> np.ndarray:
    raw = rng.lognormal(mean=0.0, sigma=sigma, size=n)
    raw = raw / raw.sum() * total_reads
    counts = np.maximum(1, np.round(raw)).astype(np.int64)
    # Adjust the largest cluster so counts.sum() == total_reads exactly.
    counts[counts.argmax()] += total_reads - counts.sum()
    assert counts.sum() == total_reads and counts.min() >= 1
    return counts


def _mutate_one_barcode(bc: str, n_reads: int, sub_rate: float,
                        ins_rate: float, del_rate: float,
                        rng: np.random.Generator) -> Counter:
    """Generate n_reads noisy copies of bc and return their collapsed counts."""
    L = len(bc)
    if n_reads == 0:
        return Counter()

    # Per-read, per-base random decisions.
    subs = rng.random((n_reads, L)) < sub_rate
    dels = rng.random((n_reads, L)) < del_rate
    inss = rng.random((n_reads, L + 1)) < ins_rate

    # Random base choices for substitutions (offset so we never pick the
    # original) and for insertions (uniform).
    sub_offset = rng.integers(1, 4, size=(n_reads, L))
    ins_base   = rng.integers(0, 4, size=(n_reads, L + 1))

    bc_idx = np.array(["ACGT".index(c) for c in bc])  # original base index per position

    out = Counter()
    for r in range(n_reads):
        parts = []
        for j in range(L):
            if inss[r, j]:
                parts.append("ACGT"[ins_base[r, j]])
            if dels[r, j]:
                continue
            if subs[r, j]:
                parts.append("ACGT"[(bc_idx[j] + sub_offset[r, j]) % 4])
            else:
                parts.append(bc[j])
        if inss[r, L]:
            parts.append("ACGT"[ins_base[r, L]])
        out["".join(parts)] += 1
    return out


def simulate(
    out_dir: str | Path,
    *,
    n_barcodes: int,
    barcode_length: int,
    n_reads: int,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    sigma: float = 1.5,
    seed: int = 0,
) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bc_rng    = np.random.default_rng(seed)
    abund_rng = np.random.default_rng(seed + 1)
    mut_rng   = np.random.default_rng(seed + 2)

    barcodes = random_barcodes(n_barcodes, barcode_length, bc_rng)
    true_counts = lognormal_abundances(n_barcodes, n_reads, sigma, abund_rng)

    # Sort descending by true count so the file matches the original benchmark layout.
    order = np.argsort(-true_counts)
    barcodes    = [barcodes[i] for i in order]
    true_counts = true_counts[order]

    with (out_dir / "true_counts.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["BC", "True Count"])
        for bc, c in zip(barcodes, true_counts):
            w.writerow([bc, int(c)])

    # Generate noisy reads, collapsing as we go.
    all_reads: Counter = Counter()
    for bc, c in zip(barcodes, true_counts):
        all_reads.update(_mutate_one_barcode(
            bc, int(c), sub_rate, ins_rate, del_rate, mut_rng,
        ))

    # Write input.csv (sorted desc by count) and shepherd_input.txt (same data, TSV).
    items = sorted(all_reads.items(), key=lambda kv: -kv[1])
    with (out_dir / "input.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode", "counts"])
        for seq, c in items:
            w.writerow([seq, c])
    with (out_dir / "shepherd_input.txt").open("w") as fh:
        for seq, c in items:
            fh.write(f"{seq}\t{c}\n")

    # Manifest so the runner can re-load the config later.
    manifest = {
        "n_barcodes":     n_barcodes,
        "barcode_length": barcode_length,
        "n_reads":        n_reads,
        "sub_rate":       sub_rate,
        "ins_rate":       ins_rate,
        "del_rate":       del_rate,
        "sigma":          sigma,
        "seed":           seed,
        "n_unique_reads": len(items),
    }
    import json
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return out_dir


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("out_dir")
    p.add_argument("--n-barcodes",     type=int,   default=10_000)
    p.add_argument("--barcode-length", type=int,   default=20)
    p.add_argument("--n-reads",        type=int,   default=1_000_000)
    p.add_argument("--sub-rate",       type=float, default=0.005)
    p.add_argument("--ins-rate",       type=float, default=0.0)
    p.add_argument("--del-rate",       type=float, default=0.0)
    p.add_argument("--sigma",          type=float, default=1.5)
    p.add_argument("--seed",           type=int,   default=0)
    args = p.parse_args()
    simulate(
        args.out_dir,
        n_barcodes=args.n_barcodes,
        barcode_length=args.barcode_length,
        n_reads=args.n_reads,
        sub_rate=args.sub_rate,
        ins_rate=args.ins_rate,
        del_rate=args.del_rate,
        sigma=args.sigma,
        seed=args.seed,
    )
