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


def structured_barcodes(n: int, template: str,
                        rng: np.random.Generator) -> list[str]:
    """Draw n distinct barcodes from a fixed-anchor design template.

    The template is a string over {A,C,G,T,N}: each 'N' is filled with a
    uniform random base, every other position is a fixed anchor kept constant
    across all barcodes (e.g. 'NNNNNNNNATGCNNNNNNNNATCGTTAA'). Diversity is
    4 ** (number of N positions), so keep enough N's that distinct barcodes
    stay separated at the clustering distance.
    """
    template = template.upper()
    bad = set(template) - set("ACGTN")
    if bad:
        raise ValueError(f"template may only contain A/C/G/T/N, got {sorted(bad)}")
    var_pos = [i for i, ch in enumerate(template) if ch == "N"]
    n_var = len(var_pos)
    if n_var == 0:
        raise ValueError("template has no variable (N) positions")
    max_distinct = 4 ** n_var
    if n > max_distinct:
        raise ValueError(f"template has {n_var} variable positions "
                         f"(<= {max_distinct} distinct barcodes) but {n} requested")
    base = list(template)
    seen: set[str] = set()
    out: list[str] = []
    while len(out) < n:
        draws = rng.integers(0, 4, size=(n - len(out), n_var))
        for row in draws:
            for k, p in enumerate(var_pos):
                base[p] = "ACGT"[row[k]]
            s = "".join(base)
            if s not in seen:
                seen.add(s)
                out.append(s)
                if len(out) == n:
                    break
    return out


def johnson_abundances(n: int, total_reads: int,
                       rng: np.random.Generator) -> np.ndarray:
    """Abundances shaped like the benchmark library of Johnson et al. (2023).

    That library mixes three exponential components -- the bulk at mean 1, a
    hundredth of a percent at mean 10, and a handful at mean 1000 -- which
    gives a much heavier head and a lighter low-count tail than a lognormal of
    the same median. The tail is what decides how many barcodes land in the
    range where no method can recover them, so matching its shape matters when
    the simulation is meant to stand in for a real library.

    Reference: Johnson MS, Venkataram S, Kryazhimskiy S (2023) J Mol Evol
    91:263-280, doi:10.1007/s00239-022-10083-z; data doi:10.5281/zenodo.7052124
    """
    n_high = max(1, round(n * 5e-5))          # 5 in 100,000
    n_mid = max(1, round(n * 1e-3))           # 100 in 100,000
    n_low = n - n_mid - n_high
    if n_low <= 0:
        raise ValueError(f"n={n} is too small for the Johnson abundance mixture")
    raw = np.concatenate([
        rng.exponential(scale=1.0, size=n_low),
        rng.exponential(scale=10.0, size=n_mid),
        rng.exponential(scale=1000.0, size=n_high),
    ])
    rng.shuffle(raw)
    raw = raw / raw.sum() * total_reads
    counts = np.maximum(1, np.round(raw)).astype(np.int64)
    counts[counts.argmax()] += total_reads - counts.sum()
    assert counts.sum() == total_reads and counts.min() >= 1
    return counts


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

    # At realistic error rates most reads carry no error at all, and every one
    # of them is the same string. Counting those directly leaves the per-base
    # assembly loop to the minority that actually needs it, which is what makes
    # simulations of tens of millions of reads tractable. The draws above are
    # untouched, so a given seed still yields exactly the same reads.
    dirty = subs.any(axis=1) | dels.any(axis=1) | inss.any(axis=1)
    n_clean = int(n_reads - int(dirty.sum()))

    out = Counter()
    if n_clean:
        out[bc] += n_clean

    # Substituted bases depend only on (position, offset), so resolve them once
    # instead of recomputing a modulo and a string index per base per read.
    sub_char = [["ACGT"[(bc_idx[j] + o) % 4] for o in range(4)] for j in range(L)]

    for r in np.flatnonzero(dirty):
        parts = []
        for j in range(L):
            if inss[r, j]:
                parts.append("ACGT"[ins_base[r, j]])
            if dels[r, j]:
                continue
            if subs[r, j]:
                parts.append(sub_char[j][sub_offset[r, j]])
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
    template: str | None = None,
    abundance: str = "lognormal",
    tie_order: str = "generation",
) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bc_rng    = np.random.default_rng(seed)
    abund_rng = np.random.default_rng(seed + 1)
    mut_rng   = np.random.default_rng(seed + 2)

    if template:
        barcodes = structured_barcodes(n_barcodes, template, bc_rng)
        barcode_length = len(template)
    else:
        barcodes = random_barcodes(n_barcodes, barcode_length, bc_rng)
    if abundance == "johnson":
        true_counts = johnson_abundances(n_barcodes, n_reads, abund_rng)
    else:
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

    # Counter preserves insertion order. Because the unmutated true sequence is
    # normally inserted before its error variants, retaining that order leaks
    # ground truth to clustering tools that use input order to break count ties.
    # Keep it as the default for reproduction of historical simulations, while
    # offering sequence order for fair comparisons with deterministic methods.
    if tie_order == "generation":
        items = sorted(all_reads.items(), key=lambda kv: -kv[1])
    elif tie_order == "sequence":
        items = sorted(all_reads.items(), key=lambda kv: (-kv[1], kv[0]))
    else:
        raise ValueError("tie_order must be 'generation' or 'sequence'")
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
        "template":       template,
        "abundance":      abundance,
        "n_unique_reads": len(items),
    }
    if tie_order != "generation":
        manifest["tie_order"] = tie_order
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
    p.add_argument("--template",       type=str,   default=None,
                   help="fixed-anchor design over {A,C,G,T,N}; N = random position")
    p.add_argument("--abundance",      type=str,   default="lognormal",
                   choices=["lognormal", "johnson"],
                   help="abundance model; 'johnson' matches the exponential "
                        "mixture of the Johnson et al. (2023) benchmark library")
    p.add_argument("--tie-order", type=str, default="generation",
                   choices=["generation", "sequence"],
                   help="ordering for equal-count reads; 'generation' reproduces "
                        "historical files, while 'sequence' avoids leaking the "
                        "simulator's truth-first construction order")
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
        template=args.template,
        abundance=args.abundance,
        tie_order=args.tie_order,
    )
