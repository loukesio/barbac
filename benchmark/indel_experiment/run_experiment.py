"""Run super_cluster2 vs Shepherd vs Starcode vs Bartender across a sweep
of noise conditions.

For each condition: simulate -> run all methods -> evaluate -> log timings.
Writes results/<condition>/{barbac_out.csv, shepherd_out.csv, starcode_out.csv,
bartender_out.csv, eval.json} and a summary.csv with all conditions
side-by-side.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
TOOLS_DIR = Path("~/Documents/Projects/Barcodes/barbac-benchmark/tools").expanduser()
SHEPHERD_DIR   = TOOLS_DIR / "Shepherd"
STARCODE_BIN   = TOOLS_DIR / "starcode" / "starcode"
BARTENDER_DIR  = TOOLS_DIR / "bartender-1.1"
BARTENDER_BIN  = BARTENDER_DIR / "bartender_single_com"
MAX_DIST = 3
BARCODE_LEN = 20

# Conditions: sub_rate held at 0.5%, indel rate swept.
CONDITIONS = [
    {"name": "sub_only",   "sub_rate": 0.005, "ins_rate": 0.000, "del_rate": 0.000},
    {"name": "low_indel",  "sub_rate": 0.005, "ins_rate": 0.005, "del_rate": 0.005},
    {"name": "mid_indel",  "sub_rate": 0.005, "ins_rate": 0.020, "del_rate": 0.020},
    {"name": "high_indel", "sub_rate": 0.005, "ins_rate": 0.050, "del_rate": 0.050},
    {"name": "nanopore",   "sub_rate": 0.010, "ins_rate": 0.080, "del_rate": 0.080},
]


@dataclass
class MethodResult:
    name:        str
    n_centroids: int
    pearson_r:   float
    fn:          int
    fn_rate:     float
    fp:          int
    fp_rate:     float
    ws:          int
    ws_rate:     float
    runtime_s:   float
    # Pure algorithm time, excluding process/interpreter boot. For
    # methods that don't distinguish (Shepherd, Starcode, Bartender) this
    # is set equal to runtime_s. For barbac we parse it out of the R
    # script's stdout so the R + package load overhead is separated.
    algo_time_s: float = 0.0


def evaluate(name: str, centroids_path: Path, true_counts_path: Path,
             runtime_s: float, max_dist: int = MAX_DIST,
             algo_time_s: float | None = None) -> MethodResult:
    result = pd.read_csv(centroids_path)
    result.columns = ["central_barcode", "sum_counts"]

    true_counts = pd.read_csv(true_counts_path)
    true_counts.columns = ["barcode", "true_count"]

    true_set  = set(true_counts["barcode"])
    true_list = true_counts["barcode"].tolist()
    n_true    = len(true_set)

    centroid_set = set(result["central_barcode"])
    fn = list(true_set - centroid_set)
    fp = list(centroid_set - true_set)

    matched = result.merge(
        true_counts, left_on="central_barcode", right_on="barcode", how="inner",
    )
    if len(matched) >= 2:
        r = float(np.corrcoef(
            np.log10(matched["sum_counts"]),
            np.log10(matched["true_count"]),
        )[0, 1])
    else:
        r = float("nan")

    ws = 0
    if fp:
        batch = 200
        for i in range(0, len(fp), batch):
            chunk = fp[i:i + batch]
            dm = process.cdist(
                chunk, true_list,
                scorer=Levenshtein.distance,
                score_cutoff=max_dist,
                dtype=np.uint8,
            )
            ws += int((dm.min(axis=1) <= max_dist).sum())

    return MethodResult(
        name=name,
        n_centroids=len(result),
        pearson_r=r,
        fn=len(fn), fn_rate=len(fn) / n_true,
        fp=len(fp), fp_rate=len(fp) / n_true,
        ws=ws,      ws_rate=ws / n_true,
        runtime_s=runtime_s,
        algo_time_s=runtime_s if algo_time_s is None else algo_time_s,
    )


def run_barbac(input_csv: Path, out_csv: Path) -> tuple[float, float]:
    """Returns (wall_seconds, algo_seconds).

    algo_seconds is measured inside the R script and excludes R boot +
    package loading. It is parsed out of the child's stdout line
    ``BARBAC_ALGO_SECONDS=<x>``.
    """
    script = HERE / "run_barbac.R"
    t0 = time.time()
    proc = subprocess.run(
        ["Rscript", str(script), str(input_csv), str(out_csv),
         str(MAX_DIST), "20"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    wall = time.time() - t0
    algo = wall
    for line in proc.stdout.decode(errors="replace").splitlines():
        if line.startswith("BARBAC_ALGO_SECONDS="):
            try:
                algo = float(line.split("=", 1)[1])
            except ValueError:
                pass
    return wall, algo


def run_shepherd(shepherd_input: Path, work_dir: Path) -> tuple[Path, float]:
    """Returns (centroids_csv, wall_seconds). Shepherd writes outputs alongside
    its input file with a fixed suffix, so we stage by copying / linking and
    then move outputs back."""
    # Shepherd writes outputs in the working directory based on the input
    # filename stem: <stem>_pb_freq.csv, <stem>_seq_clust.csv, etc.
    stem = shepherd_input.stem  # e.g. "shepherd_input"
    t0 = time.time()
    subprocess.run(
        [sys.executable, str(SHEPHERD_DIR / "shepherd_t0.py"),
         "-f", str(shepherd_input),
         "-l", str(BARCODE_LEN),
         "-eps", str(MAX_DIST)],
        check=True,
        cwd=work_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    dt = time.time() - t0
    centroids = work_dir / f"{stem}_pb_freq.csv"
    return centroids, dt


def run_starcode(shepherd_input: Path, out_csv: Path) -> float:
    """Runs starcode -d 3 -s and writes a (central_barcode, sum_counts) CSV.

    Reuses shepherd_input.txt as-is because starcode accepts the same
    TAB-separated `sequence<TAB>count` layout.
    """
    tmp_out = out_csv.with_suffix(".starcode.tsv")
    t0 = time.time()
    subprocess.run(
        [str(STARCODE_BIN),
         "-d", str(MAX_DIST),
         "-s",
         "-i", str(shepherd_input),
         "-o", str(tmp_out)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    dt = time.time() - t0
    df = pd.read_csv(tmp_out, sep="\t", header=None,
                     names=["central_barcode", "sum_counts"])
    df.to_csv(out_csv, index=False)
    return dt


def run_bartender(input_csv: Path, work_dir: Path, out_csv: Path) -> float:
    """Runs bartender_single_com -d 3.

    Bartender treats each input row as one read (with the second column
    as a UMI/line-number), so we expand `(barcode, count)` -> one row per
    read first, then invoke bartender. Output <prefix>_cluster.csv is
    reshaped into the standard (central_barcode, sum_counts) layout using
    the `time_point_1` column as the summed count.
    """
    expanded = work_dir / "bartender_expanded.csv"
    with input_csv.open() as fh, expanded.open("w") as out:
        next(fh)  # skip header
        line_id = 0
        for line in fh:
            seq, cnt = line.rstrip("\n").split(",")
            for _ in range(int(cnt)):
                line_id += 1
                out.write(f"{seq},{line_id}\n")

    prefix = work_dir / "bartender_out"
    env = os.environ.copy()
    env["PATH"] = f"{BARTENDER_DIR}:{env.get('PATH', '')}"

    t0 = time.time()
    subprocess.run(
        [str(BARTENDER_BIN),
         "-f", str(expanded),
         "-o", str(prefix),
         "-d", str(MAX_DIST)],
        check=True,
        cwd=work_dir,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    dt = time.time() - t0
    df = pd.read_csv(work_dir / "bartender_out_cluster.csv")
    df = df.rename(columns={"Center": "central_barcode",
                            "time_point_1": "sum_counts"})[
        ["central_barcode", "sum_counts"]
    ]
    df.to_csv(out_csv, index=False)
    # Clean up the expanded file - it can be ~30MB for 1M reads
    try:
        expanded.unlink()
    except FileNotFoundError:
        pass
    return dt


def run_condition(cond: dict, *, n_barcodes: int, n_reads: int, seed: int,
                  reuse_sim: bool = False) -> dict:
    cond_dir = RESULTS_DIR / cond["name"]
    cond_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {cond['name']}  sub={cond['sub_rate']}  ins={cond['ins_rate']}  del={cond['del_rate']} ===")

    # 1. Simulate (or reuse).
    from simulate import simulate
    if not reuse_sim or not (cond_dir / "input.csv").exists():
        t0 = time.time()
        simulate(
            cond_dir,
            n_barcodes=n_barcodes,
            barcode_length=BARCODE_LEN,
            n_reads=n_reads,
            sub_rate=cond["sub_rate"],
            ins_rate=cond["ins_rate"],
            del_rate=cond["del_rate"],
            seed=seed,
        )
        print(f"  sim:      {time.time()-t0:.1f}s")
    else:
        print(f"  sim:      reused")

    true_csv  = cond_dir / "true_counts.csv"
    input_csv = cond_dir / "input.csv"
    shep_input = cond_dir / "shepherd_input.txt"

    # 2. barbac.
    barbac_out = cond_dir / "barbac_out.csv"
    barbac_wall, barbac_algo = run_barbac(input_csv, barbac_out)
    barbac_eval = evaluate("barbac", barbac_out, true_csv,
                           barbac_wall, algo_time_s=barbac_algo)
    print(f"  barbac:   wall={barbac_wall:.1f}s algo={barbac_algo:.1f}s  "
          f"R={barbac_eval.pearson_r:.4f}  "
          f"FN={barbac_eval.fn} ({100*barbac_eval.fn_rate:.3f}%)  "
          f"FP={barbac_eval.fp} ({100*barbac_eval.fp_rate:.3f}%)  "
          f"WS={barbac_eval.ws} ({100*barbac_eval.ws_rate:.3f}%)")

    # 3. Shepherd.
    shep_centroids, shep_t = run_shepherd(shep_input, cond_dir)
    shep_out = cond_dir / "shepherd_out.csv"
    # normalize column names for downstream eval
    df = pd.read_csv(shep_centroids)
    df.columns = ["central_barcode", "sum_counts"]
    df.to_csv(shep_out, index=False)
    shep_eval = evaluate("shepherd", shep_out, true_csv, shep_t)
    print(f"  shepherd: {shep_t:.1f}s  R={shep_eval.pearson_r:.4f}  "
          f"FN={shep_eval.fn} ({100*shep_eval.fn_rate:.3f}%)  "
          f"FP={shep_eval.fp} ({100*shep_eval.fp_rate:.3f}%)  "
          f"WS={shep_eval.ws} ({100*shep_eval.ws_rate:.3f}%)")

    # 4. Starcode.
    starcode_out = cond_dir / "starcode_out.csv"
    starcode_t = run_starcode(shep_input, starcode_out)
    starcode_eval = evaluate("starcode", starcode_out, true_csv, starcode_t)
    print(f"  starcode: {starcode_t:.1f}s  R={starcode_eval.pearson_r:.4f}  "
          f"FN={starcode_eval.fn} ({100*starcode_eval.fn_rate:.3f}%)  "
          f"FP={starcode_eval.fp} ({100*starcode_eval.fp_rate:.3f}%)  "
          f"WS={starcode_eval.ws} ({100*starcode_eval.ws_rate:.3f}%)")

    # 5. Bartender.
    bart_out = cond_dir / "bartender_out.csv"
    bart_t = run_bartender(input_csv, cond_dir, bart_out)
    bart_eval = evaluate("bartender", bart_out, true_csv, bart_t)
    print(f"  bartender:{bart_t:.1f}s  R={bart_eval.pearson_r:.4f}  "
          f"FN={bart_eval.fn} ({100*bart_eval.fn_rate:.3f}%)  "
          f"FP={bart_eval.fp} ({100*bart_eval.fp_rate:.3f}%)  "
          f"WS={bart_eval.ws} ({100*bart_eval.ws_rate:.3f}%)")

    eval_dict = {
        "condition": cond,
        "barbac":    barbac_eval.__dict__,
        "shepherd":  shep_eval.__dict__,
        "starcode":  starcode_eval.__dict__,
        "bartender": bart_eval.__dict__,
    }
    (cond_dir / "eval.json").write_text(json.dumps(eval_dict, indent=2))
    return eval_dict


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--n-barcodes", type=int, default=10_000)
    p.add_argument("--n-reads",    type=int, default=1_000_000)
    p.add_argument("--seed",       type=int, default=42)
    p.add_argument("--only",       type=str, help="only run this condition name")
    p.add_argument("--reuse-sim",  action="store_true",
                   help="reuse existing input.csv if present")
    args = p.parse_args()

    conds = CONDITIONS
    if args.only:
        conds = [c for c in conds if c["name"] == args.only]
        if not conds:
            sys.exit(f"unknown condition {args.only!r}")

    all_results = []
    for cond in conds:
        all_results.append(run_condition(
            cond,
            n_barcodes=args.n_barcodes,
            n_reads=args.n_reads,
            seed=args.seed,
            reuse_sim=args.reuse_sim,
        ))

    # Compile summary
    rows = []
    for r in all_results:
        c = r["condition"]
        for m in ("barbac", "shepherd", "starcode", "bartender"):
            if m not in r:
                continue
            row = {"condition": c["name"], "method": m,
                   "sub_rate": c["sub_rate"], "ins_rate": c["ins_rate"], "del_rate": c["del_rate"]}
            row.update({k: r[m].get(k, r[m]["runtime_s"]) if k == "algo_time_s"
                        else r[m][k]
                        for k in ("n_centroids", "pearson_r", "fn", "fn_rate",
                                  "fp", "fp_rate", "ws", "ws_rate",
                                  "runtime_s", "algo_time_s")})
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary_path = RESULTS_DIR / "summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\nWrote {summary_path}")
    print("\n=== SUMMARY ===")
    cols = ["condition", "method", "pearson_r", "fn_rate", "fp_rate",
            "ws_rate", "runtime_s", "algo_time_s"]
    with pd.option_context("display.float_format", "{:.4f}".format):
        print(summary[cols].to_string(index=False))


if __name__ == "__main__":
    main()
