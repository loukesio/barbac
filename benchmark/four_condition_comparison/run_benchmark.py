"""Reproduce the four-condition barbac method comparison.

The benchmark covers random and fixed-anchor barcode designs, each with
substitutions only and substitutions plus low-rate insertions/deletions. Pure
substitution conditions run barbac in both Hamming and Levenshtein modes.

Raw/generated files are written below ``generated/`` and are intentionally not
tracked. The compact summary and environment manifest can be compared with the
dated snapshots stored beside this script.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
INDEL_BENCHMARK = REPO_ROOT / "benchmark" / "indel_experiment"
sys.path.insert(0, str(INDEL_BENCHMARK))

from run_experiment import (  # noqa: E402
    BARCODE_LEN,
    BARTENDER_DIR,
    SHEPHERD_DIR,
    STARCODE_BIN,
    evaluate,
    run_barbac,
    run_bartender,
    run_shepherd,
    run_starcode,
)
from simulate import simulate  # noqa: E402


ANCHOR_TEMPLATE = "NNNNNNNNATGCNNNNNNNNATCGTTAA"
CONDITIONS = (
    {
        "slug": "random_substitutions",
        "design": "random",
        "error_model": "substitutions",
        "template": None,
        "sub_rate": 0.005,
        "ins_rate": 0.0,
        "del_rate": 0.0,
    },
    {
        "slug": "random_low_indels",
        "design": "random",
        "error_model": "substitutions_low_indels",
        "template": None,
        "sub_rate": 0.005,
        "ins_rate": 0.005,
        "del_rate": 0.005,
    },
    {
        "slug": "anchored_substitutions",
        "design": "anchored",
        "error_model": "substitutions",
        "template": ANCHOR_TEMPLATE,
        "sub_rate": 0.005,
        "ins_rate": 0.0,
        "del_rate": 0.0,
    },
    {
        "slug": "anchored_low_indels",
        "design": "anchored",
        "error_model": "substitutions_low_indels",
        "template": ANCHOR_TEMPLATE,
        "sub_rate": 0.005,
        "ins_rate": 0.005,
        "del_rate": 0.005,
    },
)


def git_revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def install_current_barbac() -> tuple[tempfile.TemporaryDirectory, Path]:
    temporary = tempfile.TemporaryDirectory(prefix="barbac-four-condition-")
    library = Path(temporary.name)
    subprocess.run(
        [
            "R",
            "CMD",
            "INSTALL",
            "--preclean",
            f"--library={library}",
            str(REPO_ROOT),
        ],
        check=True,
    )
    return temporary, library


def barbac_build_id(library: Path) -> str:
    expression = (
        f'.libPaths(c("{library}", .libPaths())); '
        'library(barbac); cat(barbac:::barbac_build_id())'
    )
    return subprocess.check_output(["Rscript", "-e", expression], text=True).strip()


def result_row(design: str, error_model: str, method: str,
               distance_method: str | None, result) -> dict:
    return {
        "design": design,
        "error_model": error_model,
        "method": method,
        "distance_method": distance_method,
        "n_centroids": result.n_centroids,
        "pearson_r": result.pearson_r,
        "fn": result.fn,
        "fp": result.fp,
        "ws": result.ws,
        "algorithm_time_s": result.algo_time_s,
    }


def run_condition(condition: dict, output_root: Path, library: Path,
                  n_barcodes: int, n_reads: int, seed: int) -> list[dict]:
    condition_dir = output_root / condition["slug"]
    condition_dir.mkdir(parents=True, exist_ok=True)
    barcode_len = (
        len(condition["template"])
        if condition["template"] is not None
        else BARCODE_LEN
    )
    print(f"\n=== {condition['slug']} ===", flush=True)
    simulate(
        condition_dir,
        n_barcodes=n_barcodes,
        barcode_length=barcode_len,
        n_reads=n_reads,
        sub_rate=condition["sub_rate"],
        ins_rate=condition["ins_rate"],
        del_rate=condition["del_rate"],
        sigma=1.5,
        seed=seed,
        template=condition["template"],
        abundance="lognormal",
        tie_order="sequence",
    )

    input_csv = condition_dir / "input.csv"
    truth_csv = condition_dir / "true_counts.csv"
    shepherd_input = condition_dir / "shepherd_input.txt"
    rows: list[dict] = []

    barbac_methods = (
        ("hamming", "lv")
        if condition["ins_rate"] == 0 and condition["del_rate"] == 0
        else ("lv",)
    )
    for distance_method in barbac_methods:
        output = condition_dir / f"barbac_{distance_method}.csv"
        wall, algorithm = run_barbac(
            input_csv,
            output,
            method=distance_method,
            library_path=library,
        )
        scored = evaluate(
            f"barbac_{distance_method}", output, truth_csv, wall,
            algo_time_s=algorithm,
        )
        rows.append(result_row(
            condition["design"], condition["error_model"], "barbac",
            distance_method, scored,
        ))
        print(
            f"  barbac {distance_method}: FN={scored.fn} FP={scored.fp} "
            f"WS={scored.ws} algorithm={algorithm:.2f}s",
            flush=True,
        )

    shepherd_centroids, shepherd_time = run_shepherd(
        shepherd_input, condition_dir, barcode_len
    )
    shepherd_output = condition_dir / "shepherd.csv"
    shepherd_data = pd.read_csv(shepherd_centroids)
    shepherd_data.columns = ["central_barcode", "sum_counts"]
    shepherd_data.to_csv(shepherd_output, index=False)
    shepherd_result = evaluate(
        "shepherd", shepherd_output, truth_csv, shepherd_time
    )
    rows.append(result_row(
        condition["design"], condition["error_model"], "shepherd", None,
        shepherd_result,
    ))

    starcode_output = condition_dir / "starcode.csv"
    starcode_time = run_starcode(shepherd_input, starcode_output)
    starcode_result = evaluate(
        "starcode", starcode_output, truth_csv, starcode_time
    )
    rows.append(result_row(
        condition["design"], condition["error_model"], "starcode", None,
        starcode_result,
    ))

    bartender_output = condition_dir / "bartender.csv"
    bartender_time = run_bartender(input_csv, condition_dir, bartender_output)
    bartender_result = evaluate(
        "bartender", bartender_output, truth_csv, bartender_time
    )
    rows.append(result_row(
        condition["design"], condition["error_model"], "bartender", None,
        bartender_result,
    ))

    for scored in (shepherd_result, starcode_result, bartender_result):
        print(
            f"  {scored.name}: FN={scored.fn} FP={scored.fp} WS={scored.ws} "
            f"runtime={scored.runtime_s:.2f}s",
            flush=True,
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=HERE / "generated")
    parser.add_argument("--n-barcodes", type=int, default=10_000)
    parser.add_argument("--n-reads", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    temporary, library = install_current_barbac()
    try:
        build_id = barbac_build_id(library)
        if "hamming-refinement-v11" not in build_id:
            raise RuntimeError(f"unexpected barbac native build: {build_id}")

        rows: list[dict] = []
        for condition in CONDITIONS:
            rows.extend(run_condition(
                condition,
                args.output_dir,
                library,
                args.n_barcodes,
                args.n_reads,
                args.seed,
            ))

        summary = pd.DataFrame(rows)
        summary.to_csv(args.output_dir / "summary.csv", index=False)
        versions = {
            "barbac_build_id": build_id,
            "barbac_git_revision": git_revision(REPO_ROOT),
            "barbac_worktree_was_clean": not bool(
                subprocess.check_output(
                    ["git", "-C", str(REPO_ROOT), "status", "--porcelain"],
                    text=True,
                ).strip()
            ),
            "shepherd_git_revision": git_revision(SHEPHERD_DIR),
            "starcode_git_revision": git_revision(STARCODE_BIN.parent),
            "bartender_git_revision": git_revision(BARTENDER_DIR),
            "r_version": subprocess.check_output(
                ["Rscript", "-e", "cat(R.version.string)"], text=True
            ).strip(),
            "python_packages": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "rapidfuzz", "scipy")
            },
            "configuration": {
                "n_barcodes": args.n_barcodes,
                "n_reads": args.n_reads,
                "seed": args.seed,
                "tie_order": "sequence",
                "abundance": "lognormal",
                "sub_rate": 0.005,
                "low_ins_rate": 0.005,
                "low_del_rate": 0.005,
                "anchor_template": ANCHOR_TEMPLATE,
            },
        }
        (args.output_dir / "versions.json").write_text(
            json.dumps(versions, indent=2) + "\n"
        )
        print(f"\nWrote {args.output_dir / 'summary.csv'}")
        print(summary.to_string(index=False))
    finally:
        temporary.cleanup()


if __name__ == "__main__":
    main()
