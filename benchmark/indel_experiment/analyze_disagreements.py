"""Explain centroid disagreements between two barcode clustering results.

The regular benchmark reports aggregate FN/FP/WS counts.  This companion
script reports the individual true barcodes recovered by one method but not
the other, together with their observed counts and closest centroids.  It is
intended for diagnosing merge rules, not for scoring a method.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Hamming


def _read_two_columns(path: Path, names: tuple[str, str]) -> pd.DataFrame:
    data = pd.read_csv(path)
    if data.shape[1] != 2:
        raise ValueError(f"{path} must contain exactly two columns")
    data.columns = list(names)
    return data


def _closest(sequence: str, centroids: list[str], max_distance: int) -> tuple[str | None, int | None]:
    matches = process.extract(
        sequence,
        centroids,
        scorer=Hamming.distance,
        score_cutoff=max_distance,
        limit=1,
    )
    if not matches:
        return None, None
    centroid, distance, _ = matches[0]
    return centroid, int(distance)


def explain(
    truth_path: Path,
    input_path: Path,
    left_path: Path,
    right_path: Path,
    max_distance: int,
) -> pd.DataFrame:
    truth = _read_two_columns(truth_path, ("barcode", "true_count"))
    observed = _read_two_columns(input_path, ("barcode", "observed_count"))
    left = _read_two_columns(left_path, ("centroid", "cluster_count"))
    right = _read_two_columns(right_path, ("centroid", "cluster_count"))

    truth_set = set(truth["barcode"])
    left_set = set(left["centroid"])
    right_set = set(right["centroid"])
    observed_count = observed.set_index("barcode")["observed_count"].to_dict()
    observed_position = {
        barcode: position for position, barcode in enumerate(observed["barcode"])
    }
    true_count = truth.set_index("barcode")["true_count"].to_dict()
    left_cluster_count = left.set_index("centroid")["cluster_count"].to_dict()
    right_cluster_count = right.set_index("centroid")["cluster_count"].to_dict()
    left_centroids = left["centroid"].tolist()
    right_centroids = right["centroid"].tolist()

    rows: list[dict[str, object]] = []
    categories = (
        ("right_only", (truth_set & right_set) - left_set),
        ("left_only", (truth_set & left_set) - right_set),
        ("missed_by_both", truth_set - (left_set | right_set)),
    )
    for category, barcodes in categories:
        for barcode in sorted(barcodes):
            left_parent, left_distance = _closest(barcode, left_centroids, max_distance)
            right_parent, right_distance = _closest(barcode, right_centroids, max_distance)
            barcode_position = observed_position.get(barcode)
            left_position = observed_position.get(left_parent)
            right_position = observed_position.get(right_parent)
            rows.append(
                {
                    "category": category,
                    "barcode": barcode,
                    "true_count": int(true_count[barcode]),
                    "observed_count": int(observed_count.get(barcode, 0)),
                    "input_position": barcode_position,
                    "left_nearest": left_parent,
                    "left_distance": left_distance,
                    "left_cluster_count": left_cluster_count.get(left_parent),
                    "left_observed_count": observed_count.get(left_parent),
                    "left_input_position": left_position,
                    "input_before_left": (
                        barcode_position < left_position
                        if barcode_position is not None and left_position is not None
                        else None
                    ),
                    "lex_before_left": barcode < left_parent if left_parent else None,
                    "right_nearest": right_parent,
                    "right_distance": right_distance,
                    "right_cluster_count": right_cluster_count.get(right_parent),
                    "right_observed_count": observed_count.get(right_parent),
                    "right_input_position": right_position,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--max-distance", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = explain(
        args.truth,
        args.input,
        args.left,
        args.right,
        args.max_distance,
    )
    if args.output:
        result.to_csv(args.output, index=False)
    counts = result.groupby("category", sort=False).size()
    print(counts.to_string())
    print()
    print(
        result.groupby("category", sort=False)[["true_count", "observed_count"]]
        .describe(percentiles=[0.25, 0.5, 0.75])
        .to_string()
    )
    right_only = result[result["category"] == "right_only"]
    with_left = right_only[right_only["left_nearest"].notna()]
    if not with_left.empty:
        tied = with_left["observed_count"] == with_left["left_observed_count"]
        print()
        print("right_only tie diagnostics")
        print(f"  same observed count as left neighbor: {int(tied.sum())}/{len(with_left)}")
        print(f"  input before left neighbor: {int(with_left['input_before_left'].sum())}/{len(with_left)}")
        print(f"  lexicographically before left neighbor: {int(with_left['lex_before_left'].sum())}/{len(with_left)}")


if __name__ == "__main__":
    main()
