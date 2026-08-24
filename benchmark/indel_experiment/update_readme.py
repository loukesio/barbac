"""Update benchmark tables in README.md from recorded runner summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
METHODS = ("barbac", "shepherd", "starcode", "bartender")
LABELS = {"barbac": "barbac", "shepherd": "Shepherd",
          "starcode": "Starcode", "bartender": "Bartender"}


def selected(frame: pd.DataFrame) -> pd.DataFrame:
    wanted = frame[frame["condition"].isin(("sub_only", "low_indel"))].copy()
    missing = {(condition, method) for condition in ("sub_only", "low_indel")
               for method in METHODS} - set(zip(wanted.condition, wanted.method))
    if missing:
        raise ValueError(f"summary is missing rows: {sorted(missing)}")
    return wanted


def pct(value: float) -> str:
    return f"{100 * value:.2f}"


def unstructured_table(frame: pd.DataFrame) -> str:
    rows = [
        "| Condition | Method | Pearson R | FN% | FP% | WS% | Wall (s) | Algo (s) |",
        "|:--|:--|--:|--:|--:|--:|--:|--:|",
    ]
    frame = selected(frame).set_index(["condition", "method"])
    for condition in ("sub_only", "low_indel"):
        for index, method in enumerate(METHODS):
            row = frame.loc[(condition, method)]
            condition_label = f"**{condition}**" if index == 0 else ""
            rows.append(
                f"| {condition_label} | {LABELS[method]} | {row.pearson_r:.4f} | "
                f"{pct(row.fn_rate)} | {pct(row.fp_rate)} | {pct(row.ws_rate)} | "
                f"{row.runtime_s:.2f} | {row.algo_time_s:.2f} |"
            )
    return "\n".join(rows)


def structured_table(frame: pd.DataFrame) -> str:
    rows: list[str] = []
    frame = selected(frame).set_index(["condition", "method"])
    for condition, title in (("sub_only", "**sub_only (0% indel — representative of Illumina):**"),
                             ("low_indel", "**low_indel (0.5% ins, 0.5% del):**")):
        if rows:
            rows.append("")
        rows.extend([title, "", "| Method | Pearson R | FN% | FP% | WS% |",
                     "|:--|--:|--:|--:|--:|"])
        for method in METHODS:
            row = frame.loc[(condition, method)]
            rows.append(f"| {LABELS[method]} | {row.pearson_r:.4f} | "
                        f"{pct(row.fn_rate)} | {pct(row.fp_rate)} | {pct(row.ws_rate)} |")
    return "\n".join(rows)


def replace_block(text: str, name: str, replacement: str) -> str:
    start = f"<!-- benchmark:{name}:start -->"
    end = f"<!-- benchmark:{name}:end -->"
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError(f"README markers for {name!r} are missing or duplicated")
    before, remainder = text.split(start, 1)
    _, after = remainder.split(end, 1)
    return f"{before}{start}\n{replacement}\n{end}{after}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--structured-summary", type=Path, required=True)
    parser.add_argument("--readme", type=Path, default=REPO / "README.md")
    args = parser.parse_args()

    plain = pd.read_csv(args.summary)
    structured = pd.read_csv(args.structured_summary)
    for path, frame in ((args.summary, plain),
                        (args.structured_summary, structured)):
        if "git_commit" not in frame.columns:
            raise ValueError(
                f"{path} has no git_commit column; regenerate it with the "
                "provenance-aware benchmark runner"
            )
    commits = set(plain.git_commit.dropna()) | set(structured.git_commit.dropna())
    if len(commits) != 1:
        raise ValueError(f"summaries must come from one recorded Git commit: {commits}")

    text = args.readme.read_text()
    text = replace_block(text, "unstructured", unstructured_table(plain))
    text = replace_block(text, "structured", structured_table(structured))
    args.readme.write_text(text)
    print(f"Updated {args.readme} from commit {next(iter(commits))}")


if __name__ == "__main__":
    main()
