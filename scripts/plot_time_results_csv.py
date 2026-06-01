#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


NAME_RE = re.compile(r"^(?P<family>.+?)_(?P<edges>\d+e\d+)(?:_(?P<tag>-1|0|t))?$")


def candidate_paths(value: str) -> list[Path]:
    path = Path(value)
    if path.is_absolute():
        return [path.resolve()]

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    return [
        (Path.cwd() / path).resolve(),
        (script_dir / path).resolve(),
        (project_root / path).resolve(),
    ]


def resolve_existing_path(value: str) -> Path:
    for candidate in candidate_paths(value):
        if candidate.exists():
            return candidate
    return candidate_paths(value)[0]


def parse_graph_name(graph_name: str) -> tuple[str, str, float]:
    match = NAME_RE.match(graph_name)
    if not match:
        raise ValueError(f"Unrecognized graph_name format: {graph_name}")
    family = match.group("family")
    tag = match.group("tag") or "nothing"
    edges = float(match.group("edges"))
    return family, tag, edges


def fit_power_law(x: np.ndarray, y: np.ndarray) -> dict[str, float] | None:
    """Fit a power law y = a * x^b using OLS on log10-transformed data.

    Returns a dict with keys `a` and `b`, or `None` if insufficient data.
    """
    mask = (x > 0) & (y > 0)
    if np.count_nonzero(mask) < 2:
        return None

    x_log = np.log10(x[mask])
    y_log = np.log10(y[mask])
    slope, intercept, _, _, _ = stats.linregress(x_log, y_log)

    return {
        "a": 10 ** intercept,
        "b": slope,
    }




def plot_group_overlay(
    grouped_series: dict[str, pd.DataFrame], family: str, tag: str, output_dir: Path
) -> Path:
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = ["tab:orange", "tab:green", "tab:blue", "tab:red", "tab:purple"]

    title_tag = "original" if tag == "nothing" else tag
    plotted_any = False

    for (idx, (label, group)) in enumerate(grouped_series.items()):
        merged_rows = []
        for version in ["v1", "v2", "v3"]:
            version_df = group[["edges", version]].dropna()
            if not version_df.empty:
                merged_rows.append(version_df.rename(columns={version: "time"}))

        if not merged_rows:
            continue

        merged = pd.concat(merged_rows, ignore_index=True)
        merged = merged.groupby("edges", as_index=False)["time"].mean().sort_values("edges")
        x = merged["edges"].to_numpy(dtype=float)
        # convert milliseconds to seconds for plotting
        y = merged["time"].to_numpy(dtype=float) / 1000.0

        color = colors[idx % len(colors)]
        ax.scatter(x, y, s=5, color=color)

        fit = fit_power_law(x, y)
        if fit is not None:
            a = fit["a"]
            b = fit["b"]
            x_fit = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
            y_fit = a * np.power(x_fit, b)
            ax.plot(x_fit, y_fit, color=color, label=f"{label} ({a:.2e} $m^{{{b:.2f}}}$)")
            print(f"{family} ({title_tag}) [{label}]: a={a:.3e}, b={b:.3f}")
        else:
            ax.plot(x, y, color=color, linestyle="--", label=f"{label} trend")
            print(f"{family} ({title_tag}) [{label}]: not enough data for a fit")

        plotted_any = True

    if not plotted_any:
        plt.close(fig)
        raise ValueError(f"No plottable data for {family} / {tag}")

    ax.set_title(f"{family} ({title_tag})", fontsize=10)
    ax.set_xlabel("#edges")
    ax.set_ylabel("time (s)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.grid(True, which="both", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.legend(fontsize=10)

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = "" if tag == "nothing" else f"_{tag}"
    output_path = output_dir / f"{family}{suffix}.pdf"
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def load_groups(csv_path: Path) -> dict[tuple[str, str], pd.DataFrame]:
    df = pd.read_csv(csv_path)
    required = {"graph_name", "v1", "v2", "v3"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    rows = []
    for _, row in df.iterrows():
        family, tag, edges = parse_graph_name(str(row["graph_name"]))
        rows.append(
            {
                "family": family,
                "tag": tag,
                "edges": edges,
                "v1": pd.to_numeric(row["v1"], errors="coerce"),
                "v2": pd.to_numeric(row["v2"], errors="coerce"),
                "v3": pd.to_numeric(row["v3"], errors="coerce"),
            }
        )

    parsed = pd.DataFrame(rows)
    groups: dict[tuple[str, str], pd.DataFrame] = {}
    for (family, tag), group in parsed.groupby(["family", "tag"], sort=True):
        groups[(family, tag)] = group.reset_index(drop=True)
    return groups


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot v1/v2/v3 results from a CSV")
    parser.add_argument(
        "--input",
        default="../data/plots/Raport/GORC",
        help="Input CSV with graph_name,v1,v2,v3",
    )
    parser.add_argument(
        "--output-dir",
        default="../data/plots/Raport",
        help="Directory for generated PDFs",
    )
    parser.add_argument(
        "--family",
        default=None,
        help="Optional family filter, e.g. big_aug_gor or a comma-separated list",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Optional tag filter, e.g. original, t, 0, -1",
    )
    parser.add_argument(
        "--series",
        action="append",
        default=None,
        help="Optional labeled CSV series in the form LABEL=PATH; repeat to overlay multiple datasets",
    )
    args = parser.parse_args()

    requested_families = None
    if args.family:
        requested_families = {part.strip() for part in args.family.split(",") if part.strip()}

    output_dir = resolve_existing_path(args.output_dir)

    if args.series:
        series_specs = []
        for item in args.series:
            if "=" in item:
                label, path_value = item.split("=", 1)
                label = label.strip()
                path_value = path_value.strip()
            else:
                path_value = item.strip()
                label = Path(path_value).stem
            csv_path = resolve_existing_path(path_value)
            if not csv_path.exists():
                raise FileNotFoundError(f"Input CSV not found: {csv_path}")
            series_specs.append((label, load_groups(csv_path)))
    else:
        csv_path = resolve_existing_path(args.input)
        if not csv_path.exists():
            raise FileNotFoundError(f"Input CSV not found: {csv_path}")
        series_specs = [(Path(csv_path).stem, load_groups(csv_path))]

    outputs = []
    all_keys = sorted({key for _, groups in series_specs for key in groups.keys()})
    for family, tag in all_keys:
        if requested_families and family not in requested_families:
            continue
        if args.tag:
            requested_tag = "nothing" if args.tag == "original" else args.tag
            if tag != requested_tag:
                continue

        grouped_series = {
            label: groups[(family, tag)]
            for label, groups in series_specs
            if (family, tag) in groups
        }
        # Use the consolidated overlay plotting function for both single
        # and multiple series to avoid duplicated code paths.
        outputs.append(plot_group_overlay(grouped_series, family, tag, output_dir))

    print(f"Wrote {len(outputs)} plot(s) to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())