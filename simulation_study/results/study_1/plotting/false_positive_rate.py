from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))
from results.utils import (
    aggregate_significance,
    apply_minimal_style,
    apply_yaxes_significance_style,
    COLOR_MAP,
    LINETYPE_MAP,
    LABELSIZE,
    PANEL_SIZE, 
    MARKERS
)

CSV_PATH = Path("study_1_results_new.csv")
OUT_DIR = Path("results/study_1/plotting/false_positive_rate")
OUT_DIR.mkdir(parents=True, exist_ok=True)
CUSTOM_ORDER: List[str] = ["solid", "longdash", "dotted", "dotdash", "twodash", "dashed"]
LEGEND_SIZE = (14, 1)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    max_bmxr = df["beta_mxr"].max()
    df = df[(df["beta_ym"] == 0.0) & (df["beta_mxr"] == max_bmxr)].copy()
    allowed_em = [0.0, 0.3, 0.6, 0.9]
    allowed_conf = [0.0, 0.2, 0.4, 0.6]
    df = df[df["effect_modification"].isin(allowed_em) &
            df["confounding"].isin(allowed_conf)].copy()
    df["method"] = df["method"].apply(
        lambda m: "ML Estimator" if str(m) == "standard" else "G-Estimator"
    )
    return df


def plot_per_confounding(df: pd.DataFrame, conf: float) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    apply_minimal_style(ax)
    apply_yaxes_significance_style(ax)
    ax.axhline(0.05, color=COLOR_MAP.get("Threshold", "black"), linewidth=0.8)
    ax.set_yticks([0.05, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Sample Size (N)", fontsize=LABELSIZE)
    if conf== 0:
        ax.set_ylabel("False Positive Rate", fontsize=LABELSIZE)
    else:
        ax.set_ylabel(" ", fontsize=LABELSIZE)

    df_subset = df[df["confounding"] == conf]
    significance_df = aggregate_significance(
        df_subset, group_vars=["num_obs", "method", "effect_modification"]
    )

    effect_levels = sorted(df_subset["effect_modification"].unique())

    for (method, eff), g in significance_df.groupby(["method", "effect_modification"], sort=False):
        x = g["num_obs"].to_numpy()
        y = g["significance"].to_numpy()
        color = COLOR_MAP.get(method, "black")

        idx = effect_levels.index(eff)
        style_name = CUSTOM_ORDER[idx % len(CUSTOM_ORDER)]
        dash = LINETYPE_MAP.get(style_name, (0, (None, None)))
        marker = MARKERS[idx % len(MARKERS)]

        line, = ax.plot(
            x, y,
            label=f"{method}, {eff}",
            linewidth=1.2,
            marker=marker,
            markersize=5,
            color=color,
        )
        if dash[1] is not None:
            line.set_dashes(dash[1])

    return fig


def save_all_confounding_plots(df: pd.DataFrame) -> List[str]:
    saved: List[str] = []
    conf_levels = sorted(df["confounding"].unique())

    for conf in conf_levels:
        fig = plot_per_confounding(df, conf)
        out_path = OUT_DIR / f"false_positive_rate_conf_{conf}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        saved.append(str(out_path))

    legend = save_legend(df)
    saved.append(legend)
    return saved


def save_legend(df: pd.DataFrame) -> str:
    methods = sorted(df["method"].astype(str).unique())
    effect_levels = sorted(df["effect_modification"].unique())

    # Method handles (colors)
    method_handles = [
        Line2D([0], [0], linestyle="-", linewidth=1.2, marker="o", color=COLOR_MAP.get(m, "black"), label=m)
        for m in methods
    ]

    # Effect handles (linestyles + markers)
    effect_handles = []
    for i, eff in enumerate(effect_levels):
        style_name = CUSTOM_ORDER[i % len(CUSTOM_ORDER)]
        dash = LINETYPE_MAP.get(style_name, (0, (None, None)))
        marker = MARKERS[i % len(MARKERS)]
        h = Line2D(
            [0], [0], linestyle="-", linewidth=1.2, color="black", marker=marker, markersize=5, label=f"{eff}"
        )
        if dash[1] is not None:
            h.set_dashes(dash[1])
        effect_handles.append(h)

    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    ax.axis("off")
    leg1 = ax.legend(
        handles=method_handles,
        loc="center",
        title="Method",
        fontsize=LABELSIZE - 4,
        title_fontsize=LABELSIZE,
        bbox_to_anchor=(0.2, 0.5),
        ncol=len(method_handles),
        frameon=False,
    )
    ax.add_artist(leg1)
    leg2 = ax.legend(
        handles=effect_handles,
        title=r"Effect Modification $\delta_{ur}$",
        loc="center",
        fontsize=LABELSIZE - 4,
        title_fontsize=LABELSIZE,
        bbox_to_anchor=(0.7, 0.5),
        ncol=len(effect_handles),
        frameon=False,
    )
    ax.add_artist(leg2)

    out_path = OUT_DIR / "false_positive_rate_legend.pdf"
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(out_path)


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV not found at {CSV_PATH.resolve()}")
    df = pd.read_csv(CSV_PATH)
    df = prepare_data(df)
    saved = save_all_confounding_plots(df)
    print("Saved:")
    for p in saved:
        print(" -", p)


if __name__ == "__main__":
    main()
