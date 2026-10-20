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
from simulation_study.results.utils import (
    apply_minimal_style,
    compute_offsets,
    COLOR_MAP,
    LABELSIZE,
    PANEL_SIZE,
    LEGEND_SIZE,
    METHOD_MAP,
    MARKER_MAP
)

CSV_PATHS = [
    Path("simulation_study/results/study_1/study_1_results.csv"),
]
OUT_DIR = Path("simulation_study/results/study_1/plots/mse")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CUSTOM_ORDER: List[str] = ["solid", "longdash", "dotted", "dotdash", "twodash", "dashed"]


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["variable"] == "med"].copy()
    max_bmxr = df["interaction"].max()
    # confounding scenarios of study 1 (no mediator-outcome confounding)
    df = df[(df["true_theta_m"] == 0.0) & (df["interaction"] == max_bmxr)
            & (df["conf_med"] == 0.0)].copy()
    df["confounding"] = df["confounding"].round(3)
    df["method"] = df["method"].astype(str).map(METHOD_MAP).fillna(df["method"])
    return df


def compute_mse(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["squared_error"] = (df["true_theta_m"] - df["mean"]) ** 2
    mse_df = (
        df.groupby(["num_obs", "method"], as_index=False)
        .agg(
            mse=("squared_error", "mean"),
            sd=("squared_error", "std"),
            n=("squared_error", "count")
        )
    )
    mse_df["se"] = mse_df["sd"] / (mse_df["n"] ** 0.5)
    mse_df["ci_low"] = 1.96 * mse_df["se"]
    mse_df["ci_high"] = 1.96 * mse_df["se"]
    return mse_df


def plot_per_confounding(df: pd.DataFrame, conf: float) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    apply_minimal_style(ax)
    ax.tick_params(
        axis='both', 
        which='major', 
        length=6, 
        width=1.5, 
        labelsize=LABELSIZE - 4
    )
    ax.set_xlabel("Sample Size", fontsize=LABELSIZE)
    ax.set_xticks([250, 500, 750, 1000])
    ytick_labels = [0.00, 0.02, 0.04, 0.06, 0.08, 0.1]
    for y in ytick_labels:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=1)
    ax.set_ylim(0, 0.11)
    if conf == 0.1:
         ax.set_yticks(ytick_labels)
         ax.set_ylabel("Mean Squared Error", fontsize=LABELSIZE)
    else:
        ax.set_yticks([])
        ax.set_ylabel("")
    df_subset = df[df["confounding"] == conf]
    mse_df = compute_mse(df_subset)
    n_groups = len(METHOD_MAP)
    offsets = compute_offsets(n_groups)
    method_order = list(METHOD_MAP.values())
    methods = [m for m in method_order if m in mse_df["method"].unique()]
    mse_df["method"] = pd.Categorical(mse_df["method"], categories=methods, ordered=True)
    for idx, (method, g) in enumerate(mse_df.groupby("method", sort=True, observed=True)):
        x = g["num_obs"].to_numpy()
        y = g["mse"].to_numpy()
        yerr = [g["ci_low"].to_numpy(), g["ci_high"].to_numpy()]
        color = COLOR_MAP.get(method, "black")
        x_jittered = x + offsets[idx]
        ax.errorbar(
            x_jittered, y, yerr=yerr,
            label=method,
            color=color,
            marker=MARKER_MAP.get(method, "o"),
            markersize=8,
            linestyle="None",
            capsize=5,
            elinewidth=2.5
        )
    return fig


def save_all_confounding_plots(df: pd.DataFrame) -> List[str]:
    saved: List[str] = []
    conf_levels = sorted(df["confounding"].unique())
    for conf in conf_levels:
        fig = plot_per_confounding(df, conf)
        out_path = OUT_DIR / f"mse_conf_{conf}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        saved.append(str(out_path))
    legend = save_legend(df)
    saved.append(legend)
    return saved


def save_legend(df: pd.DataFrame) -> str:
    method_order = list(METHOD_MAP.values())
    methods = [m for m in method_order if m in df["method"].unique()]
    method_handles = [
        Line2D([0], [0], linestyle="-", linewidth=2.5, marker=MARKER_MAP.get(m, "o"), markersize=8,
               color=COLOR_MAP.get(m, "black"), label=m)
        for m in methods
    ]
    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    ax.axis("off")
    ax.legend(
        handles=method_handles,
        loc="center",
        title="Method",
        fontsize=LABELSIZE - 2,
        title_fontsize=LABELSIZE,
        ncol=len(method_handles),
        frameon=False,
    )
    out_path = OUT_DIR / "mse_legend.pdf"
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(out_path)


def main():
    missing = [p for p in CSV_PATHS if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing CSV files: {missing}")
    dfs = [pd.read_csv(p) for p in CSV_PATHS]
    df = pd.concat(dfs, ignore_index=True)
    df = prepare_data(df)
    saved = save_all_confounding_plots(df)
    print("Saved:")
    for p in saved:
        print(" -", p)


if __name__ == "__main__":
    main()