from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))
from results.utils import apply_minimal_style, COLOR_MAP
LABELSIZE = 26

CSV_PATH = Path("study_1_results.csv")
OUT_DIR = Path("results/study_1/plotting/bias_histograms")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df[(df["beta_ym"] == 0.0) & (df["effect_modification"] == 0)].copy()
    df["method"] = df["method"].replace({"standard": "ML Estimator"})
    df.loc[df["method"] != "ML Estimator", "method"] = "G-Estimator"
    return df


def get_subplot_grid(n_items: int, nrows: int = 2) -> tuple[int, int]:
    ncols = math.ceil(n_items / nrows)
    return nrows, ncols


def plot_panel(ax, sub_df: pd.DataFrame, n_obs: int, c: int):
    apply_minimal_style(ax)

    methods = sub_df["method"].unique()[::-1]
    bin_width = 0.06
    min_val, max_val = -1.2, 1.2
    bins = np.arange(min_val, max_val + bin_width, bin_width)

    for method in methods:
        msub = sub_df[sub_df["method"] == method]
        ax.hist(
            msub["bias"],
            bins=bins,
            alpha=0.4,
            color=COLOR_MAP.get(method, "grey"),
            label=method,
            weights=np.ones(len(msub)) / len(msub)
        )

    ax.set_xlim(min_val, max_val)
    ax.set_ylim(0, 0.6)

    ax.set_xlabel("CME Bias", fontsize=LABELSIZE)
    if c == 0:
        ax.set_ylabel("Fraction", fontsize=LABELSIZE)
    else:
        ax.set_ylabel(" ")
        ax.set_yticklabels([])
    ax.tick_params(axis="both", which="major", labelsize=LABELSIZE - 4)

    return ax.get_legend_handles_labels()


def save_legend_plot(handles, labels, out_path: Path):
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.axis("off")
    ax.legend(
        handles, labels,
        title="Estimator",
        loc="center",
        fontsize=LABELSIZE - 4,
        title_fontsize=LABELSIZE,
        bbox_to_anchor=(0.5, 0.5),
        ncol=2,
        frameon=False
    )
    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved legend:", out_path)


def save_histograms(df: pd.DataFrame, n_obs: int, out_dir: Path):
    df_subset = df[df["num_obs"] == n_obs]
    conf_levels = sorted(df_subset["confounding"].unique())

    # Create subfolder for this sample size
    n_out_dir = out_dir / f"n{n_obs}"
    n_out_dir.mkdir(parents=True, exist_ok=True)

    legend_handles, legend_labels = None, None

    for idx, conf in enumerate(conf_levels):
        sub = df_subset[df_subset["confounding"] == conf]

        fig, ax = plt.subplots(figsize=(6, 4))
        legend_handles, legend_labels = plot_panel(
            ax, sub, n_obs, c=idx
        )

        fig.tight_layout()
        out_path = n_out_dir / f"bias_histogram_conf{conf}_n{n_obs}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)

        print("Saved:", out_path)

    if legend_handles and legend_labels:
        legend_path = n_out_dir / f"bias_hist_legend_n{n_obs}.pdf"
        save_legend_plot(legend_handles, legend_labels, legend_path)


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV not found at {CSV_PATH.resolve()}")

    df = pd.read_csv(CSV_PATH)
    df = prepare_data(df)

    for n_obs in sorted(df["num_obs"].unique()):
        save_histograms(df, n_obs, OUT_DIR)


if __name__ == "__main__":
    main()
