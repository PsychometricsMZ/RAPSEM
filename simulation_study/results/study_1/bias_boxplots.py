from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))
from simulation_study.results.utils import (
    apply_minimal_style,
    COLOR_MAP,
    LEGEND_SIZE,
    METHOD_MAP,
    LABELSIZE
    )

NEW_LABELSIZE = LABELSIZE - 6

CSV_PATHS = [
    Path("simulation_study/results/study_1/study_1_results.csv"),
]
OUT_DIR = Path("simulation_study/results/study_1/plots/bias_boxplots")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["variable"] == "med"]
    # confounding scenarios of study 1 (no mediator-outcome confounding)
    df = df[(df["true_theta_m"] == 0.0) & (df["conf_med"] == 0.0)].copy()
    df["confounding"] = df["confounding"].round(3)
    df["method"] = df["method"].astype(str).map(METHOD_MAP).fillna(df["method"])
    return df


def plot_panel(ax, sub_df: pd.DataFrame, c: float):
    apply_minimal_style(ax)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax.tick_params(
        axis='both', 
        which='major', 
        length=6, 
        width=1.5, 
        labelsize=NEW_LABELSIZE - 2
    )
    ax.set_ylim(-0.55, 0.72)
    ytick_labels = [-0.5, -0.25, 0.0, 0.25, 0.5]
    if c == 0.1:
         ax.set_yticks(ytick_labels)
         ax.set_ylabel("Samples", fontsize=NEW_LABELSIZE)
    else:
        ax.set_yticks([])
        ax.set_ylabel("")   
    for y in ytick_labels:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=1)
    ax.axhline(0, color=COLOR_MAP.get("Threshold", "black"), linewidth=1.5)
    method_order = list(METHOD_MAP.values())
    methods = [m for m in method_order if m in sub_df["method"].unique()]
    data = [sub_df[sub_df["method"] == m]["mean"].values for m in methods]
    colors = [COLOR_MAP.get(m, "grey") for m in methods]
    bp = ax.boxplot(
        data,
        patch_artist=True,
        widths=0.5,
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5),
        flierprops=dict(marker="o", markersize=3, linestyle="none", alpha=0.4),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    for flier, color in zip(bp["fliers"], colors):
        flier.set_markerfacecolor(color)
        flier.set_markeredgecolor(color)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_MAP.get(m, "grey"), alpha=0.6)
        for m in methods
    ]
    return handles, methods


def save_legend_plot(handles, labels, out_path: Path):
    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    ax.axis("off")
    ax.legend(
        handles, labels,
        title="Method",
        loc="center",
        fontsize=LABELSIZE - 2,
        title_fontsize=LABELSIZE,
        bbox_to_anchor=(0.5, 0.5),
        ncol=3,
        frameon=False
    )
    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved legend:", out_path)


def save_boxplots(df: pd.DataFrame, n_obs: int, out_dir: Path):
    df_subset = df[df["num_obs"] == n_obs]
    conf_levels = sorted(df_subset["confounding"].unique())
    n_out_dir = out_dir / f"n{n_obs}"
    n_out_dir.mkdir(parents=True, exist_ok=True)
    legend_handles, legend_labels = None, None

    for idx, conf in enumerate(conf_levels):
        sub = df_subset[df_subset["confounding"] == conf]
        fig, ax = plt.subplots(figsize=(2, 2))
        legend_handles, legend_labels = plot_panel(ax, sub, c=conf)
        out_path = n_out_dir / f"bias_boxplot_conf_{conf}_n{n_obs}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print("Saved:", out_path)
        if idx == 1 and legend_handles and legend_labels:
            legend_path = out_dir / f"bias_boxplot_legend.pdf"
            save_legend_plot(legend_handles, legend_labels, legend_path)


def main():
    missing = [p for p in CSV_PATHS if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing CSV files: {missing}")
    dfs = [pd.read_csv(p) for p in CSV_PATHS]
    df = pd.concat(dfs, ignore_index=True)
    df = prepare_data(df)
    for n_obs in sorted(df["num_obs"].unique()):
        save_boxplots(df, n_obs, OUT_DIR)


if __name__ == "__main__":
    main()