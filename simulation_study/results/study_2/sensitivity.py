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
    COLOR_MAP,
    apply_minimal_style,
    compute_offsets,
    LABELSIZE,
    PANEL_SIZE,
    LEGEND_SIZE,
    METHOD_MAP,
    INTERACTION_PALETTE,
    MARKER_LIST
)

CSV_PATHS = [
    Path("simulation_study/results/study_2/study_2_results.csv"),
]
OUT_DIR = Path("simulation_study/results/study_2/plots/sensitivity")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_RELIABILITY = 0.90

def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["variable"] == "sensitivity"].copy()
    df = df[df["confounding"] == 0.0].copy()
    df = df[df["true_theta_m"] == 0.0].copy()
    df = df[df["interaction"] != 0.0].copy()
    df = df[df["reliability"].round(2) == TARGET_RELIABILITY].copy()
    df["method"] = df["method"].astype(str).map(METHOD_MAP).fillna(df["method"])
    df = df[df["method"] == "RAPSEM"].copy()
    df["true_theta_m"] = df["true_theta_m"].round(6)
    df["interaction"] = df["interaction"].round(6)
    return df

def compute_inv_kappa(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["individual_sensitivity"] = 1.0 / df["mean"] 
    
    agg = (
        df.groupby(["num_obs", "interaction"], as_index=False)
        .agg(
            inv_kappa_mean=("individual_sensitivity", "mean"),
            inv_kappa_sd=("individual_sensitivity", "std"),
            n=("individual_sensitivity", "count"),
        )
    )
    
    agg["inv_kappa"] = agg["inv_kappa_mean"]
    agg["se_inv_kappa"] = agg["inv_kappa_sd"] / (agg["n"] ** 0.5)
    
    agg["ci_low"] = 1.96 * agg["se_inv_kappa"]
    agg["ci_high"] = 1.96 * agg["se_inv_kappa"]
    return agg


def plot_panel(agg_df: pd.DataFrame, interaction_levels: List, show_ylabel: bool) -> plt.Figure:
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
    ytick_labels = [0.02, 0.04, 0.06, 0.08, 0.1]
    for y in ytick_labels:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=1)
    ax.set_ylim(0, 0.115)
    ax.set_yticks(ytick_labels)
    ax.set_ylabel("Robustness", fontsize=LABELSIZE)
    offsets = compute_offsets(len(interaction_levels))
    for idx, (interaction, g) in enumerate(
        agg_df.groupby("interaction", sort=True)
    ):
        x = g["num_obs"].to_numpy()
        y = g["inv_kappa"].to_numpy()
        yerr = [g["ci_low"].to_numpy(), g["ci_high"].to_numpy()]
        color = INTERACTION_PALETTE[idx % len(INTERACTION_PALETTE)]
        marker = MARKER_LIST[idx % len(MARKER_LIST) + 3]
        x_jittered = x + offsets[idx]

        ax.errorbar(
            x_jittered, y,
            yerr=yerr,
            color=color,
            marker=marker,
            markersize=8,
            capsize=5,
            elinewidth=2.5,
            linestyle="None",
        )

    return fig


def save_panels(df: pd.DataFrame) -> List[str]:
    saved: List[str] = []
    cme_levels = sorted(df["true_theta_m"].unique())
    interaction_levels = sorted(df["interaction"].unique())

    for i, cme in enumerate(cme_levels):
        df_cme = df[df["true_theta_m"] == cme].copy()
        agg_df = compute_inv_kappa(df_cme)
        show_ylabel = (i == 0)
        fig = plot_panel(agg_df, interaction_levels, show_ylabel=show_ylabel)
        out_path = OUT_DIR / f"sensitivity_inv_kappa_cme_{cme}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
        saved.append(str(out_path))

    return saved


def save_legend(interaction_levels: List) -> str:
    if not interaction_levels:
        print("Warning: No interaction levels found. Skipping legend.")
        return ""
    
    handles = [
        Line2D(
            [0], [0],
            linestyle="-",
            linewidth=2.5,
            marker=MARKER_LIST[i % len(MARKER_LIST) + 3],
            markersize=8,
            color=INTERACTION_PALETTE[i % len(INTERACTION_PALETTE)],
            label=f"{b:.3f}"
        )
        for i, b in enumerate(sorted(interaction_levels))
    ]
    
    fig, ax = plt.subplots(figsize=(3, 1.5))
    ax.axis("off")

    ax.legend(
        handles=handles,
        title=r"Interaction $\gamma_{xr}$",
        loc="center",
        bbox_to_anchor=(0.5, 0.5),
        ncol=1,
        frameon=False,
        fontsize=LABELSIZE - 2,
        title_fontsize=LABELSIZE,
    )
    out_path = OUT_DIR / "sensitivity_legend.pdf"
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

    interaction_levels = sorted(df["interaction"].unique())
    saved = save_panels(df)
    saved.append(save_legend(interaction_levels))

    print("Saved:")
    for p in saved:
        print(" -", p)


if __name__ == "__main__":
    main()