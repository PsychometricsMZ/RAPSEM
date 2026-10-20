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
    aggregate_significance,
    apply_minimal_style,
    compute_offsets,
    LABELSIZE,
    PANEL_SIZE,
    LEGEND_SIZE,
    METHOD_MAP,
    INTERACTION_PALETTE,
    MARKER_LIST,
    COLOR_MAP
)

CSV_PATHS = [
    Path("simulation_study/results/study_2/study_2_results.csv"),
]
OUT_DIR = Path("simulation_study/results/study_2/plots/power")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["variable"] == "med"].copy()
    df = df[(df["confounding"] == 0.0)].copy()
    df = df[(df["interaction"] != 0.0)].copy()
    df = df[(df["true_theta_m"] != 0.0)].copy()
    df["interaction"] = df["interaction"].astype(str).astype("category")
    df["reliability"] = df["reliability"].round(3)
    df["true_theta_m"] = df["true_theta_m"].round(3)
    df["method"] = df["method"].astype(str).map(METHOD_MAP).fillna(df["method"])
    return df


def apply_axes_style(ax: plt.Axes, show_ylabel: bool):
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
    ytick_labels = [0.2, 0.4, 0.6, 0.8, 1.0]
    for y in ytick_labels:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=1)
    ax.axhline(0.8, color=COLOR_MAP.get("Threshold", "black"), linewidth=1.5)
    ax.set_ylim(0, 1.02)
    if show_ylabel:
         ax.set_yticks(ytick_labels)
         ax.set_ylabel("Proportion Significant", fontsize=LABELSIZE)
    else:
        ax.set_yticks([])
        ax.set_ylabel("")


def plot_panel(sig_df: pd.DataFrame, interaction_levels: List[str], show_ylabel: bool) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    apply_axes_style(ax, show_ylabel=show_ylabel)
    n_groups = len(interaction_levels)
    offsets = compute_offsets(n_groups)
    for idx, (bmxr, g) in enumerate(sig_df.groupby("interaction", sort=False, observed=False)):
        color = INTERACTION_PALETTE[idx % len(INTERACTION_PALETTE)]
        marker = MARKER_LIST[idx % len(MARKER_LIST) + 3]
        x = g["num_obs"].to_numpy()
        y = g["significance"].to_numpy()
        yerr = [
            (g["significance"] - g["ci_low"]).clip(lower=0).to_numpy(),
            (g["ci_high"] - g["significance"]).clip(lower=0).to_numpy()
        ]
        x_jittered = x + offsets[idx]
        ax.errorbar(
            x_jittered, y, yerr=yerr,
            color=color,
            marker=marker,
            markersize=8,
            linestyle="None",
            capsize=5,
            elinewidth=2.5
        )

    return fig

def save_panels(df: pd.DataFrame) -> List[str]:
    saved: List[str] = []
    beta_vals = df["true_theta_m"].dropna().unique()
    try:
        beta_vals_sorted = sorted(beta_vals, key=lambda v: float(v))
    except Exception:
        beta_vals_sorted = sorted(beta_vals, key=lambda v: str(v))
    interaction_levels = list(pd.Categorical(df["interaction"]).categories.astype(str))
    for i, beta in enumerate(beta_vals_sorted, start=1):
        df_beta = df[df["true_theta_m"] == beta].copy()
        reliabilities = sorted(df_beta["reliability"].dropna().unique())
        for rel in reliabilities:
            df_rel = df_beta[df_beta["reliability"] == rel].copy()
            sig_df = aggregate_significance(df_rel, ["num_obs", "interaction"])
            fig = plot_panel(sig_df, interaction_levels, show_ylabel=(i == 1))
            out_path = OUT_DIR / f"power_cme_{beta}_rel_{rel}.pdf"
            fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
            plt.close(fig)
            saved.append(str(out_path))
    return saved


def make_legend_handles(interaction_levels: List[str]) -> List[Line2D]:
    handles = []
    for i, level in enumerate(interaction_levels):
        color = INTERACTION_PALETTE[i % len(INTERACTION_PALETTE)]
        marker = MARKER_LIST[i % len(MARKER_LIST) + 3]
        h = Line2D(
            [0], [0],
            linestyle="-",
            linewidth=2.5,
            marker=marker,
            markersize=8,
            color=color,
            label=str(level),
        )
        handles.append(h)
    return handles


def save_legend(df: pd.DataFrame) -> str:
    interaction_levels = [f"{float(x):.3f}" for x in sorted(df["interaction"].unique())]
    handles = make_legend_handles(interaction_levels)
    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    ax.axis("off")
    if handles:
        ax.legend(
            handles=handles,
            title=r"Interaction $\gamma_{xr}$",
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            ncol=min(len(handles), 6),
            frameon=False,
            fontsize=LABELSIZE - 2,
            title_fontsize=LABELSIZE,
        )
    out_path = OUT_DIR / "power_legend.pdf"
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
    saved = save_panels(df)
    legend_path = save_legend(df)
    print("Saved:")
    for p in saved:
        print(" -", p)
    print(" -", legend_path)


if __name__ == "__main__":
    main()