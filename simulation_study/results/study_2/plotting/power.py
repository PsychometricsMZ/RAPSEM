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
    LINETYPE_MAP,
    LABELSIZE,
    PANEL_SIZE,
    MARKERS,
)

CSV_PATH = Path("study_2_results.csv")
OUT_DIR = Path("results/study_2/plotting/power")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CUSTOM_ORDER: List[str] = ["dotted", "dotdash", "dashed", "solid", "longdash", "twodash"]
LEGEND_SIZE = (8, 1)
ALLOWED_RELIABILITIES = {0.4, 0.5, 0.667, 0.8}


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["method"] = df["method"].apply(lambda m: "ML Estimator" if str(m) == "standard" else "G-Estimator")
    df = df[df["method"] == "G-Estimator"].copy()
    df = df[df["num_obs"] != 100].copy()
    df = df[df["beta_ym"] != 0].copy()
    df = df[df["reliability"].isin(ALLOWED_RELIABILITIES)].copy()
    df["beta_mxr"] = df["beta_mxr"].astype(str).astype("category")
    return df


def apply_axes_style(ax: plt.Axes, show_ylabel: bool):
    apply_minimal_style(ax)
    apply_yaxes_significance_style(ax)
    ax.set_ylabel("Proportion Significant" if show_ylabel else " ", fontsize=LABELSIZE)
    ax.set_xticks([250, 500, 750, 1000])  # removed 100
    ax.set_xlabel(r"Sample size $N$", fontsize=LABELSIZE)


def plot_panel(sig_df: pd.DataFrame, beta_mxr_levels: List[str], show_ylabel: bool) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    apply_axes_style(ax, show_ylabel=show_ylabel)

    for bmxr, g in sig_df.groupby("beta_mxr", sort=False, observed=False):
        x = g["num_obs"].to_numpy()
        y = g["significance"].to_numpy()

        idx = beta_mxr_levels.index(str(bmxr)) if str(bmxr) in beta_mxr_levels else 0
        style_name = CUSTOM_ORDER[idx % len(CUSTOM_ORDER)]
        marker = MARKERS[idx % len(MARKERS)]
        dash = LINETYPE_MAP.get(style_name, (0, (None, None)))

        line, = ax.plot(
            x,
            y,
            linewidth=1.2,
            linestyle="-",
            marker=marker,
            markersize=6,
            color="black",
        )
        if dash[1] is not None:
            line.set_dashes(dash[1])

    return fig


def save_panels(df: pd.DataFrame) -> List[str]:
    saved: List[str] = []
    beta_vals = df["beta_ym"].dropna().unique()
    try:
        beta_vals_sorted = sorted(beta_vals, key=lambda v: float(v))
    except Exception:
        beta_vals_sorted = sorted(beta_vals, key=lambda v: str(v))
    beta_mxr_levels = list(pd.Categorical(df["beta_mxr"]).categories.astype(str))

    for beta in beta_vals_sorted:
        df_beta = df[df["beta_ym"] == beta].copy()
        out_dir = OUT_DIR / f"{beta}"
        out_dir.mkdir(parents=True, exist_ok=True)

        reliabilities = sorted(df_beta["reliability"].dropna().unique())
        for i, rel in enumerate(reliabilities, start=1):
            df_rel = df_beta[df_beta["reliability"] == rel].copy()
            sig_df = aggregate_significance(df_rel, ["num_obs", "beta_mxr"])
            fig = plot_panel(sig_df, beta_mxr_levels, show_ylabel=(i == 1))
            out_path = out_dir / f"power_cme_{beta}rel_{rel}.pdf"
            fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
            plt.close(fig)
            saved.append(str(out_path))

    return saved


def make_beta_handles(beta_levels: List[str]) -> List[Line2D]:
    handles = []
    for i, b in enumerate(beta_levels):
        style_name = CUSTOM_ORDER[i % len(CUSTOM_ORDER)]
        marker = MARKERS[i % len(MARKERS)]
        dash = LINETYPE_MAP.get(style_name, (0, (None, None)))
        h = Line2D(
            [0], [0],
            linestyle="-",
            linewidth=1.2,
            color="black",
            marker=marker,
            markersize=6,
            label=str(b),
        )
        if dash[1] is not None:
            h.set_dashes(dash[1])
        handles.append(h)
    return handles


def save_legend(df: pd.DataFrame) -> str:
    first_beta = df["beta_ym"].dropna().iloc[0]
    df_example = df[df["beta_ym"] == first_beta].copy()

    b_levels = list(pd.Categorical(df_example["beta_mxr"].astype(str)).categories)
    b_handles = make_beta_handles(b_levels)

    fig, ax = plt.subplots(figsize=LEGEND_SIZE)
    ax.axis("off")

    if b_handles:
        ax.legend(
            handles=b_handles,
            title=r"Interaction $\gamma_{xr}$",
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            ncol=min(len(b_handles), 6),
            frameon=False,
        )

    out_path = OUT_DIR / "power_legend.pdf"
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(out_path)


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"CSV not found at {CSV_PATH.resolve()}")
    df = pd.read_csv(CSV_PATH)
    df = prepare_data(df)
    saved = save_panels(df)
    legend_path = save_legend(df)
    print("Saved:")
    for p in saved:
        print(" -", p)
    print(" -", legend_path)


if __name__ == "__main__":
    main()
