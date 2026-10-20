from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
import numpy as np

# Setup project root and imports
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

from simulation_study.results.utils import (
    apply_minimal_style,
    compute_offsets,
    COLOR_MAP,
    MARKER_MAP,
    METHOD_MAP,
    LABELSIZE,
    PANEL_SIZE
)

CSV_PATHS = [Path("emperical_example/results/parameters.csv")]
OUT_DIR = Path("emperical_example/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Study Design Definitions
MODEL_ORDER = ["full", "reduced", "reduced_extreme"]
MODEL_LABELS = {"full": "1", "reduced": "2", "reduced_extreme": "3"}
METHOD_KEYS = ["gestimation", "standard", "lavaan"]
LEGEND_SIZE = (8, 1)

def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # Path filtering based on empirical experiment columns
    df = df[(df["dependent"] == "post_ct") & (df["path"] == "post_ti")].copy()
    # robust bootstrap standard error per imputation from the 95% percentile
    # interval (not inflated by single extreme replicates; equals the Wald
    # standard error for lavaan)
    df["se"] = (df["upper_pct"] - df["lower_pct"]) / (2 * 1.959964)

    # Rubin's Rules for pooling across multiple imputations
    pooled = df.groupby(['model_version', 'method']).agg(
        estimate=('estimate', 'mean'),
        var_within=('se', lambda x: np.mean(np.square(x))),
        var_between=('estimate', lambda x: np.var(x, ddof=1)),
        m_count=('estimate', 'count')
    ).reset_index()
    
    # Pooled Standard Error Calculation
    pooled['se_total'] = np.sqrt(
        pooled['var_within'] + (1 + 1/pooled['m_count']) * pooled['var_between']
    )
    pooled['ci_half'] = 1.96 * pooled['se_total']
    return pooled

def plot_estimates(df: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    apply_minimal_style(ax)
    
    # Logic: X-axis is Model Complexity (1, 2, 3)
    x_centers = np.arange(len(MODEL_ORDER))
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks(x_centers)
    ax.set_xticklabels([f'M{i+1}' for i in range(len(MODEL_ORDER))], fontsize=LABELSIZE)
    ax.set_xlabel("Model Complexity", fontsize=LABELSIZE)
    ax.tick_params(
        axis='x', 
        which='major', 
        bottom=False,
        top=False, 
        labelsize=LABELSIZE - 4
    )
    ax.tick_params(
        axis='y', 
        which='major', 
        length=6, 
        width=1.5,
    )
    ytick_labels = [-0.75, -0.5, -0.25, 0, 0.25, 0.5]
    for y in ytick_labels:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=1)
    ax.axhline(0, color=COLOR_MAP.get("Threshold", "black"), linewidth=1.5)
    ax.set_ylabel("$\\theta_m$", fontsize=LABELSIZE)

    # Calculate Jitter
    # NOTE: compute_offsets default spread=30.0 is for N-scales. 
    # For categorical 0, 1, 2, we use a small spread (e.g., 0.3)
    methods_present = [m for m in METHOD_KEYS if m in df["method"].unique()]
    offsets = compute_offsets(len(methods_present), spread=0.3)
    
    for idx, m_key in enumerate(methods_present):
        # Map raw CSV key ('gestimation') to display name ('G-Estimator')
        m_display = METHOD_MAP.get(m_key, m_key)
        
        # Select data and ensure it follows MODEL_ORDER
        df_sub = df[df['method'] == m_key].set_index('model_version').reindex(MODEL_ORDER)
        
        # X-position = Model Center + Jitter Offset
        x_pos = x_centers + offsets[idx]
        
        ax.errorbar(
            x_pos, 
            df_sub['estimate'].to_numpy(),
            yerr=df_sub['ci_half'].to_numpy(),
            label=m_display,
            color=COLOR_MAP.get(m_display, "black"),
            marker=MARKER_MAP.get(m_display, "o"),
            markersize=8,
            linestyle="None",
            capsize=5,
            elinewidth=2.5,
            zorder=3
        )
    
    ax.set_xlim(-0.5, len(MODEL_ORDER) - 0.5)
    return fig

def save_legend(df: pd.DataFrame) -> str:
    methods_present = [m for m in METHOD_KEYS if m in df["method"].unique()]
    
    method_handles = []
    for m_key in methods_present:
        m_display = METHOD_MAP.get(m_key, m_key)
        handle = Line2D([0], [0], linestyle="-", linewidth=2.5, 
                        marker=MARKER_MAP.get(m_display, "o"), 
                        color=COLOR_MAP.get(m_display, "black"), 
                        label=m_display, markersize=8)
        method_handles.append(handle)
    
    fig, ax = plt.subplots(figsize=(5,1.5))
    ax.axis("off")
    ax.legend(
        handles=method_handles,
        loc="center",
        title="Method",
        fontsize=LABELSIZE - 4,
        title_fontsize=LABELSIZE,
        ncol=2,
        frameon=False,
    )
    
    out_path = OUT_DIR / "mediation_comparison_legend.pdf"
    fig.savefig(out_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(out_path)

def main():
    if not CSV_PATHS[0].exists():
        raise FileNotFoundError(f"Missing CSV: {CSV_PATHS[0]}")
    
    df = pd.read_csv(CSV_PATHS[0])
    df = prepare_data(df)
    
    if df.empty:
        print("Error: No data found for the specified paths.")
        return

    # 1. Main Plot
    fig = plot_estimates(df)
    out_path = OUT_DIR / "mediation_comparison.pdf"
    fig.savefig(out_path, format='pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    # 2. Legend
    legend_path = save_legend(df)
    
    print(f"Empirical plots generated successfully:\n - {out_path}\n - {legend_path}")

if __name__ == "__main__":
    main()