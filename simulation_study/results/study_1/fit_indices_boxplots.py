from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))
from simulation_study.results.utils import (
    apply_minimal_style,
    LABELSIZE,
    COLOR_MAP  # Assuming this exists in your utils
)

OUT_DIR = Path("simulation_study/results/study_1/plots/fit_indices_boxplots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

Y_LIMITS = {
    "CFI": (0.985, 1.001),
    "TLI": (0.98, 1.01),
    "RMSEA": (-0.01, 0.065),
    "SRMR": (0.0, 0.045),
}
fit_indices = list(Y_LIMITS.keys())

CSV_PATHS = [
    Path("simulation_study/results/study_1/study_1_fit_indices.csv"),
]

def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    # confounding scenarios of study 1 (no mediator-outcome confounding)
    df = df[df["conf_med"] == 0.0].copy()
    df["confounding"] = df["confounding"].round(3)
    return df

def plot_fit_index_by_conf(df, fit_index, base_out_dir=OUT_DIR):
    out_dir = base_out_dir / fit_index
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_sizes = sorted(df["num_obs"].unique())
    conf_levels = sorted(df["confounding"].unique())
    
    # Define color for lavaan
    lavaan_color = COLOR_MAP.get("lavaan", "grey")

    for n_obs in sample_sizes:
        sub = df[df["num_obs"] == n_obs]
        if sub.empty:
            continue
            
        plt.figure(figsize=(3, 2))
        ax = plt.gca()
        apply_minimal_style(ax)
        
        # Prepare data for Matplotlib ax.boxplot (list of arrays)
        data_to_plot = [sub[sub["confounding"] == c][fit_index].dropna().values for c in conf_levels]
        
        # Hide standard X-axis ticks initially to match your style
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

        # Plotting with the requested style
        bp = ax.boxplot(
            data_to_plot,
            patch_artist=True,
            widths=0.5,
            medianprops=dict(color="black", linewidth=2),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            flierprops=dict(marker="o", markersize=3, linestyle="none", alpha=0.4),
        )

        # Applying colors and alpha
        for patch in bp["boxes"]:
            patch.set_facecolor(lavaan_color)
            patch.set_alpha(0.6)
            
        for flier in bp["fliers"]:
            flier.set_markerfacecolor(lavaan_color)
            flier.set_markeredgecolor(lavaan_color)

        # Formatting Y-axis
        ax.tick_params(
            axis='both', 
            which='major', 
            length=6, 
            width=1.5, 
            labelsize=LABELSIZE - 4
        )
        
        if n_obs == sample_sizes[0]:
            ax.set_ylabel(fit_index, fontsize=LABELSIZE)
            if fit_index in ["CFI", "TLI"]:
                ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0))
        else:
            ax.set_yticklabels([])

        # X-axis Labels
        ax.set_xticks(range(1, len(conf_levels) + 1))
        ax.set_xticklabels([f"{c:.3f}".replace('.', ',') for c in conf_levels], rotation=90)
        ax.set_xlabel(r"Confounding", fontsize=LABELSIZE)
        ax.tick_params(
            axis='both', 
            which='major', 
            length=6, 
            width=1.5, 
            labelsize=LABELSIZE - 4
        )
        
        if fit_index in Y_LIMITS:
            ax.set_ylim(Y_LIMITS[fit_index])

        out_path = out_dir / f"fit_indices_n{n_obs}.pdf"
        plt.savefig(out_path, bbox_inches="tight", facecolor="white")
        plt.close()

        print("Saved:", out_path)

def main():
    missing = [p for p in CSV_PATHS if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing CSV files: {missing}")
    dfs = [pd.read_csv(p) for p in CSV_PATHS]
    df = pd.concat(dfs, ignore_index=True)
    df = prepare_data(df)

    for fit_index in fit_indices:
        plot_fit_index_by_conf(df, fit_index)

if __name__ == "__main__":
    main()