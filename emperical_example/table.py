from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

CSV_PATHS = [Path("emperical_example/results/parameters.csv")]
OUT_DIR = Path("emperical_example/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MODEL_ORDER = ["full", "reduced", "reduced_extreme"]
MODEL_LABELS = {
    "full": "M1 (Full)",
    "reduced": "M2 (Reduced)",
    "reduced_extreme": "M3 (Reduced Extreme)",
}
METHOD_KEYS = ["gestimation", "standard", "lavaan"]
METHOD_LABELS = {
    "gestimation": "G-Estimator",
    "standard": "Standard",
    "lavaan": "\\texttt{lavaan}",
}


def robust_se(df: pd.DataFrame) -> pd.Series:
    """Robust bootstrap standard error per imputation from the 95% percentile
    interval, (q97.5 - q2.5) / (2 * 1.96); unlike the standard deviation of
    the bootstrap replicates, it is not inflated by single extreme replicates.
    For lavaan (Wald intervals) it equals the reported standard error."""
    return (df["upper_pct"] - df["lower_pct"]) / (2 * 1.959964)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[(df["dependent"] == "post_ct") & (df["path"] == "post_ti")].copy()
    df["se"] = robust_se(df)

    def pool(g: pd.DataFrame) -> pd.Series:
        m = len(g)
        est = g["estimate"].mean()

        v_within  = np.mean(g["se"] ** 2)
        v_between = g["estimate"].var(ddof=1) if m > 1 else 0.0
        v_total   = v_within + (1 + 1 / m) * v_between

        se_pooled = np.sqrt(v_total)
        ci_half   = 1.96 * se_pooled

        return pd.Series({
            "estimate": est,
            "se_total": se_pooled,
            "ci_lo":    est - ci_half,
            "ci_hi":    est + ci_half,
            "m_count":  m,
        })

    pooled = (
        df.groupby(["model_version", "method"])
          .apply(pool)
          .reset_index()
    )
    return pooled


def format_cell(estimate: float, ci_lo: float, ci_hi: float) -> str:
    return f"${estimate:.2f}$ $[{ci_lo:.2f},\\,{ci_hi:.2f}]$"


def build_latex_table(df: pd.DataFrame) -> str:
    lookup = df.set_index(["model_version", "method"])

    col_spec = "l" + "c" * len(METHOD_KEYS)
    method_header = " & ".join(
        f"\\textbf{{{METHOD_LABELS[m]}}}" for m in METHOD_KEYS
    )

    rows = []
    for model_key in MODEL_ORDER:
        model_label = MODEL_LABELS[model_key]
        cells = [f"\\textbf{{{model_label}}}"]
        for m_key in METHOD_KEYS:
            try:
                row = lookup.loc[(model_key, m_key)]
                cells.append(format_cell(row["estimate"], row["ci_lo"], row["ci_hi"]))
            except KeyError:
                cells.append("---")
        rows.append(" & ".join(cells) + " \\\\")

    rows_str = "\n        ".join(rows)

    table = rf"""
\begin{{table}}[ht!]
    \centering
    \caption{{Comparison of mediation effects ($\check{{\theta}}_m$)
             across estimation methods. Values represent point estimates
             and 95\% confidence intervals
             $\check{{\theta}} \pm 1.96\,\check{{\sigma}}_{{\text{{pooled}}}}$
             aggregated across imputations via Rubin's rules, using robust
             bootstrap standard errors derived from the percentile
             intervals within each imputation.}}
    \label{{tab:empirical_mediation}}
    \begin{{tabular}}{{{col_spec}}}
    \hline
    \textbf{{Model}} & \multicolumn{{{len(METHOD_KEYS)}}}{{c}}{{\textbf{{Estimation Method}}}} \\
    \cline{{2-{len(METHOD_KEYS)+1}}}
     & {method_header} \\
    \hline
        {rows_str}
    \hline\hline
    \end{{tabular}}
    \begin{{tablenotes}}
        \small
        \item \textit{{Note.}}
        M1\,=\,full model;\;M2\,=\,reduced model;\;M3\,=\,extreme reduced model.
        Estimates and 95\% confidence intervals are shown as
        $\hat{{\theta}}\ [\ell,\,u]$.
    \end{{tablenotes}}
\end{{table}}
""".strip()
    return table


def main():
    if not CSV_PATHS[0].exists():
        raise FileNotFoundError(f"Missing CSV: {CSV_PATHS[0]}")

    df = pd.read_csv(CSV_PATHS[0])
    df = prepare_data(df)

    if df.empty:
        print("Error: No data found for the specified filter.")
        return

    latex = build_latex_table(df)

    out_path = OUT_DIR / "mediation_table.tex"
    out_path.write_text(latex, encoding="utf-8")
    print(f"LaTeX table written to: {out_path}")

if __name__ == "__main__":
    main()
