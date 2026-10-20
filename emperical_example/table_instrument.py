from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

RESULTS_DIR = Path("emperical_example/results")
OUT_DIR     = Path("emperical_example/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CSV_PATTERN = "valid_instrument_results_lambda*.csv"
# boosts shown in the table, in this order (increasing instrument strength)
LAMBDAS = [-0.2, -0.4, -0.6, -0.8]


def load_all(results_dir: Path, pattern: str) -> pd.DataFrame:
    files = sorted(results_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found in {results_dir}"
        )
    chunks = [pd.read_csv(f) for f in files]
    df = pd.concat(chunks, ignore_index=True)
    print(f"Loaded {len(files)} file(s), {len(df)} rows total.")
    return df


def pool_rubins(df: pd.DataFrame, path: str) -> pd.DataFrame:
    df = df[df["path"] == path].copy()
    # robust bootstrap standard error per imputation from the 95% percentile
    # interval (not inflated by single extreme replicates)
    df["se"] = (df["upper_pct"] - df["lower_pct"]) / (2 * 1.959964)

    pooled = (
        df.groupby("lambda")
        .agg(
            estimate      = ("estimate",      "mean"),
            var_within    = ("se",            lambda x: np.mean(np.square(x))),
            var_between   = ("estimate",      lambda x: np.var(x, ddof=1)),
            m_count       = ("estimate",      "count"),
            interaction_f = ("interaction_f", "first"),
        )
        .reset_index()
    )

    pooled["se_total"] = np.sqrt(
        pooled["var_within"] + (1 + 1 / pooled["m_count"]) * pooled["var_between"]
    )
    pooled["ci_lo"] = pooled["estimate"] - 1.96 * pooled["se_total"]
    pooled["ci_hi"] = pooled["estimate"] + 1.96 * pooled["se_total"]

    return pooled.sort_values("lambda").reset_index(drop=True)


def fmt_est(estimate: float, ci_lo: float, ci_hi: float) -> str:
    return rf"${estimate:.2f}$ $[{ci_lo:.2f},\,{ci_hi:.2f}]$"


def build_latex_table(theta: pd.DataFrame, kappa: pd.DataFrame) -> str:
    rows = []
    for lam in LAMBDAS:
        rt = theta[np.isclose(theta["lambda"], lam)].iloc[0]
        rk = kappa[np.isclose(kappa["lambda"], lam)].iloc[0]

        kappa_val = rk["estimate"]
        rc_val    = 1.0 / kappa_val if kappa_val != 0 else float("nan")

        rows.append(
            rf"    ${lam:.1f}$ & ${rt['interaction_f']:.1f}$ & "
            rf"{fmt_est(rt['estimate'], rt['ci_lo'], rt['ci_hi'])} & "
            rf"${kappa_val:.1f}$ & "
            rf"${rc_val:.3f}$ \\"
        )

    table = rf"""
\begin{{table}}[ht!]
    \centering
    \caption{{$\check{{\theta}}_m$ estimates under artificially increased instrument strength: The parameter $\lambda$ represents the artificial boost applied to the interaction between treatment and task interest within the mediator equation, $F_{{\text{{interact}}}}$ denotes the $F$-statistic for all interaction terms, and $1/\kappa$ is the robustness parameter.}}
    \label{{tab:lambda}}
    \begin{{tabular}}{{lcccc}}
    \toprule
    $\lambda$ & $F_{{\text{{interact}}}}$ &
    $\hat{{\theta}}_m\ [95\,\%\ \text{{CI}}]$ & $\hat{{\kappa}}$ &
    $1/\hat{{\kappa}}$ \\
    \midrule
{chr(10).join(rows)}
    \bottomrule
    \end{{tabular}}
\end{{table}}
""".strip()

    return table


def main():
    df_raw = load_all(RESULTS_DIR, CSV_PATTERN)

    theta = pool_rubins(df_raw, "post_ti")
    kappa = pool_rubins(df_raw, "sensitivity")

    print("\nPooled theta (post_ti):")
    print(theta[["lambda", "interaction_f", "estimate", "ci_lo", "ci_hi"]].to_string(index=False, float_format="{:.4f}".format))
    print("\nPooled kappa (sensitivity):")
    print(kappa[["lambda", "estimate", "ci_lo", "ci_hi"]].to_string(index=False, float_format="{:.4f}".format))

    latex    = build_latex_table(theta, kappa)
    out_path = OUT_DIR / "lambda_sensitivity_table.tex"
    out_path.write_text(latex, encoding="utf-8")
    print(f"\nTable written to: {out_path}")

if __name__ == "__main__":
    main()
