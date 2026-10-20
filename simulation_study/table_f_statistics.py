def calculate_f_stat(n, r2_int, r2_main=0.35, m=2, k=5):
    """
    Calculates the F-statistic for a subset of predictors.
    n: sample size
    r2_int: partial R2 of the interaction block
    r2_main: combined R2 of treatment and covariates (default 0.35)
    m: degrees of freedom for the numerator (interaction terms)
    k: total number of predictors in full model
    """
    r2_total = r2_main + r2_int
    df1 = m
    df2 = n - k - 1
    
    numerator = r2_int / df1
    denominator = (1 - r2_total) / df2
    return numerator / denominator

def main():
    sample_sizes = [250, 500, 750, 1000]
    r2_increments = [0.01, 0.06, 0.10, 0.14]
    gamma_labels = ["0.062", "0.152", "0.196", "0.232"]

    # Start LaTeX String
    latex_output = [
        r"\begin{table}[h!]",
        r"    \centering",
        r"    \caption{Theoretical $F$-statistics for the joint test of covariate--treatment interactions ($m=2$) across simulated conditions.}",
        r"    \label{tab:f_stats_simulation}",
        r"    \begin{tabular}{lcccc}",
        r"    \hline",
        r"    \textbf{Sample Size ($n$)} & \multicolumn{4}{c}{\textbf{Interaction Effect Size ($\boldsymbol{\gamma}_{xr}$)}} \\",
        r"    \cline{2-5}",
        f"     & \\textbf{{{gamma_labels[0]}}} & \\textbf{{{gamma_labels[1]}}} & \\textbf{{{gamma_labels[2]}}} & \\textbf{{{gamma_labels[3]}}} \\\\",
        r"    \hline"
    ]

    for n in sample_sizes:
        row_values = [str(n)]
        for r2 in r2_increments:
            f_val = calculate_f_stat(n, r2)
            row_values.append(f"{f_val:.2f}")
        
        latex_output.append(f"        {'  & '.join(row_values)} \\\\")
    latex_output.append(r"    \hline\hline")
    latex_output.append(r"    \end{tabular}")
    latex_output.append(r"\end{table}")
    print("\n".join(latex_output))

if __name__ == "__main__":
    main()