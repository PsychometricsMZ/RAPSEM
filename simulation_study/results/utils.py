from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as style
from typing import List, Dict
from statsmodels.stats.proportion import proportion_confint

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]

style.use('tableau-colorblind10')
prop_cycle = plt.rcParams['axes.prop_cycle']
tableau_colors = prop_cycle.by_key()['color']
_METHOD_COLORS = [tableau_colors[0], tableau_colors[1], tableau_colors[3]]
INTERACTION_PALETTE = tableau_colors[4:8]

print(f"Using these colors: {_METHOD_COLORS}")

PANEL_SIZE = (3,3.75)
LEGEND_SIZE = (7, 0.8)
LABELSIZE = 24
MARKER_LIST = ["o", "^", "X", "*", "v", "P", "D", "s", "H", "8", "p", "d"]
METHOD_MAP = {
    "gestimation": "RAPSEM",
    "standard": "Uncorrected",
    "lavaan": "Lavaan"
}
COLOR_MAP: Dict[str, str] = {
    method: color
    for method, color in zip(METHOD_MAP.values(), _METHOD_COLORS)
} | {
    "Threshold": "black",
    "Gridline": (0.8, 0.8, 0.8)
}
MARKER_MAP: Dict[str, str] = {
    method: marker
    for method, marker in zip(METHOD_MAP.values(), MARKER_LIST)
}

def aggregate_significance(df: pd.DataFrame, group_vars: List[str]) -> pd.DataFrame:
    grouped = df.groupby(group_vars, dropna=False, observed=False).agg(
        significance=("significance", "mean"),
        n_reps=("significance", "count")
    ).reset_index()

    ci = grouped.apply(
        lambda r: proportion_confint(
            count=round(r["significance"] * r["n_reps"]),
            nobs=r["n_reps"],
            method="wilson"
        ), axis=1
    )
    grouped["ci_low"] = ci.apply(lambda x: x[0])
    grouped["ci_high"] = ci.apply(lambda x: x[1])

    return grouped.sort_values(group_vars)

def compute_offsets(n_groups: int, spread: float = 70.0) -> List[float]:
    if n_groups == 1:
        return [0.0]
    step = spread / (n_groups - 1)
    return [-spread / 2 + i * step for i in range(n_groups)]

def apply_minimal_style(ax: plt.Axes):
    ax.grid(False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_facecolor("white")
    ax.figure.set_facecolor("white")
    ax.tick_params(axis="both", which="major", labelsize=LABELSIZE - 4)    
