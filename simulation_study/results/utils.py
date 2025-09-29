from __future__ import annotations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from typing import List, Dict

PANEL_SIZE = (4,5)
LABELSIZE = 20
COLOR_MAP: Dict[str, str] = {
    "0.4": "#0f0bde",
    "0.5": "#940bde",
    "0.667": "#ce0ca4",
    "0.8": "#076d07",
    "ML Estimator": "#0bc5de",
    "G-Estimator": "#076d07",
    "Threshold": "#847e7e",
    "Gridline": (0.8, 0.8, 0.8)
}
LINETYPE_MAP: Dict[str, tuple] = {
    "dotted": (0, (1, 2)), 
    "dotdash": (0, (3, 2, 9, 2)),
    "dashed": (0, (6, 3)),
    "solid": (0, (None, None)),
    "longdash": (0, (12, 3)),
    "twodash": (0, (6, 3, 2, 3)),
}
MARKERS = ["o", "s", "D", "^", "v", "<", ">"]

def aggregate_significance(df: pd.DataFrame, group_vars: List[str]) -> pd.DataFrame:
    grouped = df.groupby(group_vars, dropna=False, observed=False)["significance"].mean().reset_index()
    return grouped.sort_values(group_vars)


def apply_minimal_style(ax: plt.Axes):
    ax.grid(False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_facecolor("white")
    ax.figure.set_facecolor("white")
    ax.tick_params(axis="both", which="major", labelsize=LABELSIZE - 4)


def apply_yaxes_significance_style(ax: plt.Axes):
    ax.axhline(0.8, color=COLOR_MAP.get("Threshold", "black"), linewidth=0.8)
    for y in [0.2, 0.4, 0.6, 1.0]:
        ax.axhline(y, color=COLOR_MAP.get("Gridline", "black"), linewidth=0.6)
    ax.set_ylim(-0.01, 1.01)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
