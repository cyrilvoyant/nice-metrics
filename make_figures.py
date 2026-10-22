"""Figures of the NICE^k demonstration, drawn from the CSV files written by demo_nice.m.

    python make_figures.py

Writes figures/*.png and figures/*.pdf. Requires numpy, pandas and matplotlib.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "axes.titlesize": 11.5, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#dfe5ec",
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False, "savefig.dpi": 200,
})
C = {"observed": "#1b2a3a", "clear_sky": "#9aa8b6", "P": "#d1495b", "SP": "#2e86ab", "AR": "#3c9d5d",
     "Pcycle": "#e0a03a"}
LABEL = {"P": "Persistence (reference)", "SP": "Smart persistence", "AR": "AR(1) on clear-sky residual",
         "Pcycle": "Cyclic persistence (same hour, previous day)"}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# 1. Signal and forecasts over three days
s = pd.read_csv(RES / "demo_series.csv")
hours = s["t"] - s["t"].iloc[0]
fig, ax = plt.subplots(figsize=(9, 3.6))
ax.fill_between(hours, s["clear_sky"], color=C["clear_sky"], alpha=0.18, lw=0, label="Clear-sky envelope")
ax.plot(hours, s["observed"], color=C["observed"], lw=2.0, label="Observed")
for k in ["P", "SP", "AR"]:
    ax.plot(hours, s[k], color=C[k], lw=1.3, ls="--" if k == "P" else "-", label=LABEL[k])
ax.set_xlabel("Hour")
ax.set_ylabel("Irradiance-like signal (W/m²)")
ax.set_xticks(range(0, 73, 12))
ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9)
save(fig, "fig1_signal_forecasts")

# 2. NICE^k of each model, persistence = 1
sc = pd.read_csv(RES / "demo_scores.csv")
models = ["SP", "AR", "Pcycle"]
metrics = ["NICE1", "NICE2", "NICE3", "NICE_Sigma"]
mlabel = ["NICE$^1$ (MAE)", "NICE$^2$ (RMSE)", "NICE$^3$ (RMCE)", "NICE$^\\Sigma$ (mean)"]
x = np.arange(len(metrics))
w = 0.26
fig, ax = plt.subplots(figsize=(8.2, 3.8))
for i, m in enumerate(models):
    vals = sc.set_index("model").loc[m, metrics].to_numpy(float)
    bars = ax.bar(x + (i - 1) * w, vals, w, color=C[m], label=LABEL[m])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
ax.axhline(1, color=C["P"], ls="--", lw=1.4)
ax.text(-0.5, 1.02, "persistence = 1", color=C["P"], ha="left", va="bottom", fontsize=9)
ax.set_xticks(x, mlabel)
ax.set_ylabel("NICE (lower is better)")
ax.set_ylim(0, 1.45)
ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9)
save(fig, "fig2_nice_scores")

# 3. Theory versus measurement for the persistence errors (200 realisations)
mc = pd.read_csv(RES / "montecarlo_persistence.csv")
fig, axes = plt.subplots(1, 3, figsize=(11, 3.9))
panels = [("MAE_P", "MAE_P_theory", None, "MAE of persistence (k = 1)"),
          ("RMSE_P", "RMSE_P_theory", None, "RMSE of persistence (k = 2)"),
          ("RMCE_P", "RMCE_P_theory", "RMCE_P_article_constant", "RMCE of persistence (k = 3)")]
for ax, (emp, th, th_art, title) in zip(axes, panels):
    lo, hi = 0, 1.05 * max(mc[emp].max(), mc[th].max(), mc[th_art].max() if th_art else 0)
    ax.plot([lo, hi], [lo, hi], color="#46586a", lw=1, ls=":", label="theory = measured")
    ax.scatter(mc[th], mc[emp], s=12, color=C["SP"], alpha=0.75,
               label=f"Gaussian theory, median ratio {np.median(mc[emp] / mc[th]):.2f}")
    if th_art:
        ax.scatter(mc[th_art], mc[emp], s=12, color=C["P"], alpha=0.6, marker="x",
                   label=f"constant as printed, median ratio {np.median(mc[emp] / mc[th_art]):.2f}")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title(title)
    ax.set_xlabel("Theory from σ and ρ(1) (W/m²)")
    ax.legend(loc="upper left", fontsize=8)
axes[0].set_ylabel("Measured on the series (W/m²)")
save(fig, "fig3_theory_vs_measured")
print("figures written to", FIG)
