#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import warnings
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Plotting
try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    SEABORN_OK = tuple(map(int, sns.__version__.split(".")[:2])) >= (0, 12)
except Exception:
    import matplotlib.pyplot as plt  # fallback only
    sns = None
    SEABORN_OK = False
    warnings.warn("Seaborn not available or < 0.12. Falling back to matplotlib-only plots.")


# -------------------- CLI args --------------------

def parse_args():
    p = argparse.ArgumentParser(description="Recreate the R plotting/analysis in Python (pandas + seaborn).")
    p.add_argument("--csv-1024", required=True, help="Path to 1024px CSV")
    p.add_argument("--csv-576",  required=True, help="Path to 576px CSV")
    p.add_argument("--out-dir",  required=True, help="Directory for outputs (PNGs/CSVs)")
    p.add_argument("--top-n",    type=int, default=10, help="Top-N for ranking tables")
    p.add_argument("--with-errorbars",
                   choices=["none","sd","ci"], default="none",
                   help="Add error bars to line plots (seaborn ≥0.12: 'sd' or 'ci'); 'none' disables.")
    return p.parse_args()


# -------------------- Config --------------------

@dataclass(frozen=True)
class ModeCfg:
    baseline: str
    change_re: str
    label: str

PATTERN_BY_MODE: Dict[str, ModeCfg] = {
    "macro": ModeCfg("dice_unref_vs_gt_macro", r"^dice_unref_vs_k\d+_avmacro$", "Overall (macro A/V)"),
    "A":     ModeCfg("dice_unref_vs_gt_A",     r"^dice_unref_vs_k\d+_A$",       "Arteries"),
    "V":     ModeCfg("dice_unref_vs_gt_V",     r"^dice_unref_vs_k\d+_V$",       "Veins"),
}

@dataclass(frozen=True)
class DiceCfg:
    baseline: str
    refined_re: str
    label: str

DICE_PAT: Dict[str, DiceCfg] = {
    "macro": DiceCfg("dice_unref_vs_gt_macro", r"^dice_k\d+_vs_gt_macro$", "Overall (macro A/V)"),
    "A":     DiceCfg("dice_unref_vs_gt_A",     r"^dice_k\d+_vs_gt_A$",     "Arteries"),
    "V":     DiceCfg("dice_unref_vs_gt_V",     r"^dice_k\d+_vs_gt_V$",     "Veins"),
}

MODES: List[str] = ["macro", "A", "V"]


# -------------------- Helpers --------------------

def add_disease_info(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["disease_info"] = pd.NA
    df.loc[df["image"].str.contains(r"_A\b", regex=True, na=False), "disease_info"] = "AMD"
    df.loc[df["image"].str.contains(r"_N\b", regex=True, na=False), "disease_info"] = "Normal"
    df.loc[df["image"].str.contains(r"_G\b", regex=True, na=False), "disease_info"] = "Glaucoma"
    df.loc[df["image"].str.contains(r"_D\b", regex=True, na=False), "disease_info"] = "Diabetic"
    df["disease_info"] = df["disease_info"].astype("category")
    return df


def add_change_cols_for_mode(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    cfg = PATTERN_BY_MODE[mode]
    df = df.copy()
    kcols = [c for c in df.columns if re.search(cfg.change_re, c)]
    if not kcols:
        raise ValueError(f"No change columns found for mode '{mode}'. Expected pattern: {cfg.change_re}")

    for c in kcols:
        m = re.match(r"^dice_unref_vs_k(\d+)_", c)
        if not m:
            continue
        K = int(m.group(1))
        df[f"change_{mode}_k{K}"] = 1.0 - df[c]
    return df


def to_long_for_mode(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    cfg = PATTERN_BY_MODE[mode]
    base_col = cfg.baseline
    if base_col not in df.columns:
        raise ValueError(f"Baseline column '{base_col}' not found for mode '{mode}'.")
    kcols = [c for c in df.columns if re.match(fr"^change_{mode}_k\d+$", c)]
    use_cols = ["image", base_col, "disease_info"] + kcols
    tmp = df.loc[:, [c for c in use_cols if c in df.columns]].copy()
    tmp = tmp.rename(columns={base_col: "baseline"})
    long_df = tmp.melt(id_vars=["image", "baseline", "disease_info"],
                       value_vars=kcols,
                       var_name="kcol", value_name="change")
    long_df["K"] = long_df["kcol"].str.extract(r"(\d+)").astype(int)
    long_df = long_df[np.isfinite(long_df["change"]) & np.isfinite(long_df["baseline"])]
    return long_df


def per_image_overall_for_mode(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    cfg = PATTERN_BY_MODE[mode]
    kcols = [c for c in df.columns if re.match(fr"^change_{mode}_k\d+$", c)]
    if not kcols:
        raise ValueError(f"No change columns for mode '{mode}'.")
    M = df[kcols].to_numpy(dtype=float)
    M_idx = M.copy()
    M_idx[~np.isfinite(M_idx)] = -np.inf
    idx = np.argmax(M_idx, axis=1)
    K_vals = np.array([int(re.search(r"(\d+)$", c).group(1)) for c in kcols])
    change_max = M[np.arange(M.shape[0]), idx]
    change_mean = np.nanmean(M, axis=1)
    out = pd.DataFrame({
        "image": df["image"].values,
        "disease_info": df.get("disease_info", pd.Series([pd.NA]*len(df))).values,
        "baseline": df[cfg.baseline].values,
        "change_max": change_max,
        "change_mean": change_mean,
        "K_at_max": K_vals[idx]
    })
    out = out.sort_values("change_max", ascending=False, ignore_index=True)
    return out


def dice_vs_gt_long(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    cfg = DICE_PAT[mode]
    base_col = cfg.baseline
    if base_col not in df.columns:
        raise ValueError(f"Missing baseline col: {base_col}")
    kcols = [c for c in df.columns if re.match(cfg.refined_re, c)]
    if not kcols:
        raise ValueError(f"No refined-vs-GT cols matching: {cfg.refined_re}")
    cols = [base_col] + kcols
    tmp = df.loc[:, ["image"] + cols].copy()
    long_df = tmp.melt(id_vars=["image"], value_vars=cols, var_name="kcol", value_name="dice")
    long_df["K"] = np.where(long_df["kcol"] == base_col, 0,
                            long_df["kcol"].str.extract(r"(\d+)").astype(int))
    long_df = long_df[np.isfinite(long_df["dice"])]
    return long_df


# -------------------- Plotting --------------------

def _maybe_seaborn_line(ax, data, x, y, title, xlabel, ylabel, error_mode):
    if sns is not None and SEABORN_OK:
        # seaborn ≥ 0.12: errorbar can be ("sd" | ("pi", 95) | "ci" | None)
        err = None
        if error_mode == "sd":
            err = "sd"
        elif error_mode == "ci":
            err = "ci"
        sns.lineplot(data=data, x=x, y=y, marker="o", errorbar=err, ax=ax)
    else:
        # matplotlib fallback (no CI)
        ax.plot(data[x], data[y], marker="o")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def mean_change_plot(long_df: pd.DataFrame, tag: str, mode_label: str, error_mode: str) -> plt.Figure:
    pdat = long_df.groupby("K", as_index=False)["change"].mean().rename(columns={"change":"mean_change"})
    fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
    _maybe_seaborn_line(ax, pdat, "K", "mean_change",
                        f"{tag}: mean change — {mode_label}",
                        "K iterations", "mean(1 − Dice(unref, refined))", error_mode)
    if sns is None or not SEABORN_OK:
        ax.set_xticks(sorted(pdat["K"].unique()))
    return fig


def scatter_change_vs_baseline(long_df: pd.DataFrame, tag: str, mode_label: str) -> plt.Figure:
    # Facet by K (one row)
    ks = sorted(long_df["K"].unique())
    if sns is not None:
        g = sns.FacetGrid(long_df, col="K", col_order=ks, col_wrap=len(ks), height=3, sharex=True, sharey=True)
        g.map_dataframe(sns.scatterplot, x="baseline", y="change", alpha=0.6, s=16)
        # Add thin OLS line per facet (no CI for consistency with R)
        try:
            g.map_dataframe(sns.regplot, x="baseline", y="change", scatter=False, ci=None, line_kws={"linewidth":0.8, "color":"black"})
        except Exception:
            pass
        g.set_axis_labels("Baseline Dice (unref vs GT)", "Change = 1 − Dice(unref, refined)")
        g.fig.suptitle(f"{tag}: change vs baseline — {mode_label}", y=1.02)
        return g.fig
    else:
        # Simple matplotlib fallback: put all K in one plot, color-coded
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        for K in ks:
            sub = long_df[long_df["K"] == K]
            ax.scatter(sub["baseline"], sub["change"], alpha=0.6, s=16, label=f"K={K}")
        ax.set_title(f"{tag}: change vs baseline — {mode_label}")
        ax.set_xlabel("Baseline Dice (unref vs GT)")
        ax.set_ylabel("Change = 1 − Dice(unref, refined)")
        ax.grid(True, alpha=0.3)
        ax.legend(title="K")
        return fig


def box_by_disease(long_df: pd.DataFrame, tag: str, mode_label: str) -> plt.Figure:
    if sns is not None:
        g = sns.catplot(
            data=long_df, x="disease_info", y="change", kind="box",
            col="K", col_wrap=4, height=3, sharey=True, fliersize=2, linewidth=1
        )
        g.set_axis_labels("Disease", "Change")
        g.fig.suptitle(f"{tag}: change by disease — {mode_label}", y=1.02)
        return g.fig
    else:
        # matplotlib fallback: one box per disease (no facets)
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
        groups = [g["change"].values for _, g in long_df.groupby("disease_info")]
        ax.boxplot(groups, vert=True)
        ax.set_xticks(range(1, len(groups)+1))
        ax.set_xticklabels([str(k) for k in long_df["disease_info"].dropna().unique()])
        ax.set_title(f"{tag}: change by disease — {mode_label}")
        ax.set_xlabel("Disease")
        ax.set_ylabel("Change")
        ax.grid(True, axis="y", alpha=0.3)
        return fig


def plot_mean_dice_progress(long_df: pd.DataFrame, tag: str, mode_label: str, error_mode: str) -> plt.Figure:
    pdat = long_df.groupby("K", as_index=False)["dice"].mean().rename(columns={"dice":"mean_dice"})
    fig, ax = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    _maybe_seaborn_line(ax, pdat, "K", "mean_dice",
                        f"{tag}: mean Dice vs GT by K — {mode_label}",
                        "K iterations (0 = unrefined vs GT)", "Mean Dice (vs GT)", error_mode)
    if sns is None or not SEABORN_OK:
        ax.set_xticks(sorted(pdat["K"].unique()))
    return fig


# -------------------- Correlations --------------------

def cor_table_for_mode(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for K, g in long_df.groupby("K"):
        sub = g[["baseline", "change"]].dropna()
        n = len(sub)
        pearson = sub["baseline"].corr(sub["change"], method="pearson") if n > 1 else np.nan
        spearman = sub["baseline"].corr(sub["change"], method="spearman") if n > 1 else np.nan
        rows.append({"K": K, "n": n, "pearson": pearson, "spearman": spearman})
    out = pd.DataFrame(rows).sort_values("K").reset_index(drop=True)
    return out


# -------------------- Main pipeline --------------------

def save_fig(fig: plt.Figure, path: str, **kwargs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def make_progress_plots(df: pd.DataFrame, tag: str, out_dir: str, error_mode: str):
    for mode in ["macro", "A", "V"]:
        cfg = DICE_PAT[mode]
        ldat = dice_vs_gt_long(df, mode)
        fig = plot_mean_dice_progress(ldat, tag, cfg.label, error_mode)
        fn = os.path.join(out_dir, f"dice_vs_gt_progress_{mode}_{tag}.png")
        save_fig(fig, fn)


def main():
    args = parse_args()
    out_dir = args.out_dir
    top_n = args.top_n
    error_mode = None if args.with_errorbars == "none" else args.with_errorbars

    os.makedirs(out_dir, exist_ok=True)

    # Load & prep
    df1024 = pd.read_csv(args.csv_1024)
    df576  = pd.read_csv(args.csv_576)

    df1024 = add_disease_info(df1024)
    df576  = add_disease_info(df576)

    # For each mode
    for mode in MODES:
        cfg = PATTERN_BY_MODE[mode]
        label = cfg.label

        # Add change columns
        d1 = add_change_cols_for_mode(df1024, mode)
        d2 = add_change_cols_for_mode(df576,  mode)

        # Drop rows with non-finite baseline or zero
        d1 = d1[np.isfinite(d1[cfg.baseline]) & (d1[cfg.baseline] != 0)]
        d2 = d2[np.isfinite(d2[cfg.baseline]) & (d2[cfg.baseline] != 0)]

        # Long format
        long1024 = to_long_for_mode(d1, mode)
        long576  = to_long_for_mode(d2, mode)

        # Mean change vs K
        fig = mean_change_plot(long1024, "1024px", label, error_mode)
        save_fig(fig, os.path.join(out_dir, f"mean_change_vs_k_{mode}_1024.png"))
        fig = mean_change_plot(long576,  "576px",  label, error_mode)
        save_fig(fig, os.path.join(out_dir, f"mean_change_vs_k_{mode}_576.png"))

        # Scatter + correlations
        fig = scatter_change_vs_baseline(long1024, "1024px", label)
        save_fig(fig, os.path.join(out_dir, f"scatter_change_vs_baseline_{mode}_1024.png"),)
        fig = scatter_change_vs_baseline(long576,  "576px",  label)
        save_fig(fig, os.path.join(out_dir, f"scatter_change_vs_baseline_{mode}_576.png"))

        cor1024 = cor_table_for_mode(long1024)
        cor576  = cor_table_for_mode(long576)
        cor1024.to_csv(os.path.join(out_dir, f"correlations_{mode}_1024px.csv"), index=False)
        cor576.to_csv(os.path.join(out_dir, f"correlations_{mode}_576px.csv"), index=False)
        print("\n", mode, " 1024px correlations:\n", cor1024.head(), sep="")
        print("\n", mode, " 576px correlations:\n",  cor576.head(), sep="")

        # Disease-wise boxplots
        fig = box_by_disease(long1024, "1024px", label)
        save_fig(fig, os.path.join(out_dir, f"box_change_by_disease_{mode}_1024px.png"))
        fig = box_by_disease(long576,  "576px",  label)
        save_fig(fig, os.path.join(out_dir, f"box_change_by_disease_{mode}_576px.png"))

        # Per-image overall change (rankings)
        overall1024 = per_image_overall_for_mode(d1, mode)
        overall576  = per_image_overall_for_mode(d2, mode)
        overall1024.to_csv(os.path.join(out_dir, f"overall_change_rank_{mode}_1024.csv"), index=False)
        overall576.to_csv(os.path.join(out_dir, f"overall_change_rank_{mode}_576.csv"), index=False)

        print(f"\nMost changed ({mode}, 1024px):\n", overall1024.head(1))
        print(f"\nMost changed ({mode}, 576px):\n",  overall576.head(1))

        # Top-N lists
        overall1024.head(top_n).to_csv(os.path.join(out_dir, f"top{top_n}_most_changed_{mode}_1024.csv"), index=False)
        overall576.head(top_n).to_csv(os.path.join(out_dir, f"top{top_n}_most_changed_{mode}_576.csv"), index=False)

    # Dice vs GT progression (K0, K4, K5, K6, ...)
    make_progress_plots(df1024, "1024px", out_dir, error_mode)
    make_progress_plots(df576,  "576px",  out_dir, error_mode)

    print(f"\nWrote outputs to: {out_dir}\n")


if __name__ == "__main__":
    main()

