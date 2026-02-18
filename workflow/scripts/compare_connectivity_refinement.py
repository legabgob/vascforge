#!/usr/bin/env python3
"""
Compare unrefined vs refined connectivity metrics.
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path


# ============================================================================
# Helper functions (pure, no snakemake dependency)
# ============================================================================

def calculate_effect_size(before, after):
    pooled_std = np.sqrt((before.std()**2 + after.std()**2) / 2)
    if pooled_std == 0:
        return 0.0
    return float((after.mean() - before.mean()) / pooled_std)


def perform_statistical_tests(before, after, metric_name):
    t_stat, t_pvalue = stats.ttest_rel(before, after)

    try:
        w_stat, w_pvalue = stats.wilcoxon(before, after)
        wilcoxon = {'statistic': float(w_stat), 'p_value': float(w_pvalue), 'significant': bool(w_pvalue < 0.05)}
    except Exception as e:
        wilcoxon = {'error': str(e)}

    effect_size = calculate_effect_size(before, after)

    return {
        'metric': metric_name,
        'n_samples': int(len(before)),
        'before_mean': float(before.mean()),
        'before_std': float(before.std()),
        'after_mean': float(after.mean()),
        'after_std': float(after.std()),
        'mean_change': float(after.mean() - before.mean()),
        'mean_change_pct': float((after.mean() - before.mean()) / before.mean() * 100) if before.mean() != 0 else 0.0,
        'paired_ttest': {'t_statistic': float(t_stat), 'p_value': float(t_pvalue), 'significant': bool(t_pvalue < 0.05)},
        'wilcoxon': wilcoxon,
        'effect_size': {
            'cohens_d': effect_size,
            'interpretation': (
                'negligible' if abs(effect_size) < 0.2 else
                'small'      if abs(effect_size) < 0.5 else
                'medium'     if abs(effect_size) < 0.8 else
                'large'
            )
        },
    }


def create_comparison_plot(df, metrics, output_path):
    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(6 * n_metrics, 5))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        before = df[f'{metric}_unrefined']
        after  = df[f'{metric}_refined']

        for i in range(len(df)):
            ax.plot([0, 1], [before.iloc[i], after.iloc[i]], 'o-', alpha=0.2, color='gray')

        ax.errorbar([0, 1], [before.mean(), after.mean()],
                    yerr=[before.std(), after.std()],
                    fmt='o-', linewidth=3, markersize=10, color='red', label='Mean +/- SD')

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Unrefined', 'Refined'])
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(metric.replace("_", " ").title() + "\nChange: " + f"{after.mean() - before.mean():.2f}")
        ax.grid(True, alpha=0.3)
        ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


# ============================================================================
# Script body — snakemake injected as global by Snakemake script: directive
# ============================================================================

print("Loading data...")
df_unrefined = pd.read_csv(snakemake.input.unrefined)
df_refined   = pd.read_csv(snakemake.input.refined)

merged = pd.merge(df_unrefined, df_refined, on='image_id', suffixes=('_unrefined', '_refined'))
print(f"Loaded {len(merged)} matched image pairs")

metrics_to_compare = snakemake.params.metrics

print("Calculating changes...")
rows = []
for _, row in merged.iterrows():
    entry = {'image_id': row['image_id']}
    for metric in metrics_to_compare:
        u = row[f'{metric}_unrefined']
        r = row[f'{metric}_refined']
        entry[f'{metric}_unrefined']  = u
        entry[f'{metric}_refined']    = r
        entry[f'{metric}_change']     = r - u
        entry[f'{metric}_change_pct'] = (r - u) / u * 100 if u != 0 else 0.0
    rows.append(entry)

comparison_df = pd.DataFrame(rows)

Path(snakemake.output.comparison).parent.mkdir(parents=True, exist_ok=True)

comparison_df.to_csv(snakemake.output.comparison, index=False)
print(f"Saved comparison to {snakemake.output.comparison}")

print("Performing statistical tests...")
stats_results = {}
for metric in metrics_to_compare:
    stats_results[metric] = perform_statistical_tests(
        merged[f'{metric}_unrefined'],
        merged[f'{metric}_refined'],
        metric
    )

with open(snakemake.output.stats, 'w') as f:
    json.dump(stats_results, f, indent=2)
print(f"Saved statistics to {snakemake.output.stats}")

print("Creating plots...")
create_comparison_plot(comparison_df, metrics_to_compare, snakemake.output.plot)
print(f"Saved plot to {snakemake.output.plot}")

print("\n" + "=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)
for metric in metrics_to_compare:
    r = stats_results[metric]
    print(f"\n{metric}:")
    print(f"  Before: {r['before_mean']:.2f} +/- {r['before_std']:.2f}")
    print(f"  After:  {r['after_mean']:.2f} +/- {r['after_std']:.2f}")
    print(f"  Change: {r['mean_change']:.2f} ({r['mean_change_pct']:.1f}%)")
    print(f"  t-test:   p={r['paired_ttest']['p_value']:.4f} {'*' if r['paired_ttest']['significant'] else 'ns'}")
    if 'p_value' in r['wilcoxon']:
        print(f"  Wilcoxon: p={r['wilcoxon']['p_value']:.4f} {'*' if r['wilcoxon']['significant'] else 'ns'}")
    print(f"  Cohen's d: {r['effect_size']['cohens_d']:.3f} ({r['effect_size']['interpretation']})")
print("=" * 70)
