#!/usr/bin/env python3
"""
Aggregate all per-image connectivity summary CSVs into one master CSV.
Unrefined data is assigned k=0, refined data keeps its k value.
"""

import pandas as pd
import re
from pathlib import Path


connectivity_out = snakemake.params.connectivity_out
output_path      = snakemake.output.master_csv

all_dfs = []

for csv_file in snakemake.input.summaries:
    csv_path = Path(csv_file)
    try:
        df = pd.read_csv(csv_path)

        # Parse path relative to connectivity_out
        # Refined:   {connectivity_out}/refined/{dataset}[/{other_dir}]/k{k}/downsampled/{res}px/summary.csv
        # Unrefined: {connectivity_out}/unrefined/{dataset}[/{other_dir}]/downsampled/{res}px/summary.csv
        rel   = csv_path.relative_to(connectivity_out)
        parts = rel.parts

        refined = parts[0] == 'refined'

        if refined:
            k_idx     = next(i for i, p in enumerate(parts) if re.fullmatch(r'k\d+', p))
            k         = int(parts[k_idx][1:])
            res       = int(parts[k_idx + 2].replace('px', ''))
            dataset   = parts[1]
            other_dir = parts[2] if k_idx == 3 else None
        else:
            ds_idx    = next(i for i, p in enumerate(parts) if p == 'downsampled')
            k         = 0
            res       = int(parts[ds_idx + 1].replace('px', ''))
            dataset   = parts[1]
            other_dir = parts[2] if ds_idx == 3 else None

        df['dataset']    = dataset
        df['other_dir']  = other_dir
        df['k']          = k
        df['resolution'] = res

        all_dfs.append(df)

    except Exception as e:
        print(f"WARNING: skipping {csv_file}: {e}")

if not all_dfs:
    raise RuntimeError("No summary CSVs could be loaded")

master = pd.concat(all_dfs, ignore_index=True)

meta_cols   = ['dataset', 'other_dir', 'k', 'resolution', 'image_id']
metric_cols = [c for c in master.columns if c not in meta_cols]
master      = master[meta_cols + metric_cols].sort_values(['dataset', 'other_dir', 'image_id', 'k'])

Path(output_path).parent.mkdir(parents=True, exist_ok=True)
master.to_csv(output_path, index=False)

print(f"Master CSV: {len(master)} rows x {len(master.columns)} columns")
print(f"Datasets:   {sorted(master['dataset'].unique())}")
print(f"k values:   {sorted(master['k'].unique())}")
print(f"Saved to:   {output_path}")
