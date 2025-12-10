# Automatically discover samples from input folder
from snakemake.io import glob_wildcards
import os

pattern1 = "/HDD/data/relabelling-project/{dataset}/seg_legacy/av/{sample}.png"
pattern2 = "/HDD/data/relabelling-project/{dataset}/{other_dir}/seg_legacy/av/{sample}.png"

datasets1, samples1 = glob_wildcards(pattern1)
datasets2, other_dirs2, samples2 = glob_wildcards(pattern2)

input_map = {}

for d, s in zip(datasets1, samples1):
    input_map[(d, s)] = f"/HDD/data/relabelling-project/{d}/seg_legacy/av/{s}.png"

for d, odir, s in zip(datasets2, other_dirs2, samples2):
    input_map[(d, s)] = f"/HDD/data/relabelling-project/{d}/{odir}/seg_legacy/av/{s}.png"

PAIRS = sorted(input_map.keys())

# Optional: expose these for rule all in Snakefile
DATASETS = sorted({d for d, _ in PAIRS})
SAMPLES  = sorted({s for _, s in PAIRS})


def find_gray_av_input(wc):
    return input_map[(wc.dataset, wc.sample)]


rule gray_to_rgb:
    input:
        gray = find_gray_av_input
    output:
        rgb = "/data/{dataset}/segs_converted/{sample}.png"
    script:
        "scripts/gray2rgb_smk.py"
