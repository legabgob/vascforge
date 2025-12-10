# Automatically discover samples from input folder
from pathlib import Path

ROOT = Path("/HDD/data/relabelling-project")

# Map (dataset, sample) -> full path to grayscale AV label
input_map = {}

for p in ROOT.rglob("seg_legacy/av/*.png"):
    # Expect something like:
    # /HDD/.../relabelling-project/{dataset}/[optional other dirs]/seg_legacy/av/{sample}.png
    parts = p.parts
    try:
        idx = parts.index("relabelling-project")
    except ValueError:
        continue  # shouldn't happen, but be safe

    if idx + 2 >= len(parts):
        continue

    dataset = parts[idx + 1]          # e.g. "FIVES" or "Fundus-AVSeg"
    sample = p.stem                   # filename without .png

    # If the same (dataset, sample) appears in multiple subdirs, last one wins.
    input_map[(dataset, sample)] = str(p)

PAIRS = sorted(input_map.keys())
DATASETS = sorted({d for d, _ in PAIRS})
SAMPLES = sorted({s for _, s in PAIRS})


def find_gray_av_input(wc):
    key = (wc.dataset, wc.sample)
    try:
        return input_map[key]
    except KeyError:
        raise ValueError(
            f"No grayscale AV input found for dataset={wc.dataset}, sample={wc.sample}. "
            f"Known pairs: {len(input_map)}"
        )

rule gray_to_rgb:
    input:
        gray = find_gray_av_input
    output:
        rgb = "data/{dataset}/segs_converted/{sample}.png"
    script:
        "scripts/gray2rgb_smk.py"
