# workflow/rules/refined_rgb_to_labels.smk
from snakemake.io import directory

LABELS_OUT_ROOT = config.get("labels_out_root", "results/refined_labels")

# Optional thresholds in config.yaml:
# refined_to_labels:
#   tb: 100
#   th: 200
#   delta: 15
#   tie_vote: 1
LBL_CFG = config.get("refined_to_labels", {})

rule refined_rgb_to_labels:
    """
    Convert refined RGB segmentations (0–255) into grayscale label maps {0,1,2,3}.
    Batch converts a whole refined directory for each (dataset,k,res).
    """
    input:
        in_dir = "results/refined/{dataset}/k{k}/downsampled/{res}px"
    output:
        out_dir = directory(f"{LABELS_OUT_ROOT}" + "/{dataset}/k{k}/downsampled/{res}px")
    benchmark:
        "benchmarks/refined_rgb_to_labels/{dataset}_k{k}_{res}px.tsv"
    params:
        name = "refined_rgb_to_labels",
        time = "01:00:00",
        mem = 8000,
        threads = 1,
        ext = ".png",
        tb = int(LBL_CFG.get("tb", 100)),
        th = int(LBL_CFG.get("th", 200)),
        delta = int(LBL_CFG.get("delta", 15)),
        tie_vote = int(LBL_CFG.get("tie_vote", 1)),  # 1=True, 0=False
    script:
        "../scripts/refined_rgb_to_labels_dir_smk.py"

