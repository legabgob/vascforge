# Automatically discover samples from input folder
from snakemake.io import glob_wildcards

# adjust pattern to your actual files
DATASETS, SAMPLES = glob_wildcards(
    "/HDD/data/relabelling-project/{dataset}/seg_legacy/av/{sample}.png"
)

rule gray_to_rgb:
    """
    Convert grayscale label images (0,1,2,3) to RGB A/V/BV mapping.
    """
    input:
        "/HDD/data/relabelling-project/{dataset}/seg_legacy/av/{sample}.png"
    output:
        "/SSD/home/gabriel/segmentation-relabelling/data/{dataset}/segs_converted/{sample}.png"
    script:
        "scripts/gray2rgb_smk.py"
