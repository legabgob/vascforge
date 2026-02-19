# workflow/rules/downsample.smk
from snakemake.io import directory

KINDS = ["segs_converted", "roi_masks_binarized"]
WIDTHS = [str(r) for r in config.get("resolutions", ["576", "1024"])]

rule downsample:
    """
    Downsample all images in a folder for a given dataset/kind/width.
    """
    input:
        in_dir = "data/{dataset}/{kind}"
    output:
        out_dir = directory("data/{dataset}/downsampled/{width}px/{kind}")
    params:
        kind = "{kind}",
        width = lambda wc: int(wc.width),
        ext = ".png",
    script:
        "../scripts/downsample_dir_smk.py"

