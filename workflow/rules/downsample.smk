# rules/downsample.smk

# Set these to whatever makes sense in your project
KINDS = ["segs_converted", "roi_masks"]
WIDTHS = [576, 1024]

rule downsample:
    """
    Downsample all images in a folder for a given dataset/kind/width.
    """
    input:
        # Directory containing the original images
        in_dir = "data/{dataset}/{kind}/"
    output:
        # Directory to store the downsampled images
        out_dir = directory("data/{dataset}/downsampled/{width}px/{kind}/")
    params:
        kind = "{kind}",
        # width is wildcard string; cast to int in the script
        width = lambda wc: int(wc.width),
        ext = ".png",
    script:
        "scripts/downsample_smk.py"

