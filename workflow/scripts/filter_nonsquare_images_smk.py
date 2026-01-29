"""
Filter out non-square images from segmentation and mask directories.
Reads a CSV to identify which images to keep, then removes non-square images
from both directories.
"""
import sys
from pathlib import Path
import pandas as pd
from PIL import Image
import shutil

def main(meta_csv, segs_dir, masks_dir, segs_out, masks_out, log_file):
    """
    Filter directories to contain only square images.
    
    Args:
        meta_csv: Path to metadata CSV with image IDs
        segs_dir: Input segmentation directory
        masks_dir: Input mask directory
        segs_out: Output segmentation directory (filtered)
        masks_out: Output mask directory (filtered)
        log_file: Path to log file
    """
    meta_csv = Path(meta_csv)
    segs_dir = Path(segs_dir)
    masks_dir = Path(masks_dir)
    segs_out = Path(segs_out)
    masks_out = Path(masks_out)
    log_file = Path(log_file)
    
    # Create output directories
    segs_out.mkdir(parents=True, exist_ok=True)
    masks_out.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Read metadata
    df = pd.read_csv(meta_csv)
    id_col = 'id' if 'id' in df.columns else df.columns[0]
    
    # Get all image files from segmentation directory
    seg_files = sorted(segs_dir.glob("*.png"))
    
    removed = []
    copied = []
    
    with open(log_file, 'w') as log:
        log.write(f"Processing {len(seg_files)} images from {segs_dir}\n")
        log.write(f"Output: segs={segs_out}, masks={masks_out}\n\n")
        
        for seg_file in seg_files:
            # Check if image is square
            with Image.open(seg_file) as img:
                width, height = img.size
                
            if width != height:
                removed.append((seg_file.name, width, height))
                log.write(f"REMOVED (non-square): {seg_file.name} ({width}x{height})\n")
            else:
                # Copy both segmentation and corresponding mask
                mask_file = masks_dir / seg_file.name
                
                # Copy segmentation
                shutil.copy2(seg_file, segs_out / seg_file.name)
                
                # Copy mask if it exists
                if mask_file.exists():
                    shutil.copy2(mask_file, masks_out / seg_file.name)
                    copied.append(seg_file.name)
                    log.write(f"KEPT (square): {seg_file.name} ({width}x{height})\n")
                else:
                    log.write(f"WARNING: Mask not found for {seg_file.name}\n")
        
        log.write(f"\n=== SUMMARY ===\n")
        log.write(f"Total images: {len(seg_files)}\n")
        log.write(f"Square images kept: {len(copied)}\n")
        log.write(f"Non-square images removed: {len(removed)}\n")
        
        if removed:
            log.write(f"\nRemoved images:\n")
            for name, w, h in removed:
                log.write(f"  {name}: {w}x{h}\n")
    
    print(f"Filtered {len(seg_files)} images:")
    print(f"  Kept (square): {len(copied)}")
    print(f"  Removed (non-square): {len(removed)}")
    
    if removed:
        print(f"\nRemoved non-square images:")
        for name, w, h in removed[:5]:  # Show first 5
            print(f"  {name}: {w}x{h}")
        if len(removed) > 5:
            print(f"  ... and {len(removed) - 5} more")

if __name__ == "__main__":
    # When run as a Snakemake script, use the snakemake object
    # When run standalone, use command-line args
    if 'snakemake' in globals():
        main(
            snakemake.input.meta,
            snakemake.input.segs,
            snakemake.input.masks,
            snakemake.output.segs,
            snakemake.output.masks,
            snakemake.log[0]
        )
    else:
        if len(sys.argv) != 7:
            print("Usage: filter_nonsquare_images_smk.py <meta_csv> <segs_dir> <masks_dir> <segs_out> <masks_out> <log_file>")
            sys.exit(1)
        main(*sys.argv[1:])
