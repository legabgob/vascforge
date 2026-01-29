# workflow/rules/filter_downsampled_dirs.smk
"""
Filter downsampled directories to contain only square images based on meta_filtered_square.csv
"""
from snakemake.io import directory

rule filter_downsampled_dirs:
    input:
        meta_square = "data/{dataset}/meta/meta_filtered_square.csv",
        segs = "data/{dataset}/downsampled/{res}px/segs_converted",
        masks = "data/{dataset}/downsampled/{res}px/roi_masks_binarized",
    output:
        segs = directory("data/{dataset}/downsampled/{res}px/segs_converted_square"),
        masks = directory("data/{dataset}/downsampled/{res}px/roi_masks_binarized_square"),
    log:
        "logs/filter_downsampled_dirs/{dataset}_{res}px.log"
    run:
        import pandas as pd
        import shutil
        from pathlib import Path
        
        # Read CSV to get list of square images
        df = pd.read_csv(input.meta_square)
        id_col = 'id' if 'id' in df.columns else df.columns[0]
        square_ids = set(df[id_col].astype(str))
        
        # Create output directories
        segs_out = Path(output.segs)
        masks_out = Path(output.masks)
        segs_out.mkdir(parents=True, exist_ok=True)
        masks_out.mkdir(parents=True, exist_ok=True)
        
        segs_in = Path(input.segs)
        masks_in = Path(input.masks)
        
        copied = 0
        skipped = 0
        
        with open(log[0], 'w') as logfile:
            logfile.write(f"Filtering based on {input.meta_square}\n")
            logfile.write(f"Square IDs in CSV: {len(square_ids)}\n\n")
            
            for seg_file in sorted(segs_in.glob("*.png")):
                # Extract ID from filename (remove extension)
                img_id = seg_file.stem
                
                if img_id in square_ids:
                    # Copy segmentation
                    shutil.copy2(seg_file, segs_out / seg_file.name)
                    
                    # Copy corresponding mask
                    mask_file = masks_in / seg_file.name
                    if mask_file.exists():
                        shutil.copy2(mask_file, masks_out / seg_file.name)
                        copied += 1
                        logfile.write(f"KEPT: {seg_file.name}\n")
                    else:
                        logfile.write(f"WARNING: Mask missing for {seg_file.name}\n")
                else:
                    skipped += 1
                    logfile.write(f"SKIPPED: {seg_file.name} (not in square CSV)\n")
            
            logfile.write(f"\n=== SUMMARY ===\n")
            logfile.write(f"Copied: {copied}\n")
            logfile.write(f"Skipped: {skipped}\n")
        
        print(f"Filtered {wildcards.dataset} at {wildcards.res}px: kept {copied}, skipped {skipped}")
