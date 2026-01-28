##!/usr/bin/env python3
"""
Standalone VascX feature extraction - extracted from run_full_pipeline.py
Usage: python vascx_extract_standalone.py <dataset_dir> <output_file> [--n-jobs N]
"""
import argparse
from pathlib import Path
from vascx.fundus.loader import RetinaLoader
from vascx.utils.analysis import extract_in_parallel


def main():
    parser = argparse.ArgumentParser(description="Extract VascX features from a dataset")
    parser.add_argument("ds_path", type=Path, help="Path to dataset directory (with original/, av/, meta.csv)")
    parser.add_argument("output_file", type=Path, help="Output TSV file path")
    parser.add_argument("--n-jobs", type=int, default=64, help="Number of parallel jobs")
    parser.add_argument("--av-subfolder", type=str, default="av", help="Name of AV subfolder")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.ds_path.exists():
        raise ValueError(f"Dataset path does not exist: {args.ds_path}")
    
    if not (args.ds_path / "original").exists():
        raise ValueError(f"Missing 'original/' in {args.ds_path}")
    
    if not (args.ds_path / args.av_subfolder).exists():
        raise ValueError(f"Missing '{args.av_subfolder}/' in {args.ds_path}")
    
    # Create output directory
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading dataset from {args.ds_path}")
    loader = RetinaLoader.from_folder(args.ds_path, av_subfolder=args.av_subfolder)
    
    print(f"Extracting features with {args.n_jobs} jobs...")
    res = extract_in_parallel(loader.to_dict(), "bergmann", n_jobs=args.n_jobs)
    
    print(f"Saving results to {args.output_file}")
    res.to_csv(args.output_file, sep='\t', na_rep='NaN', index=True)
    
    print(f"✓ Extracted features for {len(res)} images")
    print(f"✓ Feature columns: {len(res.columns)}")


if __name__ == "__main__":
    main()
