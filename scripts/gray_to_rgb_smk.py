import os
import argparse
import numpy as np
from PIL import Image


def gray_labels_to_rgb(arr: np.ndarray) -> np.ndarray:
    # Convert grayscale labels to RGB colors

    if arr.ndim != 2:
        raise ValueError("Input array must be 2-dimensional")

    arr = arr.astype(np.int32, copy=False) # Make sure we work in integers

    H, W = arr.shape
    rgb = np.zeros((H, W, 3), dtype=np.uint8)

    # Define each class
    arteries = (arr == 1)
    veins = (arr == 2)
    crossings = (arr == 3)

    rgb[arteries, 0] = 255  # Red channel for arteries
    rgb[arteries, 2] = 255
    rgb[veins, 1] = 255     # Green channel for veins
    rgb[veins, 2] = 255
    rgb[crossings] = [255, 255, 255] # White for crossings

    return rgb

def convert(input_path:str, output_path:str):
    # Load the grayscale image
    img = Image.open(input_path).convert('L')
    arr = np.array(img)
    rgb = gray_labels_to_rgb(arr)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    Image.fromarray(rgb, mode="RGB").save(output_path)
    print(f"Saved {output_path}")

# ---------- Snakemake entry point ----------

# This file is executed by Snakemake with a 'snakemake' object in scope.
input_path = str(snakemake.input[0])
output_path = str(snakemake.output[0])

convert(input_path, output_path)

