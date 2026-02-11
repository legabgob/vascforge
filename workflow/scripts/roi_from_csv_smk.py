import csv
import ast
import numpy as np
from pathlib import Path
from PIL import Image


def make_circular_mask(h, w, cx, cy, r):
    yy, xx = np.ogrid[:h, :w]
    return (xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2


def parse_row(row, header):
    name = row[0].strip()
    success = row[1].strip().lower() in ("true", "1", "yes")

    if "bounds" in header:
        meta = ast.literal_eval(row[2].strip())
        H, W = meta["hw"]
        cx, cy = meta["center"]
        R = float(meta["radius"])
        return name, success, int(H), int(W), float(cx), float(cy), float(R)

    idx = {col: header.index(col) for col in ("h", "w", "cy", "cx", "radius")}
    H = int(float(row[idx["h"]]))
    W = int(float(row[idx["w"]]))
    cy = float(row[idx["cy"]])
    cx = float(row[idx["cx"]])
    R = float(row[idx["radius"]])
    return name, success, H, W, cx, cy, R


# ----------------- Snakemake Entry ------------------

csv_file   = Path(snakemake.input.csv)
out_dir    = Path(snakemake.output.out_dir)
img_dir    = Path(snakemake.params.img_dir) if snakemake.params.get("img_dir") else None
img_ext    = snakemake.params.get("img_ext", ".png")
skip_false = snakemake.params.get("skip_if_flag_false", False)

out_dir.mkdir(parents=True, exist_ok=True)

with open(csv_file, "r") as f:
    reader = csv.reader(f)
    header = next(reader)

    for row in reader:
        name, success, H, W, cx, cy, R = parse_row(row, header)

        if skip_false and not success:
            continue

        mask = make_circular_mask(H, W, cx, cy, R).astype(np.uint8) * 255

        # Optional resizing
        if img_dir:
            img_path = img_dir / f"{name}{img_ext}"
            if img_path.exists():
                with Image.open(img_path) as im:
                    W2, H2 = im.size
                if (H2, W2) != (H, W):
                    mask = np.array(
                        Image.fromarray(mask, mode="L").resize((W2, H2), Image.NEAREST),
                        dtype=np.uint8,
                    )

        out_path = out_dir / f"{name}.png"
        Image.fromarray(mask, mode="L").save(out_path)
        print(f"[OK] {out_path}")
