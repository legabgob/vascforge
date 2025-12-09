def find_meta_like_csv(wc):
    base = f"./{wc.dataset}/seg_legacy"
    for name in ("meta.csv", "bounds.csv"):
        p = f"{base}/{name}"
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"No meta.csv or bounds.csv in {base}/")

rule copy_meta_csv:
    input:
        src = find_meta_like_csv
    output:
        dst = "data/{dataset}/meta/meta.csv"
    run:
        os.makedirs(os.path.dirname(output.dst), exist_ok=True)
        shell("cp {input.src} {output.dst}")

