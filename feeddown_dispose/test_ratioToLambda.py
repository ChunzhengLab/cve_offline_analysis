from pathlib import Path
import pandas as pd

BASE = {
    "1.0": Path("Proton_ratioToLambda1.0/finalise_feeddown_dispose_default.csv"),
    "0.9": Path("Proton_ratioToLambda0.9/finalise_feeddown_dispose_default.csv"),
    "0.8": Path("Proton_ratioToLambda0.8/finalise_feeddown_dispose_default.csv"),
}

keys = ["diff_type", "diff_bin", "pair_type"]

for r, f in BASE.items():
    df = pd.read_csv(f)
    df = df[df["centrality"] <= 60]          # 同样的过滤
    grp = df.groupby(keys).size()
    print(f"\nratio {r}: {len(df)} rows;  {len(grp)} unique (dT,dB,pT)")
    print(grp.head(10))                      # 先看前10行
