#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot delta (或 gamma) 并在右侧子图绘制 Balow ratio
Balow 公式: |Δ_test − Δ_0.9| / sqrt(σ_test² + σ_0.9²)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============ 参数 ============
BASE_DIRS = {
    "1.0": Path("Proton_1.0"),
    "0.9": Path("Proton_0.9"),
    "0.8": Path("Proton_0.8"),
}
CSV_NAME   = "finalise_feeddown_dispose_default.csv"
OUT_DIR    = Path("plots")
CENT_MAX   = 60
DRAW_GAMMA = False
SHIFT_MAP  = {"0.8": -2, "0.9": 0, "1.0": 2}
FLOAT_ATOL = 1e-9
# ==============================

OUT_DIR.mkdir(exist_ok=True, parents=True)

# ---------- 读取 ----------
frames = {}
for r, folder in BASE_DIRS.items():
    fp = folder / CSV_NAME
    if not fp.exists():
        raise FileNotFoundError(fp)
    df = pd.read_csv(fp)
    df["ratio"] = r
    df = df[df["centrality"] <= CENT_MAX]
    frames[r] = df
    print(f"[INFO] ratio {r}: {len(df)} rows")

# 合并
data = pd.concat(frames.values(), ignore_index=True)

# 三份数据共同组合
keys = ["diff_type", "diff_bin", "pair_type"]
def combos(df): return set(map(tuple, df[keys].round(8).to_numpy()))
common = set.intersection(*(combos(df) for df in frames.values()))
print(f"[INFO] 共 {len(common)} 个组合\n")

def sel(df, c):
    dtyp, dbin, ptyp = c
    return df[
        (df["diff_type"] == dtyp)&
        (np.isclose(df["diff_bin"], dbin, atol=FLOAT_ATOL))&
        (df["pair_type"] == ptyp)
    ]

def plot_combo(sub, combo, val, err):
    # ------------ Figure ------------
    fig, (axL, axR) = plt.subplots(
        ncols=2, figsize=(12,4), constrained_layout=True,
        gridspec_kw={"width_ratios":[3,2]}
    )

    # 左：原曲线
    for ratio, g in sub.groupby("ratio"):
        g = g.sort_values("centrality")
        x = g["centrality"] + SHIFT_MAP[ratio]
        axL.errorbar(x, g[val], yerr=g[err],
                     marker="o", ls="-", capsize=3,
                     label=f"ratio {ratio}")
    axL.set_xlabel("Centrality (shifted)")
    axL.set_ylabel(val)
    dtyp, dbin, ptyp = combo
    axL.set_title(f"{val}: {dtyp}, {dbin}, {ptyp}")
    axL.grid(True, ls="--", alpha=.3)
    axL.legend()

    # 右：Balow
    def_d  = sub[sub["ratio"]=="0.9"].sort_values("centrality")
    centr  = def_d["centrality"].values

    for ratio in ("0.8","1.0"):
        test = sub[sub["ratio"]==ratio].sort_values("centrality")
        if test.empty: continue
        balow = np.abs(test[val].values - def_d[val].values) / (test[err].values - def_d[err].values)**2
        axR.plot(centr, balow, "s-", label=f"{ratio} vs 0.9")
    axR.set_xlabel("Centrality")
    axR.set_ylabel("Balow ratio")
    axR.set_title(r"$|\Delta-\Delta_{0.9}|/\sqrt{\sigma^2 - \sigma_{0.9}^2}$")
    axR.grid(True, ls="--", alpha=.3)
    axR.legend()

    # 保存
    fname = f"{val}_balow_dT-{dtyp}_dB-{dbin}_pT-{ptyp}.pdf"
    fig.suptitle("")   # 防止标题重叠
    fig.savefig(OUT_DIR/fname)
    plt.close(fig)
    print(f"[OK] saved {fname}")

plotted = 0
for combo in sorted(common):
    subset = sel(data, combo)
    if subset.empty: continue
    plot_combo(subset, combo, "delta", "delta_err"); plotted += 1
    if DRAW_GAMMA:
        plot_combo(subset, combo, "gamma", "gamma_err"); plotted += 1

print(f"\n✓ 完成，共生成 {plotted} 张 PDF，保存在 {OUT_DIR.resolve()}")
