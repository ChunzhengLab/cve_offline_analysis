#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Barlow systematics  ▶  smooth (per-category exp-fit & poly2) ▶  plot
可选 --original 把 default 行的 (delta,delta_err,gamma,gamma_err) 并入结果。
仅保留 centrality ≤ 60 (%) 的数据。
"""

import argparse, re, sys
from pathlib import Path
from typing import List
import numpy as np, pandas as pd

import matplotlib
matplotlib.use("Agg")            # 无显示环境
import matplotlib.pyplot as plt  # noqa: E402


# ═════════════════════ 辅助 ═════════════════════
def centre_of(label) -> float:
    if isinstance(label, (int, float)):
        return float(label)
    m = re.match(r"(\d*\.?\d+)\D+(\d*\.?\d+)", str(label))
    if m:
        a, b = map(float, m.groups())
        return 0.5 * (a + b)
    try:
        return float(label)
    except Exception:
        return np.nan


def exp_fit_centered(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    if len(y) < 3:
        return y.copy()
    shift = 0.0
    if np.any(y <= 0):
        shift = -y.min() + 1e-6
        y = y + shift
    ln_y = np.log(y)
    B, lnA = np.polyfit(x, ln_y, 1)
    return np.exp(lnA) * np.exp(B * x) - shift


def safe_tag(v) -> str:
    return "NA" if (v is None or (isinstance(v, str) and v.strip() == "")) else re.sub(r"[^\w\-.]", "_", str(v))


# ═════════════════════ 系统误差 + 平滑 ═════════════════════
def calc_syst(sys_df: pd.DataFrame, default_df: pd.DataFrame | None) -> pd.DataFrame:
    group_cols: List[str] = ["centrality", "pair_type"]
    if "diff_type" in sys_df and not sys_df["diff_type"].fillna("").eq("").all():
        group_cols.append("diff_type")
    if "diff_bin" in sys_df and not sys_df["diff_bin"].fillna("").eq("").all():
        group_cols.append("diff_bin")

    sys_df[group_cols] = sys_df[group_cols].fillna("")
    if default_df is not None:
        default_df[group_cols] = default_df[group_cols].fillna("")

    # ── (1) 计算 Δ/Γ 系统误差 ───────────────────────────────
    recs = []
    for key, g in sys_df.groupby(group_cols):
        rec = {c: v for c, v in zip(group_cols, key if isinstance(key, tuple) else (key,))}
        rec["delta_syst_err"] = np.sqrt((g.loc[g["is_delta_syssrc"], "delta_diff"] ** 2).sum())
        rec["gamma_syst_err"] = np.sqrt((g.loc[g["is_gamma_syssrc"], "gamma_diff"] ** 2).sum())
        recs.append(rec)
    res = pd.DataFrame(recs)

    # ── (2) 如果有 original，把默认值并进来 ────────────────
    if default_df is not None:
        default_keep = default_df[default_df["source"] == "default"][group_cols +
            ["delta", "delta_err", "gamma", "gamma_err"]]
        res = res.merge(default_keep, on=group_cols, how="left")

    # ── (3) 组合列列表，用于分块独立平滑 ───────────────────
    combo_cols = ["pair_type"]
    if "diff_type" in res.columns:
        combo_cols.append("diff_type")
    if "diff_bin" in res.columns:
        combo_cols.append("diff_bin")

    # 初始化平滑列
    for col in ["delta_syst_err_expfit", "gamma_syst_err_expfit",
                "delta_syst_err_pol2", "gamma_syst_err_pol2"]:
        res[col] = np.nan

    # ── (4) 在每个组合内单独做 exp-fit & poly2 ─────────────
    for _, sub_idx in res.groupby(combo_cols, sort=False).groups.items():
        idx = list(sub_idx)
        sub = res.loc[idx].copy()
        sub["c_mid"] = sub["centrality"].apply(centre_of)
        sub = sub.sort_values("c_mid")
        x = sub["c_mid"].to_numpy(float)

        for which in ["delta", "gamma"]:
            y = sub[f"{which}_syst_err"].to_numpy(float)
            res.loc[sub.index, f"{which}_syst_err_expfit"] = exp_fit_centered(x, y)
            res.loc[sub.index, f"{which}_syst_err_pol2"] = (
                np.polyval(np.polyfit(x, y, 2), x) if len(sub) >= 3 else y)

    return res


# ═════════════════════ 绘图 ═════════════════════
def plot_one(sub: pd.DataFrame, particle: str,
             pair: str, difftype: str, diffbin: str,
             out_dir: Path, which: str):
    sub = sub.copy()
    sub["c_mid"] = sub["centrality"].apply(centre_of)
    sub = sub.dropna(subset=["c_mid"]).sort_values("c_mid")
    if sub.empty:
        return
    x = sub["c_mid"].to_numpy()
    y_raw = sub[f"{which}_syst_err"].to_numpy()
    y_exp = sub[f"{which}_syst_err_expfit"].to_numpy()
    y_pol2 = sub[f"{which}_syst_err_pol2"].to_numpy()

    plt.figure(figsize=(6, 4))
    plt.plot(x, y_raw, "o", label="raw")
    plt.plot(x, y_exp, "-", label="exp-fit")
    plt.plot(x, y_pol2, "--", label="pol2")
    plt.xlabel("Centrality (%)")
    plt.ylabel(f"$\delta$ syst err" if which == "delta" else "$\gamma$ syst err")

    plt.title(f"{particle}: {pair}/{difftype}/{diffbin}")
    plt.legend(fontsize=8)
    plt.tight_layout()
    fname = f"{which}_err_vs_centrality_{particle}_{safe_tag(pair)}_{safe_tag(difftype)}_{safe_tag(diffbin)}.pdf"
    plt.savefig(out_dir / fname)
    plt.close()
    print(f"[PLOT] {fname}")


def make_plots(res: pd.DataFrame, particle: str, out_dir: Path):
    combo_cols = ["pair_type"]
    if "diff_type" in res.columns: combo_cols.append("diff_type")
    if "diff_bin"  in res.columns: combo_cols.append("diff_bin")

    for combo_vals, sub in res.groupby(combo_cols, dropna=False, sort=False):
        if not isinstance(combo_vals, tuple): combo_vals = (combo_vals,)
        while len(combo_vals) < 3: combo_vals += ("",)
        plot_one(sub, particle, *combo_vals, out_dir=out_dir, which="delta")
        plot_one(sub, particle, *combo_vals, out_dir=out_dir, which="gamma")


# ═════════════════════ 主程序 ═════════════════════
def main():
    ap = argparse.ArgumentParser(description="Barlow systematics (per-category fitting, optional default merge).")
    ap.add_argument("--input", required=True, help="syssrc_barlow_finalise_sys_<particle>.csv")
    ap.add_argument("--original", help="finalise_sys_<particle>.csv (含 default 行，可选)")
    ap.add_argument("--output", help="输出汇总 CSV")
    ap.add_argument("--verbose", action="store_true", help="打印统计摘要")
    args = ap.parse_args()

    inp = Path(args.input).resolve()
    if not inp.exists():  print(f"[ERROR] 输入 {inp} 不存在"); sys.exit(1)
    if args.original and not Path(args.original).exists():
        print(f"[ERROR] original {args.original} 不存在"); sys.exit(1)

    sys_df = pd.read_csv(inp)
    need = ["centrality", "pair_type", "delta_diff", "gamma_diff",
            "is_delta_syssrc", "is_gamma_syssrc"]
    if any(c not in sys_df.columns for c in need):
        print("[ERROR] 输入缺必要列"); sys.exit(1)

    # ≤ 60 %
    sys_df["c_mid"] = sys_df["centrality"].apply(centre_of)
    sys_df = sys_df[sys_df["c_mid"] <= 60].drop(columns="c_mid")
    if sys_df.empty: print("[ERROR] centrality ≤ 60 的数据为空"); sys.exit(1)

    default_df = None
    if args.original:
        df_orig = pd.read_csv(args.original)
        if "source" in df_orig.columns and "default" in df_orig["source"].values:
            default_df = df_orig[df_orig["source"] == "default"]
        else:
            print("[WARN] original 文件没有 'default' 行，忽略。")

    particle = inp.stem.replace("syssrc_barlow_finalise_sys_", "")
    res = calc_syst(sys_df, default_df)

    out_csv = Path(args.output) if args.output else inp.parent / f"finalise_syst_{particle}.csv"
    res.to_csv(out_csv, index=False)
    print(f"[CSV] {out_csv}")

    if args.verbose: print(res.describe())

    # Create plots_sys_smoothing directory
    plots_dir = out_csv.parent / "plots_sys_smoothing"
    plots_dir.mkdir(exist_ok=True)
    
    make_plots(res, particle, plots_dir)


if __name__ == "__main__":
    main()
