#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FLOAT_ATOL = 1e-9


def load_data(file, tag, cent_max=60):
    df = pd.read_csv(file)
    df = df[df['centrality'] <= cent_max].copy()
    df['tag'] = tag
    return df


def select_common(df1, df2):
    keys = ['diff_type', 'diff_bin', 'pair_type']
    def combos(df):
        return set(map(tuple, df[keys].round(8).to_numpy()))
    return set.intersection(combos(df1), combos(df2))


def sel(df, combo):
    dtyp, dbin, ptyp = combo
    return df[
        (df['diff_type'] == dtyp) &
        (np.isclose(df['diff_bin'], dbin, atol=FLOAT_ATOL)) &
        (df['pair_type'] == ptyp)
    ]


def plot_combo_grouped(df1, df2, dtyp, dbin, val, err, output_dir):
    for pair_group in [['OS', 'SS'], ['Del']]:
        fig, (axL, axR) = plt.subplots(
            ncols=2, figsize=(12, 4), constrained_layout=True,
            gridspec_kw={"width_ratios": [3, 2]}
        )

        for ptyp in pair_group:
            df1c = sel(df1, (dtyp, dbin, ptyp)).sort_values('centrality')
            df2c = sel(df2, (dtyp, dbin, ptyp)).sort_values('centrality')

            if df1c.empty or df2c.empty:
                continue

            style = '-' if ptyp == 'OS' else '--' if ptyp == 'SS' else ':'
            label_suffix = f" ({ptyp})"

            # 左图
            for df, label, shift in [(df1c, 'w NUA', -1), (df2c, 'wo NUA', 1)]:
                x = df['centrality'] + shift
                axL.errorbar(x, df[val], yerr=df[err], fmt='o', linestyle=style, capsize=3,
                             label=f'{label}{label_suffix}')

            # 右图
            centr = df1c['centrality'].values
            ratio = np.abs(df2c[val].values - df1c[val].values) / np.sqrt(
                df2c[err].values**2 + df1c[err].values**2
            )
            axR.plot(centr, ratio, 's-', linestyle=style, label=f'{ptyp} Barlow')

        axL.set_xlabel('Centrality (%)')
        axL.set_ylabel(val)
        axL.set_title(f'{val} - {dtyp}, {dbin}, {pair_group}')
        axL.legend()
        axL.grid(True, ls='--', alpha=0.3)

        axR.axhline(1, ls='--', color='gray', alpha=0.6)
        axR.set_xlabel('Centrality (%)')
        axR.set_ylabel('Barlow ratio')
        axR.set_title(r"$|x_2 - x_1| / \sqrt{\sigma_1^2 + \sigma_2^2}$")
        axR.grid(True, ls='--', alpha=0.3)
        axR.legend()

        fname = f'{val}_GroupedBarlow_dT-{dtyp}_dB-{dbin}_pG-{"_".join(pair_group)}.pdf'
        plt.savefig(output_dir / fname)
        plt.close(fig)
        print(f'  ✓ Group plot: {val}, {dtyp}, {dbin}, {pair_group}')


def main():
    parser = argparse.ArgumentParser(description='Compare two systematic error files using Barlow method')
    parser.add_argument('--ref', required=True, help='Reference CSV (e.g. 0.9)')
    parser.add_argument('--test', required=True, help='Test CSV (e.g. 1.0)')
    parser.add_argument('--output', default='plots_compare', help='Output folder')
    parser.add_argument('--draw-gamma', action='store_true', help='Also draw gamma plots')
    parser.add_argument('--max-cent', type=int, default=60)
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    df1 = load_data(args.ref, 'ref', args.max_cent)
    df2 = load_data(args.test, 'test', args.max_cent)

    common = select_common(df1, df2)
    print(f'[INFO] Found {len(common)} common combinations')

    val_errs = [('delta', 'delta_err')]
    if args.draw_gamma:
        val_errs.append(('gamma', 'gamma_err'))

    seen = set()
    total_plots = 0

    for combo in sorted(common):
        dtyp, dbin, ptyp = combo
        key = (dtyp, dbin)
        if key in seen:
            continue
        seen.add(key)
        for val, err in val_errs:
            plot_combo_grouped(df1, df2, dtyp, dbin, val, err, output_dir)
            total_plots += 1

    print(f'\n✓ Done. {total_plots} grouped plots saved in {output_dir.resolve()}')


if __name__ == '__main__':
    main()
