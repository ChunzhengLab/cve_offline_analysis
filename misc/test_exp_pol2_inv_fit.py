import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import argparse

# 拟合函数：指数 + n阶多项式
def exp_poly_func(x, *params):
    A, k = params[0], params[1]
    poly_params = params[2:]
    poly = sum(c * (x ** i) for i, c in enumerate(poly_params))
    return A * np.exp(k * x) + poly

# 简单指数，用来做预拟合
def simple_exp(x, A, k):
    return A * np.exp(k * x)

# 解析参数
parser = argparse.ArgumentParser(description="Fit invariant mass distributions with exp + poly")
parser.add_argument('--input', type=str, required=True, help='Input CSV file')
parser.add_argument('--signal-min', type=float, default=1.103, help='Signal region min')
parser.add_argument('--signal-max', type=float, default=1.128, help='Signal region max')
parser.add_argument('--poly-order', type=int, default=2, help='Polynomial order')
parser.add_argument('--normalize', action='store_true', help='Normalize counts to max=1')
parser.add_argument('--output', type=str, required=True, help='Output PDF filename')
parser.add_argument('--maxfev', type=int, default=1000000, help='Maximum function evaluations for curve_fit')
args = parser.parse_args()

# 读取数据
df = pd.read_csv(args.input)

# 获取所有 diff_type, diff_bin, pair_type
diff_types = df['diff_type'].unique()

# 创建 PDF
with PdfPages(args.output) as pdf:
    for diff_type in diff_types:
        diff_bin_values = df[df['diff_type'] == diff_type]['diff_bin'].unique()
        for diff_bin in diff_bin_values:
            pair_types = df[
                (df['diff_type'] == diff_type) &
                (df['diff_bin'] == diff_bin)
            ]['pair_type'].unique()

            for pair_type in pair_types:
                cent_values = sorted(
                    df[
                        (df['diff_type'] == diff_type) &
                        (df['diff_bin'] == diff_bin) &
                        (df['pair_type'] == pair_type)
                    ]['centrality'].unique()
                )

                ncols, nrows = 3, 3
                fig, axes = plt.subplots(nrows, ncols, figsize=(15, 12))
                axes = axes.flatten()

                for idx, cent in enumerate(cent_values):
                    if idx >= len(axes):
                        break

                    subdf = df[
                        (df['diff_type'] == diff_type) &
                        (df['diff_bin'] == diff_bin) &
                        (df['pair_type'] == pair_type) &
                        (df['centrality'] == cent)
                    ]

                    x = subdf['inv_mass_center'].values
                    y = subdf['inv_mass_counts'].values

                    if args.normalize:
                        y = y / np.max(y)

                    # Sideband for background fit
                    side_mask = (x < args.signal_min) | (x > args.signal_max)
                    x_fit = x[side_mask]
                    y_fit = y[side_mask]

                    ax = axes[idx]
                    ax.plot(x, y, 'ko', markersize=4, label='Data')

                    # 绘制信号区域
                    ax.axvspan(args.signal_min, args.signal_max, color='yellow', alpha=0.3, label='Signal Region')

                    # Fitting
                    try:
                        p0_exp, _ = curve_fit(simple_exp, x_fit, y_fit, p0=[np.max(y_fit), -5], maxfev=args.maxfev)
                        p0 = list(p0_exp) + [0.0] * (args.poly_order + 1)
                        popt, _ = curve_fit(exp_poly_func, x_fit, y_fit, p0=p0, maxfev=args.maxfev)

                        # 绘制拟合曲线
                        x_plot = np.linspace(min(x), max(x), 500)
                        y_plot = exp_poly_func(x_plot, *popt)
                        ax.plot(x_plot, y_plot, 'r-', label='Fit')

                        # 显示拟合参数
                        param_text = f"exp({popt[1]:.2g}x)"
                        for i, c in enumerate(popt[2:]):
                            if abs(c) > 1e-10:
                                param_text += f" + {c:.2g}x^{i}"
                        ax.text(0.5, 0.95, param_text, ha='center', va='top', transform=ax.transAxes, fontsize=8)

                    except Exception as e:
                        ax.text(0.5, 0.5, 'Fit Failed', ha='center', va='center', color='red', fontsize=12)

                    ax.set_title(f'Centrality {cent}%')
                    ax.grid(True)

                # Hide extra subplots
                for j in range(idx+1, len(axes)):
                    axes[j].axis('off')

                fig.suptitle(f'diff_type={diff_type}, diff_bin={diff_bin}, pair_type={pair_type}', fontsize=16)
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                pdf.savefig(fig)
                plt.close(fig)

print(f"全部完成！结果保存在 {args.output}")