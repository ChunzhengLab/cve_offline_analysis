import pandas as pd
import matplotlib.pyplot as plt

# 给定的 fs/fb 比例
ratios = {
    "00-10": 1.6,
    "10-20": 2.25,
    "20-30": 3.4,
    "30-40": 4.4,
    "40-50": 5.5,
    "50-60": 6.0
}

# 计算 fs, fb, fs^2, fb^2, 和 2*fs*fb
data = []
for label, r in ratios.items():
    fs = r / (1 + r)
    fb = 1 / (1 + r)
    fs2 = fs**2
    fb2 = fb**2
    cross = 2 * fs * fb
    data.append([label, fs, fb, fs2, cross, fb2, r])

# 创建 DataFrame
df = pd.DataFrame(data, columns=["Centrality(%)", "fs", "fb", "fs^2", "2fsfb", "fb^2", "fs/fb"])
df.set_index("Centrality(%)", inplace=True)

# 创建三个子图
fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(18, 5))

# 图1：fs/fb 比例
df[["fs/fb"]].plot(kind="bar", ax=ax0, legend=False)
ax0.set_title("fs/fb Ratio by Centrality(%)")
ax0.set_ylabel("fs/fb")
ax0.set_xlabel("Centrality(%)")

# 图2：堆叠 fs 和 fb
df[["fs", "fb"]].plot(kind="bar", stacked=True, ax=ax1)
ax1.set_title("Stacked fs and fb by Centrality(%)")
ax1.set_ylabel("Value")
ax1.set_xlabel("Centrality(%)")
ax1.legend(title="Component")

# 图3：堆叠 fs²、2fsfb 和 fb²
df[["fs^2", "2fsfb", "fb^2"]].plot(kind="bar", stacked=True, ax=ax2)
ax2.set_title("Proportion of fs^2, 2fsfb, and fb^2")
ax2.set_ylabel("Proportion")
ax2.set_xlabel("Centrality(%)")
ax2.legend(title="Component")

plt.tight_layout()
plt.savefig("binomial_plot.png")
