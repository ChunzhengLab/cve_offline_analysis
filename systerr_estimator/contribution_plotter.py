import pandas as pd
import matplotlib.pyplot as plt
import re
import os
import numpy as np

def parse_data_block(lines, block_id, pair_type):
    header = lines[1].strip().split(', ')
    delta_start = lines.index('--- Delta contributions ---\n') + 1
    gamma_start = lines.index('--- Gamma contributions ---\n')

    # 跳过 header 行
    delta_lines = [line for line in lines[delta_start:gamma_start] if not line.startswith("Centrality")]
    gamma_lines = [line for line in lines[gamma_start + 1:] if not line.startswith("Centrality")]

    def parse_contributions(lines, source_type):
        data = []
        for line in lines:
            parts = line.strip().split(', ')
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            
            # 排除中心度60-70%的点
            centrality_str = row["Centrality"]
            if '%' in centrality_str:
                # 处理类似 "60-70%" 的格式
                centrality_parts = centrality_str.replace('%', '').split('-')
                if len(centrality_parts) == 2:
                    cent_low = float(centrality_parts[0])
                    cent_high = float(centrality_parts[1])
                    if cent_low >= 60 and cent_high <= 70:
                        continue  # 跳过60-70%的点
                else:
                    cent_val = float(centrality_parts[0])
                    if cent_val >= 60 and cent_val <= 70:
                        continue  # 跳过60-70%的点
            else:
                # 处理纯数字格式
                cent_val = float(centrality_str)
                if cent_val >= 60 and cent_val <= 70:
                    continue  # 跳过60-70%的点
            
            row["SourceType"] = source_type
            row["PairType"] = pair_type
            row["BlockID"] = block_id
            data.append(row)
        return data

    data = parse_contributions(delta_lines, "Delta") + parse_contributions(gamma_lines, "Gamma")
    df = pd.DataFrame(data)
    for col in header[1:]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df["CentralityNum"] = df["Centrality"].apply(lambda c: float(c.replace('%', '').split('-')[0]) if '%' in c else float(c))
    return df, header[1:]


def parse_full_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    all_data = []
    i = 0
    current_block = None
    while i < len(lines):
        line = lines[i].strip()
        if re.match(r'^\[.*\], \[.*\]$', line):
            current_block = line.replace('[', '').replace(']', '').replace(', ', '_')  # e.g. DEta_0.5
        elif line.startswith('=== Pair Type:'):
            pair_type = line.split(':')[1].replace('=', '').strip()
            block_lines = lines[i:i+50]
            if '--- Delta contributions ---\n' in block_lines and '--- Gamma contributions ---\n' in block_lines:
                try:
                    df_block, _ = parse_data_block(block_lines, current_block, pair_type)
                    all_data.append(df_block)
                except Exception as e:
                    print(f"Error processing block at line {i}: {e}")
        i += 1

    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

def save_plots_as_pdf(df_all, output_dir="plots_pdf"):
    os.makedirs(output_dir, exist_ok=True)
    grouped = df_all.groupby(["BlockID", "PairType"])
    for (block_id, pair_type), df in grouped:
        df_delta = df[df["SourceType"] == "Delta"]
        df_gamma = df[df["SourceType"] == "Gamma"]
        if df_delta.empty or df_gamma.empty:
            continue

        variables = [col for col in df.columns if col not in ["Centrality", "CentralityNum", "SourceType", "PairType", "BlockID"]]

        # pivot为宽表，确保每个中心度唯一
        df_delta_pivot = df_delta.groupby("CentralityNum")[variables].mean()
        df_gamma_pivot = df_gamma.groupby("CentralityNum")[variables].mean()

        # 获取一个大的colormap
        cmap = plt.get_cmap('tab20')
        num_colors = len(variables)
        colors = [cmap(i % cmap.N) for i in range(num_colors)]

        fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

        # Delta line+marker
        for idx, var in enumerate(variables):
            axes[0].plot(df_delta_pivot.index, df_delta_pivot[var], label=var, marker='o', color=colors[idx])
        axes[0].set_xlabel("Centrality (%)")
        axes[0].set_title("Delta Contributions")
        axes[0].grid(True)
        axes[0].set_xticks(df_delta_pivot.index)
        axes[0].set_xticklabels(df_delta_pivot.index, rotation=45)
        axes[0].legend(fontsize='small', ncol=2)

        # Gamma line+marker
        for idx, var in enumerate(variables):
            axes[1].plot(df_gamma_pivot.index, df_gamma_pivot[var], label=var, marker='s', linestyle='--', color=colors[idx])
        axes[1].set_xlabel("Centrality (%)")
        axes[1].set_title("Gamma Contributions")
        axes[1].grid(True)
        axes[1].set_xticks(df_gamma_pivot.index)
        axes[1].set_xticklabels(df_gamma_pivot.index, rotation=45)
        axes[1].legend(fontsize='small', ncol=2)

        axes[0].set_ylabel("Contribution")
        plt.suptitle(f"{block_id} — Pair Type: {pair_type}", fontsize=14)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        pdf_path = os.path.join(output_dir, f"{block_id}_{pair_type}.pdf")
        plt.savefig(pdf_path)
        plt.close()
        print(f"Saved: {pdf_path}")

# === 主程序 ===
if __name__ == "__main__":
    input_path = "contrib_formatted_Hadron.txt"  # 替换为你的实际路径
    df_all = parse_full_file(input_path)
    save_plots_as_pdf(df_all, output_dir="plots_pdf_hadron")
