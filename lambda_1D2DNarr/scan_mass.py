import os
import subprocess
import shutil

# === 配置参数 ===
base_script = "fit_obvs.py"  # 主脚本文件名
task = "default"             # ✅ 必须是default，不能修改
dataset = "LHC18q"           # 按实际情况填写（或使用自动推断也可以）

# === 精确计算bin宽度和range_mass列表 ===
bin_width = (1.115683 + 0.02 - (1.115683 - 0.02)) / 30
epsilon = 1e-6
bin_counts = [15, 10, 8, 6, 4, 3, 2, 1]
range_list = [round(n * bin_width + epsilon, 6) for n in bin_counts]

# === 遍历每个 range_mass，运行 narrMass 和 1Dfit ===
for r in range_list:
    r_str = f"{r:.6f}".replace(".", "p")
    output_dir = f"./outputs_scan/r{r_str}"
    os.makedirs(output_dir, exist_ok=True)

    for mode in ["narrMass", "1Dfit","2Dfit"]:
        print(f"\n>>> Running: range_mass={r:.6f}, mode={mode}")

        cmd = [
            "python", base_script,
            "-o", output_dir,
            "-t", task,  # 保持为default
            "-d", dataset,
            "--fit-mode", mode,
            "--range_mass", str(r)
        ]
        result = subprocess.run(cmd)

        # 检查命令是否成功执行
        if result.returncode != 0:
            print(f"警告: 命令执行失败，跳过重命名")
            continue

        # 运行后立即重命名文件
        particle_type = "Lambda"
        result_dir = os.path.join(output_dir, dataset, particle_type)

        # 原始文件名（fit_obvs.py生成的）
        original_csv = os.path.join(result_dir, f"fit_obvs_{task}.csv")
        original_root = os.path.join(result_dir, f"fit_obvs_plots_{task}.root")

        # 新文件名（平铺到outputs_scan目录）
        new_csv = os.path.join("./outputs_scan", f"fit_obvs_{task}_{mode}_range{r_str}.csv")
        new_root = os.path.join("./outputs_scan", f"fit_obvs_plots_{task}_{mode}_range{r_str}.root")

        # 重命名CSV文件
        if os.path.exists(original_csv):
            shutil.move(original_csv, new_csv)
            print(f"✅ 重命名CSV: {os.path.basename(new_csv)}")
        else:
            print(f"⚠️  未找到CSV文件: {original_csv}")

        # 重命名ROOT文件
        if os.path.exists(original_root):
            shutil.move(original_root, new_root)
            print(f"✅ 重命名ROOT: {os.path.basename(new_root)}")
        else:
            print(f"⚠️  未找到ROOT文件: {original_root}")

print("\n🎉 所有扫描完成！")
print("📁 输出文件格式示例：")
print("   fit_obvs_default_narrMass_range0p020000.csv")
print("   fit_obvs_default_1Dfit_range0p020000.csv")
print("   fit_obvs_plots_default_narrMass_range0p020000.root")
print("   fit_obvs_plots_default_1Dfit_range0p020000.root")
