# Feeddown扣除前后观测量比较绘图

本目录包含用于比较feeddown扣除前后观测量变化的绘图脚本。

## 文件说明

- `plot_feeddown_comparison.py`: 核心绘图脚本，提供feeddown比较功能
- `run_feeddown_comparison.py`: 便捷运行脚本，使用默认配置
- `README_feeddown_plots.md`: 本说明文件

## 功能特点

1. **对比绘图**: 在同一图中显示feeddown扣除前后的观测量
2. **多种观测量**: 支持delta和gamma两种观测量类型
3. **全面比较**: 显示SS、OS、Del三种配对类型的比较
4. **相对变化**: 计算并显示feeddown扣除的相对影响
5. **多差分类型**: 支持Intg、SPt、DEta等不同差分类型

## 使用方法

### 方法1：使用便捷脚本（推荐）

```bash
# 直接运行，使用默认配置
python run_feeddown_comparison.py
```

这将自动：
- 使用默认的数据路径
- 生成所有差分类型的比较图
- 输出到 `plots/feeddown_comparison/` 目录

### 方法2：使用核心脚本

```bash
# 基本用法
python plot_feeddown_comparison.py \
    --original ../dataset_merger/Merged/Proton/finalise_default.csv \
    --feeddown ./Proton/finalise_feeddown_dispose_default.csv \
    --value_type delta \
    --output_dir ./plots/feeddown_comparison

# 生成所有差分类型的比较图
python plot_feeddown_comparison.py \
    --original ../dataset_merger/Merged/Proton/finalise_default.csv \
    --feeddown ./Proton/finalise_feeddown_dispose_default.csv \
    --value_type delta \
    --all_types \
    --output_dir ./plots/feeddown_comparison

# 生成特定差分类型的比较图
python plot_feeddown_comparison.py \
    --original ../dataset_merger/Merged/Proton/finalise_default.csv \
    --feeddown ./Proton/finalise_feeddown_dispose_default.csv \
    --value_type delta \
    --diff_type Intg \
    --diff_bin 0.5 \
    --output_dir ./plots/feeddown_comparison
```

## 参数说明

- `--original`: 原始数据文件路径
- `--feeddown`: Feeddown扣除后数据文件路径
- `--value_type`: 观测量类型 (`delta` 或 `gamma`)
- `--diff_type`: 差分类型 (如 `DEta`, `SPt`, `Intg`)
- `--diff_bin`: 差分bin值
- `--output_dir`: 输出目录
- `--all_types`: 生成所有差分类型的比较图

## 输出说明

生成的图片包含四个子图：
1. **SS比较**: 显示同符号配对的feeddown前后对比
2. **OS比较**: 显示异符号配对的feeddown前后对比
3. **Del比较**: 显示Del (OS-SS) 的feeddown前后对比
4. **相对变化**: 显示feeddown扣除的相对影响 (Feeddown/Original - 1)

## 数据路径

默认数据路径：
- 原始数据: `../dataset_merger/Merged/Proton/finalise_default.csv`
- Feeddown数据: `./Proton/finalise_feeddown_dispose_default.csv`

## 示例

```bash
# 生成delta的所有比较图
python run_feeddown_comparison.py

# 或者使用核心脚本生成特定图
python plot_feeddown_comparison.py \
    --original ../dataset_merger/Merged/Proton/finalise_default.csv \
    --feeddown ./Proton/finalise_feeddown_dispose_default.csv \
    --value_type gamma \
    --diff_type SPt \
    --diff_bin 1.5
```

## 注意事项

1. 确保数据文件存在且格式正确
2. 输出目录会自动创建
3. 图片以PDF格式保存，便于学术使用
4. 相对变化图有助于量化feeddown扣除的影响