#!/usr/bin/env python3
import argparse
import os
import sys
import re
import subprocess
import ROOT
import glob
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional

def extract_info_from_filename(filename: str) -> Tuple[str, str]:
    """从ROOT文件名提取数据集和粒子类型"""
    # 提取数据集(18q或18r)
    dataset_match = re.search(r'18[qr]', filename)
    dataset = f"LHC{dataset_match.group(0)}" if dataset_match else "unknown"

    # 提取粒子类型
    particle_types = ["Proton", "Pion", "Hadron"]
    particle_type = "unknown"
    for pt in particle_types:
        if pt in filename:
            particle_type = pt
            break

    return dataset, particle_type

def get_tasks_from_root(root_file: str) -> List[str]:
    """从ROOT文件中提取可用的任务列表"""
    try:
        f = ROOT.TFile.Open(root_file)
        if not f or f.IsZombie():
            print(f"错误：无法打开ROOT文件：{root_file}")
            return []

        # 获取所有的目录名作为task
        tasks = []
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if isinstance(obj, ROOT.TDirectoryFile):
                task_name = obj.GetName()
                # 检查是否有 ListResults_{task_name}
                if obj.Get(f"ListResults_{task_name}"):
                    tasks.append(task_name)

        f.Close()
        return tasks
    except Exception as e:
        print(f"获取任务列表时出错：{str(e)}")
        return []

def print_root_info(root_file: str) -> None:
    """打印ROOT文件的基本信息"""
    try:
        f = ROOT.TFile.Open(root_file)
        if not f or f.IsZombie():
            print(f"错误：无法打开ROOT文件：{root_file}")
            return

        # 基本文件信息
        file_size = os.path.getsize(root_file) / (1024 * 1024)  # MB

        # 从文件名提取信息
        filename = os.path.basename(root_file)
        dataset, particle_type = extract_info_from_filename(filename)

        # 获取任务列表
        tasks = get_tasks_from_root(root_file)

        # 打印信息
        print("\n========== ROOT文件信息 ==========")
        print(f"文件路径: {root_file}")
        print(f"文件大小: {file_size:.2f} MB")
        print(f"推断数据集: {dataset}")
        print(f"推断粒子类型: {particle_type}")
        print(f"可用任务数: {len(tasks)}")
        print(f"可用任务列表: {', '.join(tasks)}")
        print("=================================\n")

        f.Close()
    except Exception as e:
        print(f"获取ROOT文件信息时出错：{str(e)}")

def run_command(cmd: List[str], task_name: str = "") -> bool:
    """执行命令并打印输出"""
    cmd_str = " ".join(cmd)
    prefix = f"[{task_name}] " if task_name else ""
    print(f"\n{prefix}执行命令: {cmd_str}")

    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        # 实时打印输出
        for line in process.stdout:
            print(f"{prefix}{line.strip()}")

        process.wait()

        if process.returncode != 0:
            print(f"{prefix}命令执行失败，返回码: {process.returncode}")
            return False

        return True
    except Exception as e:
        print(f"{prefix}执行命令时出错: {str(e)}")
        return False

def run_analysis_pipeline(
    input_root: str,
    output_dir: str,
    task: str,
    skip_flatten: bool = False,
    skip_inv_mass: bool = False,
    skip_obvs: bool = False,
    skip_finalise: bool = False,
    flatten_mode: str = "eff_cali",
    fs_force_non_neg: bool = True,
    fs_min: Optional[float] = None,
    fs_max: Optional[float] = None,
    range_mass: float = 0.0014,
    no_fit: bool = False,
    reso_file: str = "./resolutions.csv"
) -> bool:
    """执行完整的分析流程

    Args:
        range_mass: 质量窗口范围半宽度，表示中心质量左右各±该值（默认0.0014 GeV/c²）
    """
    # 构建基本参数
    filename = os.path.basename(input_root)
    dataset, particle_type = extract_info_from_filename(filename)

    if not dataset or not particle_type:
        print(f"错误：无法从文件名 {filename} 提取数据集或粒子类型信息")
        return False

    # 创建输出目录
    task_output_dir = os.path.join(output_dir, dataset, particle_type)
    os.makedirs(task_output_dir, exist_ok=True)

    success = True
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. flatten_data.py
    if not skip_flatten:
        flatten_cmd = [
            sys.executable,
            os.path.join(script_dir, "flatten_data.py"),
            "-i", input_root,
            "-t", task,
            "-o", output_dir,
            "-m", flatten_mode
        ]
        if not run_command(flatten_cmd, f"{task}/flatten"):
            print(f"[{task}] flatten_data.py 执行失败")
            success = False
            if not (skip_inv_mass and skip_obvs and skip_finalise):
                print(f"[{task}] 由于flatten_data.py失败，后续步骤将被跳过")
                return False

    # 2. fit_inv_mass.py
    if not skip_inv_mass and success:
        inv_mass_cmd = [
            sys.executable,
            os.path.join(script_dir, "fit_inv_mass.py"),
            "-i", output_dir,
            "-d", dataset,
            "-p", particle_type,
            "-t", task,
            "-o", output_dir
        ]

        # 默认不启用fs-force-non-neg
        if fs_force_non_neg:
            inv_mass_cmd.append("--fs-force-non-neg")

        if not run_command(inv_mass_cmd, f"{task}/fit_inv_mass"):
            print(f"[{task}] fit_inv_mass.py 执行失败")
            success = False
            if not (skip_obvs and skip_finalise):
                print(f"[{task}] 由于fit_inv_mass.py失败，后续步骤将被跳过")
                return False

    # 3. fit_obvs.py
    if not skip_obvs and success:
        obvs_cmd = [
            sys.executable,
            os.path.join(script_dir, "fit_obvs.py"),
            "-i", output_dir,
            "-d", dataset,
            "-p", particle_type,
            "-t", task,
            "-o", output_dir
        ]

        if fs_min is not None:
            obvs_cmd.extend(["--fs-min", str(fs_min)])
        if fs_max is not None:
            obvs_cmd.extend(["--fs-max", str(fs_max)])

        # 添加质量窗口参数
        obvs_cmd.extend(["--range-mass", str(range_mass)])
        if no_fit:
            obvs_cmd.append("--no-fit")

        if not run_command(obvs_cmd, f"{task}/fit_obvs"):
            print(f"[{task}] fit_obvs.py 执行失败")
            success = False
            if not skip_finalise:
                print(f"[{task}] 由于fit_obvs.py失败，后续步骤将被跳过")
                return False

    # 4. finalise.py
    if not skip_finalise and success:
        finalise_cmd = [
            sys.executable,
            os.path.join(script_dir, "finalise.py"),
            "-i", output_dir,
            "-d", dataset,
            "-p", particle_type,
            "-t", task,
            "-o", output_dir,
            "-r", reso_file
        ]

        if not run_command(finalise_cmd, f"{task}/finalise"):
            print(f"[{task}] finalise.py 执行失败")
            success = False
            return False

    print(f"\n[{task}] 分析流程{'全部' if success else '部分'}完成")
    return success

def main():
    """主函数，解析命令行参数并执行相应操作"""
    parser = argparse.ArgumentParser(
        description="ALICE CVE offline analysis tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 基本参数
    parser.add_argument("-i", "--input-root", type=str, default="../train_output/AnalysisResults_CVE2025_18q_TPC_Proton.root",
                       help="输入ROOT文件路径")
    parser.add_argument("-o", "--output-dir", type=str, default="./",
                       help="输出目录")
    parser.add_argument("-s", "--show-tasks", action="store_true",
                       help="只显示可用任务，不执行分析")

    # 分析模式
    parser.add_argument("-m", "--mode", type=str, choices=["default", "all", "list"], default="default",
                       help="分析模式：default(只分析default任务), all(分析所有任务), list(分析指定任务列表)")
    parser.add_argument("-l", "--task-list", type=str, nargs="+",
                       help="要分析的任务列表，仅当mode=list时有效")

    # 步骤控制
    parser.add_argument("--skip-flatten", action="store_true",
                       help="跳过flatten_data.py步骤")
    parser.add_argument("--skip-inv-mass", action="store_true",
                       help="跳过fit_inv_mass.py步骤")
    parser.add_argument("--skip-obvs", action="store_true",
                       help="跳过fit_obvs.py步骤")
    parser.add_argument("--skip-finalise", action="store_true",
                       help="跳过finalise.py步骤")

    # 分析参数
    parser.add_argument("--flatten-mode", type=str, default="eff_cali",
                       choices=["raw_mass", "eff_cali"],
                       help="扁平化模式：raw_mass-使用原始质量分布，eff_cali-使用delta3的bin entries作为计数")
    parser.add_argument("--fs-force-non-neg", action="store_true", default=False,
                       help="强制信号分数非负（默认不启用）")
    parser.add_argument("--fs-min", type=float, default=None,
                       help="信号分数最小值过滤")
    parser.add_argument("--fs-max", type=float, default=None,
                       help="信号分数最大值过滤")
    parser.add_argument("--range-mass", type=float, default=0.0014,
                       help="质量窗口范围半宽度，表示中心质量左右各±该值 (默认: 0.0014 GeV/c²)")
    parser.add_argument("--no-fit", action="store_true",
                       help="不进行拟合，直接在质量窗口范围内计算加权平均值")
    parser.add_argument("-r", "--reso-file", type=str, default="./resolutions.csv",
                       help="分辨率文件路径")

    args = parser.parse_args()

    # 检查输入文件
    if not os.path.exists(args.input_root):
        print(f"错误：输入文件 {args.input_root} 不存在")
        return 1

    # 打印ROOT文件信息
    print_root_info(args.input_root)

    # 获取任务列表
    tasks = get_tasks_from_root(args.input_root)

    # 如果只是显示任务，则退出
    if args.show_tasks:
        print("\n可用任务列表：")
        for i, task in enumerate(tasks, 1):
            print(f"{i}. {task}")
        return 0

    # 检查mode和task_list参数
    if args.mode == "list" and not args.task_list:
        print("错误：mode=list时必须使用--task-list指定要分析的任务")
        return 1

    # 确定要分析的任务
    tasks_to_analyze = []
    if args.mode == "default":
        if "default" in tasks:
            tasks_to_analyze = ["default"]
        else:
            print("警告：未找到'default'任务，将使用第一个可用任务")
            tasks_to_analyze = [tasks[0]] if tasks else []
    elif args.mode == "all":
        tasks_to_analyze = tasks
    elif args.mode == "list":
        # 验证所有指定的任务都存在
        for task in args.task_list:
            if task not in tasks:
                print(f"警告：指定的任务 '{task}' 在ROOT文件中不存在")
        # 过滤有效任务
        tasks_to_analyze = [task for task in args.task_list if task in tasks]

    if not tasks_to_analyze:
        print("错误：没有找到可分析的任务")
        return 1

    print(f"\n将分析以下任务: {', '.join(tasks_to_analyze)}")

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    # 执行分析
    results = {}
    for task in tasks_to_analyze:
        print(f"\n{'='*30}")
        print(f"开始分析任务: {task}")
        print(f"{'='*30}")

        success = run_analysis_pipeline(
            input_root=args.input_root,
            output_dir=args.output_dir,
            task=task,
            skip_flatten=args.skip_flatten,
            skip_inv_mass=args.skip_inv_mass,
            skip_obvs=args.skip_obvs,
            skip_finalise=args.skip_finalise,
            flatten_mode=args.flatten_mode,
            fs_force_non_neg=args.fs_force_non_neg,
            fs_min=args.fs_min,
            fs_max=args.fs_max,
            range_mass=args.range_mass,
            no_fit=args.no_fit,
            reso_file=args.reso_file
        )

        results[task] = success

    # 打印总结
    print("\n\n=========== 分析结果总结 ===========")
    all_success = True
    for task, success in results.items():
        status = "成功" if success else "失败"
        print(f"任务 {task}: {status}")
        if not success:
            all_success = False

    return 0 if all_success else 1

if __name__ == "__main__":
    sys.exit(main())
