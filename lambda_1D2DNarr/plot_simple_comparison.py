import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 设置样式
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12
sns.set_style("whitegrid")

# 定义要分析的文件和方法
METHOD_FILES = [
    {'filename': './LHC18q/Lambda/fit_obvs_2Dfit_default.csv', 'method': '2D Fit'},
    {'filename': './LHC18q/Lambda/fit_obvs_1Dfit_default.csv', 'method': '1D Fit'},
    {'filename': './LHC18q/Lambda/fit_obvs_narrMass_default.csv', 'method': 'Narrow Mass'}
]

# 更鲜艳的颜色方案
COLORS = {
    '2D Fit': '#FF2D2D',      # 鲜红色
    '1D Fit': '#2D7FFF',      # 鲜蓝色
    'Narrow Mass': '#00CC44'   # 鲜绿色
}

# 线型样式
LINE_STYLES = {
    '2D Fit': '-',
    '1D Fit': '--',
    'Narrow Mass': '-.'
}

# 标记样式
MARKERS = {
    '2D Fit': 'o',
    '1D Fit': 's',
    'Narrow Mass': '^'
}

def load_data():
    """加载数据"""
    all_data = []
    for entry in METHOD_FILES:
        try:
            df = pd.read_csv(entry['filename'])
            df['method'] = entry['method']
            # 筛选中心度 < 60的数据
            df = df[df['centrality'] < 60].copy()
            all_data.append(df)
            print(f"成功加载 {entry['filename']}: {len(df)} 行数据")
        except FileNotFoundError:
            print(f"警告: 文件 {entry['filename']} 未找到")
            continue
    
    if not all_data:
        raise FileNotFoundError("未找到任何数据文件！")
    
    return pd.concat(all_data, ignore_index=True)

def create_comparison_plot(data):
    """创建对比图"""
    # 指定要显示的粒子对类型
    pair_types = ['LambdaLambda', 'LambdaLambdaBar', 'LambdaBarLambdaBar']
    methods = data['method'].unique()
    
    # 创建2x3的子图布局
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 不设置总标题
    
    # 上排：Delta 参数
    for i, pair_type in enumerate(pair_types):
        ax = axes[0, i]
        pair_data = data[data['pair_type'] == pair_type]
        
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue
                
            centrality = method_data['centrality']
            color = COLORS[method]
            linestyle = LINE_STYLES[method]
            marker = MARKERS[method]
            
            ax.errorbar(centrality, method_data['delta_ss'], 
                       yerr=method_data['delta_ss_err'],
                       fmt=marker + linestyle, 
                       color=color, 
                       alpha=0.9, 
                       linewidth=3, 
                       markersize=8,
                       label=f'{method}', 
                       capsize=5,
                       capthick=2)
        
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel('Delta', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.4, linewidth=1)
        
        # 在每个子图添加图例
        ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # 下排：Gamma 参数
    for i, pair_type in enumerate(pair_types):
        ax = axes[1, i]
        pair_data = data[data['pair_type'] == pair_type]
        
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue
                
            centrality = method_data['centrality']
            color = COLORS[method]
            linestyle = LINE_STYLES[method]
            marker = MARKERS[method]
            
            ax.errorbar(centrality, method_data['rawgamma_ss'], 
                       yerr=method_data['rawgamma_ss_err'],
                       fmt=marker + linestyle, 
                       color=color, 
                       alpha=0.9, 
                       linewidth=3, 
                       markersize=8,
                       label=f'{method}', 
                       capsize=5,
                       capthick=2)
        
        ax.set_title(f'{pair_type}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Centrality (%)', fontsize=14)
        ax.set_ylabel('Gamma', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.4, linewidth=1)
        
        # 在每个子图添加图例
        ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig('lambda_correlation_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('lambda_correlation_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_two_panel_plot(data):
    """创建左右两个子图的对比图，每个子图都包含LambdaLambdaBar和LambdaBarLambdaBar"""
    # 指定要显示的粒子对类型
    pair_types = ['LambdaLambdaBar', 'LambdaBarLambdaBar']
    methods = data['method'].unique()
    
    # 创建1x2的子图布局
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图：Delta 参数
    ax = axes[0]
    for pair_type in pair_types:
        pair_data = data[data['pair_type'] == pair_type]
        
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue
                
            centrality = method_data['centrality']
            color = COLORS[method]
            linestyle = LINE_STYLES[method]
            marker = MARKERS[method]
            
            # 为不同的粒子对类型使用不同的标记填充
            if pair_type == 'LambdaLambdaBar':
                markerfacecolor = color
                markeredgecolor = color
                label_suffix = ' (ΛΛ̄)'
            else:  # LambdaBarLambdaBar
                markerfacecolor = 'white'
                markeredgecolor = color
                label_suffix = ' (Λ̄Λ̄)'
            
            ax.errorbar(centrality, method_data['delta_ss'], 
                       yerr=method_data['delta_ss_err'],
                       fmt=marker + linestyle, 
                       color=color,
                       markerfacecolor=markerfacecolor,
                       markeredgecolor=markeredgecolor,
                       markeredgewidth=2,
                       alpha=0.9, 
                       linewidth=3, 
                       markersize=8,
                       label=f'{method}{label_suffix}', 
                       capsize=5,
                       capthick=2)
    
    ax.set_title('Delta Parameters', fontsize=16, fontweight='bold')
    ax.set_xlabel('Centrality (%)', fontsize=14)
    ax.set_ylabel('Delta', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.4, linewidth=1)
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # 右图：Gamma 参数
    ax = axes[1]
    for pair_type in pair_types:
        pair_data = data[data['pair_type'] == pair_type]
        
        for method in methods:
            method_data = pair_data[pair_data['method'] == method]
            if len(method_data) == 0:
                continue
                
            centrality = method_data['centrality']
            color = COLORS[method]
            linestyle = LINE_STYLES[method]
            marker = MARKERS[method]
            
            # 为不同的粒子对类型使用不同的标记填充
            if pair_type == 'LambdaLambdaBar':
                markerfacecolor = color
                markeredgecolor = color
                label_suffix = ' (ΛΛ̄)'
            else:  # LambdaBarLambdaBar
                markerfacecolor = 'white'
                markeredgecolor = color
                label_suffix = ' (Λ̄Λ̄)'
            
            ax.errorbar(centrality, method_data['rawgamma_ss'], 
                       yerr=method_data['rawgamma_ss_err'],
                       fmt=marker + linestyle, 
                       color=color,
                       markerfacecolor=markerfacecolor,
                       markeredgecolor=markeredgecolor,
                       markeredgewidth=2,
                       alpha=0.9, 
                       linewidth=3, 
                       markersize=8,
                       label=f'{method}{label_suffix}', 
                       capsize=5,
                       capthick=2)
    
    ax.set_title('Gamma Parameters', fontsize=16, fontweight='bold')
    ax.set_xlabel('Centrality (%)', fontsize=14)
    ax.set_ylabel('Gamma', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.4, linewidth=1)
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig('lambda_correlation_two_panel.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('lambda_correlation_two_panel.png', dpi=300, bbox_inches='tight')
    plt.show()

def print_data_summary(data):
    """打印数据摘要"""
    print("\n" + "="*50)
    print("数据摘要")
    print("="*50)
    
    print(f"总数据行数: {len(data)}")
    print(f"分析方法: {', '.join(data['method'].unique())}")
    print(f"粒子对类型: {', '.join(data['pair_type'].unique())}")
    print(f"中心度范围: {data['centrality'].min():.1f} - {data['centrality'].max():.1f}")
    
    print("\n各方法数据分布:")
    for method in data['method'].unique():
        count = len(data[data['method'] == method])
        print(f"  {method}: {count} 行")

def main():
    """主函数"""
    try:
        print("正在加载数据...")
        data = load_data()
        print_data_summary(data)
        
        print("\n正在生成对比图...")
        create_comparison_plot(data)
        
        print("\n正在生成两面板对比图...")
        create_two_panel_plot(data)
        
        print("\n✓ 分析完成！")
        print("生成的文件:")
        print("- lambda_correlation_comparison.pdf")
        print("- lambda_correlation_comparison.png")
        print("- lambda_correlation_two_panel.pdf")
        print("- lambda_correlation_two_panel.png")
        
    except Exception as e:
        print(f"错误: {e}")
        print("\n请确保以下文件在当前目录下:")
        for entry in METHOD_FILES:
            print(f"- {entry['filename']}")
        
        import traceback
        print("\n详细错误信息:")
        traceback.print_exc()

if __name__ == "__main__":
    main()