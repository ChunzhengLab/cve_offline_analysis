import ROOT
import os

# 输入输出文件路径
file1_path = "../efficiency_study/scripts/eff_pt_calib_cent.root"  # 包含 fListNUE
file2_path = "nua.root"  # 包含 TH2D (nua_corrections 目录中)
output_path = "eff_pt_calib_cent.root"

# 支持多个标签
tags = ["18q", "18r"]  # 直接指定要处理的标签
print(f"将处理以下标签: {tags}")

# 粒子类型和中心度配置
particle_list = ["poshadron", "neghadron", "proton", "antiproton", "lambda", "antilambda"]
cent_bins = list(range(7))  # 0-6 共7个中心度区间
print(f"使用标签: {tags}")

# 标记是否只有cent0
only_cent0 = True


# ---------- 1. 处理第一个文件 ----------
f1 = ROOT.TFile.Open(file1_path, "READ")
fListNUE = f1.Get("fListNUE")
if not fListNUE or not fListNUE.InheritsFrom("TList"):
    raise RuntimeError("找不到 TList 'fListNUE'")

graphs = []
for obj in fListNUE:
    if obj.InheritsFrom("TGraphErrors"):
        g = ROOT.TGraph(obj)
        g.SetName(obj.GetName())
        graphs.append(g)
f1.Close()                       # 这里关掉没问题，graph 已完全脱离

# ---------- 2. 处理第二个文件 ----------
f2 = ROOT.TFile.Open(file2_path, "READ")
hists = []

# 获取 nua_corrections 目录
nua_dir = f2.Get("nua_corrections")
if not nua_dir or not nua_dir.InheritsFrom("TDirectory"):
    print(f"警告: 在 {file2_path} 中找不到目录 'nua_corrections'")
    print("尝试直接从根目录读取...")
    nua_dir = f2  # 如果没有 nua_corrections 目录，就直接从根目录读取

# 统计处理的直方图数量
hist_counts = {"poshadron": 0, "neghadron": 0, "proton": 0, "antiproton": 0, "lambda": 0, "antilambda": 0, "其他": 0}
cent_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
print("开始读取NUA直方图...")

processed_hists = set()  # 用于跟踪已处理的直方图名称

# 首先，尝试从nua_corrections目录获取所有中心度的直方图
# 收集所有找到的cent0直方图，用于后续可能的复制
cent0_hists = {}

for tag in tags:  # 遍历所有标签
    for p in particle_list:
        for c in cent_bins:
            # 构造期望的直方图名称（根据make_nua.py的命名模式）
            hist_name = f"nua_pt_{p}_{tag}_cent{c}"
            obj = nua_dir.Get(hist_name)
            
            if obj and obj.InheritsFrom("TH2"):
                print(f"找到直方图: {hist_name}")
                ROOT.gROOT.cd()
                h = obj
                if not h:
                    print(f"❌ Clone 失败: {hist_name}"); continue
                h.SetName(hist_name)
                h.SetDirectory(0)
                hists.append(h)
                processed_hists.add(hist_name)
                
                # 更新统计
                hist_counts[p] += 1
                cent_counts[c] += 1
                
                # 如果是cent0直方图，保存起来供后续可能的复制
                if c == 0:
                    cent0_hists[(p, tag)] = h  # 使用元组(p, tag)作为键
                elif only_cent0:
                    # 如果找到了非cent0的直方图，则不仅仅只有cent0
                    only_cent0 = False
            else:
                print(f"警告: 未找到直方图 {hist_name}")

# 如果需要，也可以遍历目录中的所有对象作为备份
print("\n开始扫描目录中的所有直方图...")
for key in nua_dir.GetListOfKeys():
    obj = key.ReadObj()
    name = obj.GetName()
    
    # 如果已经处理过这个直方图，则跳过
    if name in processed_hists:
        continue
        
    if obj.InheritsFrom("TH2"):
        print(f"处理额外直方图: {name}")
        ROOT.gROOT.cd()
        h = obj
        if not h:
            print(f"❌ Clone 失败: {name}"); continue
        h.SetName(name)
        h.SetDirectory(0)
        hists.append(h)
        processed_hists.add(name)
        
        # 统计粒子类型
        counted = False
        for particle in hist_counts.keys():
            if particle in name:
                hist_counts[particle] += 1
                counted = True
                break
        if not counted:
            hist_counts["其他"] += 1
            
        # 统计中心度区间
        for cent in cent_counts.keys():
            if f"_cent{cent}" in name:
                cent_counts[cent] += 1
                break

# 检查是否只有cent0，如果是，复制cent0直方图到其他中心度
if only_cent0 and len(cent0_hists) > 0:
    print("\n检测到只有中心度0的直方图，将复制到其他中心度...")
    
    for (p, tag), h_cent0 in cent0_hists.items():
        for c in cent_bins[1:]:  # 跳过cent0
            # 构造新直方图名称
            new_hist_name = f"nua_pt_{p}_{tag}_cent{c}"
            
            # 如果已经处理过这个直方图，则跳过
            if new_hist_name in processed_hists:
                continue
                
            print(f"复制 cent0 → cent{c} 直方图: {new_hist_name}")
            
            # 创建新直方图（复制cent0的）
            ROOT.gROOT.cd()
            h_new = h_cent0.Clone(new_hist_name)
            h_new.SetTitle(h_new.GetTitle().replace("Centrality 0-7%", f"Centrality {c}-{c+1}%"))
            h_new.SetDirectory(0)
            
            # 添加到列表
            hists.append(h_new)
            processed_hists.add(new_hist_name)
            
            # 更新统计
            hist_counts[p] += 1
            cent_counts[c] += 1

# 总结处理结果
print("\n处理结果汇总:")
expected_total = len(particle_list) * len(cent_bins) * len(tags)
actual_total = len(processed_hists)
print(f"应处理直方图总数: {expected_total}")
print(f"实际处理直方图总数: {actual_total}")

if actual_total < expected_total:
    print("\n缺失的直方图:")
    for p in particle_list:
        for tag in tags:
            for c in cent_bins:
                hist_name = f"nua_pt_{p}_{tag}_cent{c}"
                if hist_name not in processed_hists:
                    print(f"  - {hist_name}")

# 现在可以马上关 file2 —— 直方图已脱离
f2.Close()

# ---------- 3. 写入新文件 ----------
out = ROOT.TFile.Open(output_path, "RECREATE")
fListNUENUA = ROOT.TList()
fListNUENUA.SetName("fListNUENUA")

# 先添加TGraph (NUE)
print(f"\n添加 {len(graphs)} 个TGraph (NUE)到输出文件")
for g in graphs:
    fListNUENUA.Add(g)

# 再添加TH2 (NUA)
print(f"添加 {len(hists)} 个TH2 (NUA)到输出文件")
if len(hists) == 0:
    print("警告: 没有找到任何NUA直方图! 检查输入文件和目录结构。")
for h in hists:
    fListNUENUA.Add(h)

# 写入输出文件
fListNUENUA.Write("fListNUENUA", ROOT.TObject.kSingleKey)
out.Close()

print(f"完成：{len(graphs)} 个 TGraph + {len(hists)} 个 TH2 已写入 {output_path}")
print("\n各粒子类型的NUA直方图统计:")
for particle, count in hist_counts.items():
    if count > 0:
        print(f"  - {particle}: {count} 个直方图")

print("\n各中心度区间的NUA直方图统计:")
for cent, count in cent_counts.items():
    print(f"  - 中心度{cent}: {count} 个直方图")
    # 检查每个中心度是否有期望数量的直方图（每个粒子类型一个）
    expected = len(particle_list) * len(tags)  # 应该有6个粒子类型 * 2个标签
    if count != expected:
        print(f"    警告: 中心度{cent}只有{count}个直方图，应该有{expected}个")
        # 列出此中心度缺少的粒子类型
        for p in particle_list:
            for tag in tags:
                hist_name = f"nua_pt_{p}_{tag}_cent{cent}"
                if hist_name not in processed_hists:
                    print(f"      缺少: {p} 和 {tag}")

print(f"\n输出文件大小: {os.path.getsize(output_path)/1024/1024:.2f} MB")
