import ROOT
import os

def extract_tpc_qvector(input_file):
    """Extracts Qx/Qy 3D histograms from input ROOT file and writes to a new file."""
    # 打开输入文件
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        print(f"Failed to open {input_file}")
        return

    lst = f.Get("default/ResultsList_default")
    if not lst:
        print(f"Cannot find ResultsList_default in {input_file}")
        return

    # 提取所需直方图
    h3qx = lst.FindObject("QxTPCRunCentVz")
    h3qy = lst.FindObject("QyTPCRunCentVz")
    if not h3qx or not h3qy:
        print(f"Missing Q-vector histograms in {input_file}")
        return

    # 自动生成输出文件名
    tag = os.path.basename(input_file).split("_")[-1].replace(".root", "")
    output_file = f"TPCQVectorMean{tag}.root"

    # 写出新文件
    fOut = ROOT.TFile(output_file, "RECREATE")
    lOut = ROOT.TList()
    lOut.SetName("fListTPC")
    lOut.SetOwner(True)
    lOut.Add(h3qx)
    lOut.Add(h3qy)
    lOut.Write("fListTPC", ROOT.TObject.kSingleKey)
    fOut.Close()

    print(f"Wrote output: {output_file}")

extract_tpc_qvector("data_18r.root")
extract_tpc_qvector("data_18q.root")
