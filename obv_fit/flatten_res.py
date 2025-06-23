import ROOT
import pandas as pd
import re

input_file = "../rootfiles/Resolutions.root"
output_csv = "resolutions.csv"

f = ROOT.TFile.Open(input_file)
keys = [k.GetName() for k in f.GetListOfKeys()]
records = []

for key in keys:
    # 解析 plane 和 period
    m = re.match(r'hRes([A-Za-z0-9]+)(\d{2}[qr])', key)
    if not m:
        continue
    plane = m.group(1)
    period = m.group(2)
    h = f.Get(key)
    nbins = h.GetNbinsX()
    for i in range(1, nbins+1):
        centrality = h.GetXaxis().GetBinCenter(i)
        resolution = h.GetBinContent(i)
        resolution_err = h.GetBinError(i)
        records.append(dict(
            period=period,
            plane=plane,
            centrality=centrality,
            resolution=resolution,
            resolution_err=resolution_err
        ))

df = pd.DataFrame(records)
df.to_csv(output_csv, index=False)
print(f"已输出: {output_csv}")