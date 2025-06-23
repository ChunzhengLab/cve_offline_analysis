import ROOT
import numpy as np
import matplotlib.pyplot as plt

input_root = "../train_output/AnalysisResults_CVE2025_18q_TPC_Proton.root"

PAIR_TYPE = "LambdaProton"    # You can set to others
DIFF_TYPE = "Intg"             # "SPt", "DEta", "Intg"

PAIR_MAP = {
    "LambdaProton": "0",
    "LambdaProtonBar": "1",
    "LambdaBarProton": "2",
    "LambdaBarProtonBar": "3"
}

def find_obj(lst, name_key):
    for i in range(lst.GetSize()):
        obj = lst.At(i)
        name = obj.GetName()
        if name_key in name:
            return obj
    raise RuntimeError(f"Cannot find object containing '{name_key}'")

def get_hist_and_profile(lst, pair_type, diff_type):
    suffix = PAIR_MAP[pair_type]
    tried = []
    hist = None
    for prefix in ["fHist3", "fHist2"]:
        hname = f"{prefix}{pair_type}Mass{diff_type}_{suffix}"
        tried.append(hname)
        try:
            hist = find_obj(lst, hname)
            break
        except RuntimeError:
            hist = None
    if hist is None:
        raise RuntimeError(f"Tried: {tried}, but found no matching hist object!")
    # Profile名字永远是3D
    profname = f"fProfile3DDelta{pair_type}Mass{diff_type}_{suffix}"
    profile = find_obj(lst, profname)
    return hist, profile

# Open ROOT file
f = ROOT.TFile.Open(input_root)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open ROOT file: {input_root}")
default_dir = f.Get("default")
lst = default_dir.Get("ListResults_default")

h3, p3d = get_hist_and_profile(lst, PAIR_TYPE, DIFF_TYPE)

nx, ny, nz = h3.GetNbinsX(), h3.GetNbinsY(), h3.GetNbinsZ()
invmass_bins = np.array([h3.GetZaxis().GetBinCenter(iz) for iz in range(1, nz+1)])

hist_counts = np.zeros(nz)
entry_counts = np.zeros(nz)

for ix in range(1, nx+1):
    for iy in range(1, ny+1):
        for iz in range(1, nz+1):
            hist_counts[iz-1] += h3.GetBinContent(ix, iy, iz)
            entry_counts[iz-1] += p3d.GetBinEntries(p3d.GetBin(ix, iy, iz))

# Normalize
def norm(arr):
    s = np.sum(arr)
    return arr / s if s > 0 else arr

hist_counts_norm = norm(hist_counts)
entry_counts_norm = norm(entry_counts)

def mean_sigma(x, y):
    mean = np.sum(x * y)
    var = np.sum(((x - mean)**2) * y)
    sigma = np.sqrt(var)
    return mean, sigma

hist_mean, hist_sigma = mean_sigma(invmass_bins, hist_counts_norm)
entry_mean, entry_sigma = mean_sigma(invmass_bins, entry_counts_norm)

# Plot
plt.figure(figsize=(8,5))
plt.plot(invmass_bins, hist_counts_norm, 'o-', color='red', label='TH (Histogram counts)')
plt.plot(invmass_bins, entry_counts_norm, 's--', color='blue', label='TProfile (Bin entries)')

# Add mean/sigma as a text box
textstr = (
    f"Inv Mass:        mean = {hist_mean:.6f}\n"
    f"                 sigma = {hist_sigma:.6f}\n"
    f"Eff. Calibrated: mean = {entry_mean:.6f}\n"
    f"                 sigma = {entry_sigma:.6f}"
)
plt.gca().text(
    0.98, 0.98, textstr,
    fontsize=10, ha='right', va='top',
    transform=plt.gca().transAxes,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray')
)

plt.xlabel("Invariant Mass (GeV/c²)")
plt.ylabel("Normalized counts")
plt.title(f"{PAIR_TYPE}, {DIFF_TYPE} (Integrated, normalized)")
plt.legend()
plt.grid(True)
plt.tight_layout()
output_pdf = f"output_{PAIR_TYPE}_{DIFF_TYPE}.pdf"
plt.savefig(output_pdf)
print(f"Saved to {output_pdf}")
