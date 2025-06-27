import pandas as pd
import matplotlib.pyplot as plt

# 1. 画pT积分后的So/B vs centrality
ptint_csv = 'signal_sum_vs_centrality_inv_mass_fit_result_ptint.csv'
df_ptint = pd.read_csv(ptint_csv)

plt.figure(figsize=(8,6))
for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
    sub = df_ptint[df_ptint['particle'] == particle]
    plt.plot(sub['centrality'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)
plt.xlabel('Centrality (%)')
plt.ylabel('Signal/Background Ratio')
plt.title('Signal/Background Ratio vs Centrality (pT-integrated)')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('SoB_vs_centrality_ptint.pdf')
plt.close()

# 2. 画微分表（只画centrality=45）So/B vs pT
csv = 'signal_sum_vs_centrality_inv_mass_fit_result.csv'
df = pd.read_csv(csv)
plt.figure(figsize=(8,6))
for particle, marker, color in zip(['Lambda', 'LambdaBar'], ['o', 's'], ['C0', 'C1']):
    sub = df[(df['particle'] == particle) & (df['centrality'] == 35.0)]
    plt.plot(sub['pT_mean'], sub['signal_bkg_ratio'], marker=marker, color=color, label=particle, linewidth=2)
plt.xlabel('pT (GeV/c)')
plt.ylabel('Signal/Background Ratio')
plt.title('Signal/Background Ratio vs pT (Centrality=30-40%)')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('SoB_vs_pT_cent35.pdf')
plt.close()
