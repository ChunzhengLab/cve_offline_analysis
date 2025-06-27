import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 常数定义
m_lambda = 1.115683
m_p = 0.938272
m_pi = 0.139570

def random_lorentz_vector(n, mass, p_range=(0.3, 2.0)):
    """生成 boost 后的 Λ 四动量"""
    p_mag = np.random.uniform(*p_range, size=n)
    theta = np.arccos(2 * np.random.rand(n) - 1)
    phi = 2 * np.pi * np.random.rand(n)
    px = p_mag * np.sin(theta) * np.cos(phi)
    py = p_mag * np.sin(theta) * np.sin(phi)
    pz = p_mag * np.cos(theta)
    E = np.sqrt(p_mag**2 + mass**2)
    return np.stack([E, px, py, pz], axis=1)

def two_body_decay_boosted(parent_lv, m1, m2):
    """两体衰变 + boost 回 lab 系"""
    n = parent_lv.shape[0]
    M = np.sqrt(np.maximum(parent_lv[:, 0]**2 - np.sum(parent_lv[:, 1:]**2, axis=1), 0))
    E1_cm = (M**2 + m1**2 - m2**2) / (2 * M)
    p_cm = np.sqrt(np.maximum(E1_cm**2 - m1**2, 0))

    theta = np.arccos(2 * np.random.rand(n) - 1)
    phi = 2 * np.pi * np.random.rand(n)
    px = p_cm * np.sin(theta) * np.cos(phi)
    py = p_cm * np.sin(theta) * np.sin(phi)
    pz = p_cm * np.cos(theta)
    p1_cm = np.stack([px, py, pz], axis=1)
    p2_cm = -p1_cm

    def boost(p_cm, m, parent_lv):
        E_cm = np.sqrt(np.sum(p_cm**2, axis=1) + m**2)
        beta = parent_lv[:, 1:] / parent_lv[:, 0:1]
        b2 = np.sum(beta**2, axis=1)
        gamma = 1.0 / np.sqrt(1 - b2)
        bp = np.sum(beta * p_cm, axis=1)
        gamma2 = (gamma - 1) / b2
        factor = gamma2[:, None] * bp[:, None]
        E_lab = gamma * E_cm + bp
        p_lab = p_cm + factor * beta + gamma[:, None] * beta * E_cm[:, None]
        return np.concatenate([E_lab[:, None], p_lab], axis=1)

    return boost(p1_cm, m1, parent_lv), boost(p2_cm, m2, parent_lv)

def smear(p4, mass, sigma=0.0005):
    p4_sm = p4.copy()
    p4_sm[:, 1:] += np.random.normal(0, sigma, size=(len(p4), 3))
    p_sq = np.sum(p4_sm[:, 1:]**2, axis=1)
    p4_sm[:, 0] = np.sqrt(p_sq + mass**2)
    return p4_sm

def combine_mass(p1, p2):
    E = p1[:, 0] + p2[:, 0]
    px = p1[:, 1] + p2[:, 1]
    py = p1[:, 2] + p2[:, 2]
    pz = p1[:, 3] + p2[:, 3]
    M2 = E**2 - px**2 - py**2 - pz**2
    return np.sqrt(np.maximum(M2, 0))

# 设置事件数
N = 1000
lambda_lvs = random_lorentz_vector(N, m_lambda)
p_lv, pi_lv = two_body_decay_boosted(lambda_lvs, m_p, m_pi)
p_lv = smear(p_lv, m_p)
pi_lv = smear(pi_lv, m_pi)

# 1. 真 Lambda 质量谱
true_mass = combine_mass(p_lv, pi_lv)
plt.hist(true_mass, bins=100, range=(1.10, 1.13), color='blue', alpha=0.7)
plt.xlabel("Mass (GeV)")
plt.ylabel("Counts")
plt.title("True Λ Mass (same p+π⁻)")
plt.tight_layout()
plt.savefig("true_lambda_mass.png")
plt.clf()

# 2. 背景谱（不同来源组合）
fake_mass_list = []
for i in range(N):
    for j in range(N):
        if i != j:
            m = combine_mass(p_lv[i:i+1], pi_lv[j:j+1])[0]
            fake_mass_list.append(m)
fake_mass = np.array(fake_mass_list)

plt.hist(fake_mass, bins=100, range=(1.10, 1.13), color='gray', alpha=0.7, label='Background')
plt.xlabel("Mass (GeV)")
plt.ylabel("Counts")
plt.title("Fake Λ Mass (p+π⁻ from different Λ)")
plt.tight_layout()
plt.savefig("fake_lambda_mass.png")
plt.clf()

# 3. Λ–Λ 2D 分布（排除 i==j）
mass1 = []
mass2 = []
for i in range(N):
    for j in range(N):
        if i != j:
            m1 = combine_mass(p_lv[i:i+1], pi_lv[i:i+1])[0]
            m2 = combine_mass(p_lv[j:j+1], pi_lv[j:j+1])[0]
            mass1.append(m1)
            mass2.append(m2)
mass1 = np.array(mass1)
mass2 = np.array(mass2)

plt.figure(figsize=(6, 5))
plt.hist2d(mass1, mass2, bins=100, range=[[1.10, 1.13], [1.10, 1.13]], cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel("Mass1 (GeV)")
plt.ylabel("Mass2 (GeV)")
plt.title("Λ–Λ 2D Mass (i ≠ j)")
plt.tight_layout()
plt.savefig("lambda_lambda_2D_mass.png")
plt.clf()

# 4. 投影到 mass1（横轴）
plt.hist(mass1, bins=100, range=(1.10, 1.13), color='purple', alpha=0.6, label='Projection from 2D')
plt.xlabel("Mass1 (GeV)")
plt.ylabel("Counts")
plt.title("Projection of Λ–Λ 2D onto Mass1")
plt.legend()
plt.tight_layout()
plt.savefig("lambda_lambda_mass1_projection.png")
plt.clf()

# 5. 与真实谱对比（mass1投影 + true mass）
plt.hist(mass1, bins=100, range=(1.10, 1.13), color='purple', alpha=0.6, label='Λ–Λ 2D Projection')
plt.hist(true_mass, bins=100, range=(1.10, 1.13), color='blue', alpha=0.5, label='True Λ')
plt.xlabel("Mass (GeV)")
plt.ylabel("Counts")
plt.title("Λ–Λ Projection vs True Λ")
plt.legend()
plt.tight_layout()
plt.savefig("lambda_lambda_projection_vs_true.png")
plt.clf()