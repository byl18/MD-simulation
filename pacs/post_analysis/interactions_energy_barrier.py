import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = ['Arial']  
# mpl.rcParams['font.weight'] = 'bold'  
# mpl.rcParams['font.size'] = 20
# 加载拓扑
gleap_target = md.load("../../2_md-preparation/gleap.target.gro")
atom_indices_ligand = gleap_target.topology.select("resid 33")
atom_indices_rna = gleap_target.topology.select("resid 0 to 32")
atom_index_ligand_min = min(atom_indices_ligand)

# 重新加载pacs_interactions
pacs_interactions = np.zeros((20, len(atom_indices_ligand), len(atom_indices_rna), 120))
for trial in range(20):
    with open(f"./pickle-cache/interactions_pi-stacking/pacs_interactions_{trial}.pickle", "rb") as f:
        obj = pickle.load(f)
        nonzero_indices = obj[0]
        nonzero_values = obj[1]
    pacs_interactions[trial][nonzero_indices] = nonzero_values




# 指定出口残基U12 (RNA残基索引)
exit_residue_index = 11  # U12
atom_indices_exit = gleap_target.topology.select(f"resid {exit_residue_index}")

print(f"Exit residue: {gleap_target.topology.residue(exit_residue_index)}")
print(f"Exit atom indices: {atom_indices_exit}")

# 从pacs_interactions中提取配体到U12的π-π堆积
# pacs_interactions shape: (20 trials, ligand_atoms, rna_atoms, 120 distance_bins)
# 提取U12对应的原子列
pacs_exit_interactions = pacs_interactions[:, :, atom_indices_exit, :]  # (20, ligand_atoms, n_exit_atoms, 120)

# 对所有配体原子和U12原子求和
total_exit_interactions = pacs_exit_interactions.sum(axis=(0, 1, 2))  # (120,)

# 距离轴 (0-12 Å, 步长0.1 Å)
distance_bins = np.arange(120) / 10

# 归一化
total_exit_interactions = total_exit_interactions 

# 计算能量 (PMF)
kBT = 0.596  # kcal/mol at 300K
E_exit = -kBT * np.log(total_exit_interactions + 1e-10)
E_exit = E_exit - E_exit[-1]  # 解离态（最远端）设为0

# 画图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# 左图：π-π堆积频率 vs 距离
ax1.plot(distance_bins, total_exit_interactions, 'r-', linewidth=2.5)
ax1.set_xlabel('Distance to U12 (Å)', fontsize=14)
ax1.set_ylabel('π-π Stacking Frequency', fontsize=14)
ax1.set_title(f'π-π Stacking: PRF: U12', fontsize=14)
ax1.grid(alpha=0.3)

# 右图：能量壁垒
ax2.plot(distance_bins, E_exit, 'b-', linewidth=2.5)
ax2.set_xlabel('Distance to U12 (Å)', fontsize=14)
ax2.set_ylabel('PMF (kcal/mol)', fontsize=14)
ax2.set_title(f'Exit Barrier via U12', fontsize=14)

# 标记关键位置
binding_site_idx = np.argmin(E_exit[:50])  # 结合位点 (短距离区域)
barrier_idx = np.argmax(E_exit[50:80]) + 50  # 出口壁垒 (5-8 Å区域)

ax2.axvline(x=distance_bins[binding_site_idx], color='green', linestyle='--', alpha=0.7, 
            label=f'Binding site ({distance_bins[binding_site_idx]:.1f} Å)')
ax2.axvline(x=distance_bins[barrier_idx], color='red', linestyle='--', alpha=0.7,
            label=f'Exit barrier ({distance_bins[barrier_idx]:.1f} Å)')
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# 标注能量值
binding_energy = E_exit[binding_site_idx]
barrier_height = E_exit[barrier_idx] - binding_energy
ax2.annotate(f'ΔG_bind = {binding_energy:.2f} kcal/mol',
             xy=(distance_bins[binding_site_idx], binding_energy),
             xytext=(distance_bins[binding_site_idx] - 2, binding_energy - 2.5),
             fontsize=11, color='green', weight='bold',
             arrowprops=dict(arrowstyle='->', color='green'))
ax2.annotate(f'E_barrier = {barrier_height:.2f} kcal/mol',
             xy=(distance_bins[barrier_idx], E_exit[barrier_idx]),
             xytext=(distance_bins[barrier_idx] + 1.5, E_exit[barrier_idx] - 1),
             fontsize=11, color='red', weight='bold',
             arrowprops=dict(arrowstyle='->', color='red'))

ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('exit_barrier_U12.png', dpi=300, bbox_inches='tight')
plt.show()
