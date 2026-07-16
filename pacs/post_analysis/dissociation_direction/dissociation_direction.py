import mdtraj as md
import numpy as np
from collections import Counter
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt 
import os, re

def count_cycles(trial_dir: str) -> int:
    pat = re.compile(r"^cycle(\d+)$")
    cycles = []
    for name in os.listdir(trial_dir):
        full = os.path.join(trial_dir, name)
        m = pat.match(name)
        if m and os.path.isdir(full):
            cycles.append(int(m.group(1)))
    return len(cycles)


traj = md.load("/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/trial001/prd.target.trjcat-all.pbc.skip10.xtc", top="/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/2_md-preparation/gleap.target.gro")

lig = traj.topology.select("resid 33 and not element H")
rna = traj.topology.select("resid 0 to 32 and not element H")

# 1) 临时距离 d0：ligand COM 到 RNA COM
com_lig = md.compute_center_of_mass(traj.atom_slice(lig))
com_rna = md.compute_center_of_mass(traj.atom_slice(rna))
d0 = np.linalg.norm(com_lig - com_rna, axis=1)  # nm

# 2) 选“bound 候选帧”：取距离最小的那一部分帧（比如最小 20%）
q = 20  # 可改 10/20/30
thr0 = np.percentile(d0, q)
bound_frames = np.where(d0 <= thr0)[0]
traj_b = traj[bound_frames]

# 3) 在 bound 帧里找 pocket 残基：ligand 周围 4.5 Å 内出现频率高的 RNA 残基
cutoff = 0.45  # nm (4.5 Å)
neighbors = md.compute_neighbors(traj_b, cutoff=cutoff, query_indices=lig, haystack_indices=rna)

res_count = Counter()
for atoms in neighbors:
    res_in_frame = {traj.topology.atom(a).residue.index for a in atoms}
    res_count.update(res_in_frame)

n = len(traj_b)
res_freq = {res: cnt/n for res, cnt in res_count.items()}

pocket_res = sorted([res for res, f in res_freq.items() if f >= 0.2])

print("Pocket residues (mdtraj residue.index):", pocket_res)
for r in pocket_res:
    residue = traj.topology.residue(r)
    print(r, residue.name, residue.resSeq)

# 4) 用 pocket 残基定义正式距离 d：ligand COM 到 pocket COM
sel = " or ".join([f"resid {r}" for r in pocket_res]) + " and not element H"
pocket_atoms = traj.topology.select(sel + " and not element H")

com_pocket = md.compute_center_of_mass(traj.atom_slice(pocket_atoms))
d = np.linalg.norm(com_lig - com_pocket, axis=1)  # 这就是你后续用的 d(t)






top = "/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/2_md-preparation/gleap.target.gro"
ref = md.load(top)
rna_fit = ref.topology.select("resid 0 to 32 and not element H")
lig_sel = "resid 33 and not element H"

# pocket_res 先用你跑出来的结果（最好用所有trial统计的）
pocket_res = pocket_res  # 你已有
sel = " or ".join([f"resid {r}" for r in pocket_res]) + " and not element H"

u_exits = []
for trial in range(1, 21):
    print("dealing" + str(trial))
    xtc = f"/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/trial{trial:03}/prd.target.trjcat-all.pbc.xtc"
    trial_dir = f"/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/trial{trial:03d}"
    n_cycle = count_cycles(trial_dir)
    start_frame = max(0, (n_cycle - 10) * 10000)   # ✅ 这里定义
    stride = 10
    start_in_strided = start_frame // stride
    traj = md.load(xtc, top=top, stride=stride)[start_in_strided:]
    traj.superpose(ref, atom_indices=rna_fit)
    lig = traj.topology.select(lig_sel)
    pocket_atoms = traj.topology.select(sel)
    com_lig = md.compute_center_of_mass(traj.atom_slice(lig))
    com_poc = md.compute_center_of_mass(traj.atom_slice(pocket_atoms))
    r = com_lig - com_poc
    d = np.linalg.norm(r, axis=1)
    u = r / d[:, None]
    bound_th, unbound_th = 1.0, 2.5
    bound_idx = np.where(d <= bound_th)[0]
    unbound_idx = np.where(d >= unbound_th)[0]
    if len(bound_idx)==0 or len(unbound_idx)==0:
        print("trial", trial, "no full transition")
        u_exits.append([np.nan, np.nan, np.nan])
        continue
    start = bound_idx[0]
    end_candidates = unbound_idx[unbound_idx > start]
    if len(end_candidates)==0:
        print("trial", trial, "no unbound after bound")
        u_exits.append([np.nan, np.nan, np.nan])
        continue
    end = end_candidates[0]
    seg = slice(start, end+1)
    d_exit, width = 1.8, 0.1
    mask = (d[seg] >= d_exit-width) & (d[seg] <= d_exit+width)
    if mask.sum() < 3:
        k = np.argmin(np.abs(d[seg] - d_exit)) + start
        u_exit = u[k]
    else:
        u_exit = u[seg][mask].mean(axis=0)
        u_exit = u_exit / np.linalg.norm(u_exit)
    u_exits.append(u_exit)

u_exits = np.array(u_exits)  # 转为 numpy 数组
np.save("/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/post_analysis/dissociation_direction/u_exits.npy", u_exits)
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['Nimbus Roman']  
mpl.rcParams['font.weight'] = 'bold'  
mpl.rcParams['font.size'] = 20
U = np.array(u_exits)
valid = ~np.isnan(U).any(axis=1)
U2 = U[valid]

# 聚类：先试 2~4 类
labels = KMeans(n_clusters=2, random_state=0).fit_predict(U2)
colors = np.where(labels == 1, "red", "blue")  # 1红，0蓝

# 画角度散点
phi   = np.arctan2(U2[:,1], U2[:,0])
theta = np.arccos(np.clip(U2[:,2], -1, 1))

plt.figure(figsize=(6, 5))
plt.scatter(phi, theta, c=colors)
for i, (p,t) in enumerate(zip(phi,theta)):
    plt.text(p, t, str(np.where(valid)[0][i]+1), fontsize=12)

plt.xlabel("phi (azimuth)",fontweight='bold')
plt.ylabel("theta (polar)",fontweight='bold')
plt.title("Exit direction clustering (10 trials)",fontweight='bold',pad=20)
plt.tight_layout()
plt.savefig('/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/post_analysis/figures/dissociation_direction.png',dpi=200)


########################################################
#####################TwoPCA#############################
########################################################
'''
sel = " or ".join([f"resid {r}" for r in pocket_res]) + " and not element H"
pocket_atoms = ref.topology.select(sel)

# pocket COM（nm）
pocket_com = md.compute_center_of_mass(ref.atom_slice(pocket_atoms))[0]  # (3,)

# RNA 坐标（nm -> Å 更直观）
xyz = ref.xyz[0] * 10.0
p0 = pocket_com * 10.0  # Å

# 两个代表方向（单位向量），箭头长度（Å）
L = 24.0   # 30 Å = 3 nm
uA_ = uA / np.linalg.norm(uA)
uB_ = uB / np.linalg.norm(uB)

p0 = pocket_com * 10.0
pA = p0 + (uA_ * L)
pB = p0 + (uB_ * L)

pml = f"""
load ref_target.pdb, ref
hide everything, ref
show cartoon, ref and resi 1-33
color gray80, ref

# 起点/终点 pseudoatoms
pseudoatom p0, pos=[{p0[0]:.3f},{p0[1]:.3f},{p0[2]:.3f}]
pseudoatom pA, pos=[{pA[0]:.3f},{pA[1]:.3f},{pA[2]:.3f}]
pseudoatom pB, pos=[{pB[0]:.3f},{pB[1]:.3f},{pB[2]:.3f}]

# 用 distance 先画“箭杆”（简单稳定）
distance dirA, p0, pA
distance dirB, p0, pB
set dash_width, 4
set dash_gap, 0
color red, dirA
color blue, dirB

# 让球点明显一些
show spheres, p0 pA pB
set sphere_scale, 0.4
color yellow, p0
color red, pA
color blue, pB

bg_color white
"""
open("draw_arrows.pml","w").write(pml)
print("Wrote draw_arrows.pml")
'''

########################################################
#######################AllPCA###########################
########################################################
sel = " or ".join([f"resid {r}" for r in pocket_res]) + " and not element H"
pocket_atoms = ref.topology.select(sel)
pocket_com = md.compute_center_of_mass(ref.atom_slice(pocket_atoms))[0]  # nm

# nm -> Å
p0 = pocket_com * 10.0

# ===== 2) 10 个方向向量 U（单位化），箭头长度 =====
U = np.asarray(U)          # (10,3)
labels = np.asarray(labels)  # (10,)

# 单位化，避免数值问题
U_norm = np.linalg.norm(U, axis=1)
U_unit = U / U_norm[:, None]

L = 24.0  # Å

# ===== 3) 生成 pml =====
pml_lines = []
pml_lines.append("load ref_target.pdb, ref")
pml_lines.append("hide everything, ref")
pml_lines.append("show cartoon, ref")
pml_lines.append("color gray80, ref")
pml_lines.append("bg_color white")
pml_lines.append("")
pml_lines.append("# pocket center")
pml_lines.append(f"pseudoatom p0, pos=[{p0[0]:.3f},{p0[1]:.3f},{p0[2]:.3f}]")
pml_lines.append("show spheres, p0")
pml_lines.append("set sphere_scale, 0.5, p0")
pml_lines.append("color yellow, p0")
pml_lines.append("")
pml_lines.append("# arrows for 10 trials")
pml_lines.append("set dash_width, 4")
pml_lines.append("set dash_gap, 0")
pml_lines.append("set label_size, 20")


for i in range(U_unit.shape[0]):
    color = "red" if labels[i] == 1 else "blue"
    p_end = p0 + U_unit[i] * L
    name_end = f"pEnd{i+1:02d}"
    name_line = f"dir{i+1:02d}"
    pml_lines.append(f"pseudoatom {name_end}, pos=[{p_end[0]:.3f},{p_end[1]:.3f},{p_end[2]:.3f}]")
    pml_lines.append(f"distance {name_line}, p0, {name_end}")
    pml_lines.append(f"color {color}, {name_line}")
    pml_lines.append(f"show spheres, {name_end}")
    pml_lines.append(f"set sphere_scale, 0.25, {name_end}")
    pml_lines.append(f"color {color}, {name_end}")
    pml_lines.append(f"label {name_end}, \"{i+1}\"")

pml_lines.append("set label_distance_digits, 0")
pml_lines.append("hide labels, dir*")
pml = "\n".join(pml_lines) + "\n"
open("draw_arrows_all.pml", "w").write(pml)
print("Wrote draw_arrows_all.pml")
