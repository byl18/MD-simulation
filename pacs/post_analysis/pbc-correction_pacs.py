import os
import mdtraj as md
import numpy as np
import pickle
import matplotlib.pyplot as plt
import concurrent.futures


N_TIMESTEP_IN_REPLICA = 201
N_REPLICAS_IN_CYCLE = 64
CUTOFF_DISTANCE = 4

BASE_DIR = "/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2"

gleap_target = md.load(f"{BASE_DIR}/../2_md-preparation/gleap.target.gro")
atom_indices_rna    = gleap_target.topology.select("resid  0 to 32")
atom_indices_ligand = gleap_target.topology.select("resid 33")


def worker(trial):
    new_xtc = md.load(f"{BASE_DIR}/trial{str(trial+1).zfill(3)}/prd.target.trjcat-all.pbc.xtc", top=gleap_target)
    new_xtc = new_xtc.image_molecules()
    new_xtc = new_xtc.superpose(gleap_target, atom_indices=atom_indices_rna, ref_atom_indices=atom_indices_rna)
    new_coms_ligand  = md.compute_center_of_mass(new_xtc.atom_slice(atom_indices_ligand))
    new_coms_rna      = md.compute_center_of_mass(new_xtc.atom_slice(atom_indices_rna))
    new_com_vectors   = new_coms_rna - new_coms_ligand
    new_com_distances = np.linalg.norm(new_com_vectors, axis=1)
    pickle.dump(new_com_vectors,   open(f"./pickle-cache/com-features/pacs_com_vectors_{trial}.pickle",   "wb"))
    pickle.dump(new_com_distances, open(f"./pickle-cache/com-features/pacs_com_distances_{trial}.pickle", "wb"))
    print(trial, "new frames", new_xtc.n_frames)
    plt.plot(new_com_distances, linewidth=1, label=f"new_com_distances: {len(new_com_distances)}")
    plt.title(trial)
    plt.legend()
    plt.savefig('/gs/bs/tga-KitaoLab-NexusYtec/yilan/projects/PreQ1/trial01/pacs2/post_analysis/figures/' + str(trial+1) + '.png',dpi=200)
    plt.close()


if not os.path.exists("./pickle-cache/com-features"):
    os.mkdir("./pickle-cache/com-features")
# for _ in concurrent.futures.ProcessPoolExecutor(max_workers=10).map(worker, range(NUM_TRIALS)):
#     pass
for trial in range(6,8):
    print('===============' + str(trial+1) +'===============' )
    worker(trial)
