import mdtraj as md
from pathlib import Path

trial_dir = Path("trial002")
out_xtc = "pacs_traceback_rnafit.xtc"
out_pdb = "pacs_traceback_rnafit.pdb"

RNA_ATOMS = list(range(1062))

top_pdb = (
    "/gs/bs/tga-Kitao-Lab/yilan/projects/PreQ1/"
    "trial01/pacs2/trial002/cycle000/replica001/rmmol_top.pdb"
)



def read_cv_ranked(path: Path):
    rec = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split()
        rec.append((int(parts[1]), int(parts[3]), float(parts[5])))
    return rec


cycles = sorted(trial_dir.glob("cycle*"))
if not cycles:
    raise RuntimeError("No cycle* found")

trace = []  # (cycle_dir, replica, frame)

last_cycle = cycles[-1]
records = read_cv_ranked(last_cycle / "summary" / "cv_ranked.log")

current_rank = 1
current_replica, current_frame, _ = records[0]

trace.append((last_cycle, current_replica, current_frame))

for cycle in reversed(cycles[:-1]):
    records = read_cv_ranked(cycle / "summary" / "cv_ranked.log")
    rank = current_replica  # 1-based
    if rank < 1 or rank > len(records):
        raise IndexError(
            f"{cycle.name}: rank {rank} out of range (1..{len(records)})"
        )
    current_replica, current_frame, _ = records[rank - 1]
    trace.append((cycle, current_replica, current_frame))


trace.reverse()

# -------------------------------------------------------------
# build trajectory (frame-aware)

traj_list = []
time_offset = 0.0
ref_traj = None

for cycle_dir, replica, frame in trace:
    rep_dir = cycle_dir / f"replica{replica:03d}"
    xtc = rep_dir / "prd_rmmol.xtc"
    traj = md.load(str(xtc), top=top_pdb)
    if frame < traj.n_frames - 1:
        traj = traj[:frame]
    else:
        traj = traj[:-1]
    if ref_traj is None:
        ref_traj = traj[0]
    traj.superpose(ref_traj, atom_indices=RNA_ATOMS)
    # continuous time (optional)
    if traj.n_frames > 1:
        dt = traj.time[1] - traj.time[0]
    else:
        dt = 1.0
    traj.time += time_offset
    time_offset = traj.time[-1] + dt
    traj_list.append(traj)
    print(
        f"{cycle_dir.name}: replica {replica}, "
        f"frame {frame}, kept {traj.n_frames} frames"
    )


full = md.join(traj_list)
full.save_xtc(out_xtc)
full[0].save_pdb(out_pdb)

print("\nSaved:")
print(out_xtc)
print(out_pdb)
print("Total frames:", full.n_frames)
