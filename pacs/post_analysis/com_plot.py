import numpy as np
import pickle
import matplotlib.pyplot as plt
from pathlib import Path

N_TRIALS = 20
WINDOW_SIZE = 200 * 64
all_data = {}

for trial in range(N_TRIALS):
    pickle_path = f"./pickle-cache/com-features/pacs_com_distances_{trial}.pickle"
    try:
        with open(pickle_path, 'rb') as f:
            com_distances = pickle.load(f)
        # 第一段截取200帧，后面每个窗口64*200
        cv_values = []
        # 第一段: 前200帧
        if len(com_distances) >= 200:
            cv_values.append(np.max(com_distances[:200]))
        # 后面每隔WINDOW_SIZE取一个窗口中的最大值
        for i in range(200, len(com_distances), WINDOW_SIZE):
            window_end = min(i + WINDOW_SIZE, len(com_distances))
            cv_values.append(np.max(com_distances[i:window_end]))
        # 如果cv_value大于3后，出现后一段cv_value小于前面的，就不要后面的了
        truncated = False
        for j in range(1, len(cv_values)):
            if cv_values[j-1] > 3 and cv_values[j] < cv_values[j-1]:
                cv_values = cv_values[:j]
                truncated = True
                break
        all_data[trial] = cv_values
        print(f"Trial {trial:03d}: Loaded {len(com_distances)} points, extracted {len(cv_values)} CV values (truncated: {truncated})")
    except FileNotFoundError:
        print(f"Warning: Trial {trial:03d} pickle file not found at {pickle_path}")
    except Exception as e:
        print(f"Error loading Trial {trial:03d}: {e}")

max_cycles = max([len(cv_values) for cv_values in all_data.values()]) if all_data else 0

fig, ax = plt.subplots(figsize=(12, 8))
colors = plt.cm.tab20(np.linspace(0, 1, 20))

for i, (trial_num, cv_values) in enumerate(sorted(all_data.items())):
    cycles_x = range(1, len(cv_values) + 1)
    ax.plot(cycles_x, cv_values, marker='o', markersize=3,
            linewidth=1.5, label=f"Trial {trial_num+1:03d}",
            color=colors[i])

ax.set_xlabel("Cycle", fontsize=14)
ax.set_ylabel("CV Value (COM Distance)", fontsize=14)
ax.set_title("CV Value vs Cycle for 20 Trials", fontsize=16)
ax.set_xlim(0.5, max_cycles + 0.5)
ax.grid(True, alpha=0.3)
ax.legend(loc='best', ncol=2, fontsize=8)

plt.tight_layout()
plt.savefig("figures/cv_curves_all_trials.png", dpi=300)
