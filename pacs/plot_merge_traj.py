import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re

# 设置中文字体（可选）
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def extract_cv_from_log(log_file):
    """从 cv_ranked.log 提取第一行的 cv 值"""
    try:
        with open(log_file, 'r') as f:
            first_line = f.readline().strip()
            # 匹配格式: replica 9 frame 3 cv 0.27962905480777533
            match = re.search(r'cv\s+([0-9.e-]+)', first_line)
            if match:
                return float(match.group(1))
            else:
                return None
    except Exception as e:
        print(f"Error reading {log_file}: {e}")
        return None

# 存储所有 trial 的数据
all_data = {}
max_cycles = 0

# 循环 20 个 trial
for trial_num in range(1, 21):
    trial_str = f"{trial_num:03d}"
    trial_dir = Path(f"trial{trial_str}")
    
    if not trial_dir.exists():
        print(f"Warning: {trial_dir} not found")
        continue
    
    # 找到所有 cycle 目录
    cycles = sorted(trial_dir.glob("cycle*"))
    cv_values = []
    
    for cycle in cycles:
        log_file = cycle / "summary" / "cv_ranked.log"
        if log_file.exists():
            cv = extract_cv_from_log(log_file)
            if cv is not None:
                cv_values.append(cv)
            else:
                print(f"Warning: Cannot extract cv from {log_file}")
                cv_values.append(np.nan)
        else:
            print(f"Warning: {log_file} not found")
            cv_values.append(np.nan)
    
    if cv_values:
        all_data[trial_num] = cv_values
        max_cycles = max(max_cycles, len(cv_values))
        print(f"Trial {trial_str}: {len(cv_values)} cycles")

# 绘图
fig, ax = plt.subplots(figsize=(12, 8))

# 使用颜色映射
colors = plt.cm.tab20(np.linspace(0, 1, len(all_data)))

for i, (trial_num, cv_values) in enumerate(sorted(all_data.items())):
    cycles_x = range(1, len(cv_values) + 1)
    # 只画到该 trial 的最大 cycle，后面留空
    ax.plot(cycles_x, cv_values, marker='o', markersize=3, 
            linewidth=1.5, label=f"Trial {trial_num:03d}", 
            color=colors[i])

ax.set_xlabel("Cycle", fontsize=14)
ax.set_ylabel("CV Value", fontsize=14)
ax.set_title("CV Value vs Cycle for 20 Trials", fontsize=16)
ax.set_xlim(0.5, max_cycles + 0.5)
ax.grid(True, alpha=0.3)
ax.legend(loc='best', ncol=2, fontsize=8)

plt.tight_layout()
plt.savefig("cv_curves_all_trials.png", dpi=300)
plt.savefig("cv_curves_all_trials.pdf")

print(f"\nTotal trials plotted: {len(all_data)}")
print(f"Max cycles: {max_cycles}")