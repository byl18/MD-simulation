from cProfile import label
import os
import matplotlib.pyplot as plt


cvs = {}
max_trial = 0
max_cycle = 0
for trial in range(1, 100):
    trial_dir = "/gs/bs/tga-Kitao-Lab/yilan/projects/PreQ0/trial01/pacs2/trial" + str(trial).zfill(3)
    if not os.path.exists(trial_dir):
        continue
    max_trial = trial
    cv = []
    for cycle in range(0, 1000):
        cycle_dir = trial_dir + "/cycle" + str(cycle).zfill(3)
        if not os.path.exists(cycle_dir):
            break
        cv_ranked_path = cycle_dir + "/summary/cv_ranked.log"
        if not os.path.exists((cv_ranked_path)):
            break
        max_cycle = max(max_cycle, cycle)
        top_cv_value = float(open(cv_ranked_path).readline().split()[5])
        cv.append(top_cv_value)
    cvs[trial] = cv
    plt.plot(cv, label=f"trial{trial}")
plt.legend()
plt.savefig("cvs.png")


tsv = "cycle," + ",".join([str(val) for val in cvs.keys()]) + "\n"
for cycle in range(max_cycle + 1):
    line = str(cycle) + ","
    for trial, cv in cvs.items():
        if len(cv) <= cycle:
            line += ","
        else:
            line += f"{cv[cycle]},"
    tsv += line + "\n"

open("cvs.csv", "w").write(tsv)
