import pickle
import os
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


N_TIMESTEP_IN_REPLICA = 201
TEMPERATURE = 300

pickle_dir = "./pickle-cache/binding-free-energy_1d"
if not os.path.exists(pickle_dir):
    os.mkdir(pickle_dir)

pacs_com_distances = []
pacs_com_vectors = []
for trial in range(20):
    pickle_com_distances = pickle.load(open(f"./pickle-cache/com-features/pacs_com_distances_{trial}.pickle", "rb"))
    pacs_com_distances.append(pickle_com_distances.reshape((-1, N_TIMESTEP_IN_REPLICA, 1)))
    pickle_com_vectors = pickle.load(open(f"./pickle-cache/com-features/pacs_com_vectors_{trial}.pickle", "rb"))
    pacs_com_vectors.append(pickle_com_vectors.reshape((-1, N_TIMESTEP_IN_REPLICA, 3)))

trajs = pacs_com_distances

from sklearn.cluster import KMeans
import matplotlib


def cluster_trajs(trajs, n_clusters):
    clustering_model = KMeans(
        n_clusters=n_clusters,
        init="k-means++",
        n_init=10,
        max_iter=int(1e10),
        tol=1e-10,
        verbose=0,
        random_state=1,
        copy_x=True,
        algorithm="lloyd"
    )
    trajs_concat = np.concatenate(trajs)
    clustered_trajs_concat = clustering_model.fit_predict(trajs_concat)
    return clustered_trajs_concat.reshape(len(trajs), N_TIMESTEP_IN_REPLICA), clustering_model.cluster_centers_


def plot_clusters(ax, trajs, clustered_trajs, cluster_centers):
    trajs_concat = np.concatenate(trajs)
    clustered_trajs_concat = np.concatenate(clustered_trajs)
    ax.hist(trajs_concat, bins=n_clusters*5)
    for cluster_center in cluster_centers:
        ax.axvline(x=cluster_center, color="black", linewidth=0.2)

from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov.tools import estimation
from deeptime.util.validation import implied_timescales
import deeptime.plots


def construct_msm(clustered_trajs, lagtimes):
    msm_models = {}
    stationary_distributions = {}
    for lagtime in lagtimes:
        count_estimator = TransitionCountEstimator(
            lagtime=lagtime,
            count_mode="sliding",
            n_states=None,
            sparse=False,
        )
        count_estimator.fit(clustered_trajs)
        count_model = count_estimator.fetch_model()
        msm_constructor = MaximumLikelihoodMSM(
            reversible=True,
            stationary_distribution_constraint=None,
            sparse=False,
            allow_disconnected=False,
            maxiter=int(1e10),
            maxerr=1e-10,
            connectivity_threshold=0,
            transition_matrix_tolerance=1e-10,
            lagtime=None,
            use_lcc=False,
        ) 
        try:
            msm_constructor.fit(count_model)
        except:
            continue
        msm_models[lagtime] = msm_constructor.fetch_model()
        largest_connected_set = estimation.largest_connected_set(count_model.count_matrix, directed=True)
        stationary_distributions[lagtime] = np.zeros(np.max(clustered_trajs)+1)
        stationary_distributions[lagtime][largest_connected_set] = msm_models[lagtime].stationary_distribution
    return msm_models, stationary_distributions


def plot_implied_timescales(ax, msm_models):
    deeptime.plots.plot_implied_timescales(
        implied_timescales(list(msm_models.values())),
        n_its=10,
        ax=ax,
        process=None,
        show_mle=True,
        show_sample_mean=True,
        show_sample_confidence=True,
        show_cutoff=True,
        sample_confidence=0.95,
        colors=None,
    )
    ax.set_xlim(0, np.max(list(msm_models.keys())))
    ax.set_yscale("log")


n_clusters_candidates = [60, 40, 20]

if False:
    clustered_trajs = [None for _ in range(len(trajs))]
    cluster_centers = [None for _ in range(len(trajs))]
    for trial in range(len(trajs)):
        clustered_trajs[trial] = {}
        cluster_centers[trial] = {}
        for i_clusters, n_clusters in enumerate(n_clusters_candidates):
            print(trial, n_clusters)
            clustered_trajs[trial][n_clusters], cluster_centers[trial][n_clusters] = cluster_trajs(trajs[trial], n_clusters)
    pickle.dump(clustered_trajs, open(f"{pickle_dir}/clustered_trajs.pickle", "wb"))
    pickle.dump(cluster_centers, open(f"{pickle_dir}/cluster_centers.pickle", "wb"))

clustered_trajs = pickle.load(open(f"{pickle_dir}/clustered_trajs.pickle", "rb"))
cluster_centers = pickle.load(open(f"{pickle_dir}/cluster_centers.pickle", "rb"))


fig_cluster, ax_cluster = plt.subplots(len(trajs), len(n_clusters_candidates))
fig_cluster.set_size_inches(4*len(n_clusters_candidates), 4*len(trajs))

for trial in range(len(trajs)):
    for i_clusters, n_clusters in enumerate(n_clusters_candidates):
        plot_clusters(ax_cluster[trial][i_clusters], trajs[trial], clustered_trajs[trial][n_clusters], cluster_centers[trial][n_clusters])
        ax_cluster[trial][i_clusters].set_title(f"{trial}-{n_clusters}")

fig_cluster.tight_layout()
fig_cluster.savefig("figures/clusters2_all.png", dpi=300)
plt.close(fig_cluster)


if False:
    msm_models = [None for _ in range(len(trajs))]
    stationary_distributions = [None for _ in range(len(trajs))]
    for trial in range(len(trajs)):
        msm_models[trial] = {}
        stationary_distributions[trial] = {}
        for i_clusters, n_clusters in enumerate(n_clusters_candidates):
            print(trial, n_clusters)
            msm_models[trial][n_clusters], stationary_distributions[trial][n_clusters] = construct_msm(clustered_trajs[trial][n_clusters], range(2, int(N_TIMESTEP_IN_REPLICA/2), 2))
    pickle.dump(msm_models, open(f"{pickle_dir}/msm_models.pickle", "wb"))
    pickle.dump(stationary_distributions, open(f"{pickle_dir}/stationary_distributions.pickle", "wb"))

msm_models = pickle.load(open(f"{pickle_dir}/msm_models.pickle", "rb"))
stationary_distributions = pickle.load(open(f"{pickle_dir}/stationary_distributions.pickle", "rb"))


fig_its, ax_its = plt.subplots(len(trajs), len(n_clusters_candidates))
fig_its.set_size_inches(4*len(n_clusters_candidates), 4*len(trajs))

for trial in range(len(trajs)):
    for i_clusters, n_clusters in enumerate(n_clusters_candidates):
        plot_implied_timescales(ax_its[trial][i_clusters], msm_models[trial][n_clusters])
        ax_its[trial][i_clusters].set_title(f"{trial}-{n_clusters}")

fig_its.tight_layout()
fig_its.savefig("figures/its_all.png", dpi=300)
plt.close(fig_its)

import scipy.constants
from scipy.spatial import ConvexHull


def calculate_pmf(stationary_distributions):
    pmf = -np.log(stationary_distributions/max(stationary_distributions))
    pmf *= scipy.constants.Boltzmann
    pmf *= TEMPERATURE
    pmf *= scipy.constants.Avogadro
    pmf /= scipy.constants.calorie
    pmf /= 1000
    return pmf


def calculate_pmf_plot(cluster_centers, pmf, unbound_threshold):
    com_distances_cluster_centers = np.linalg.norm(cluster_centers, axis=1)
    indices_order = np.argsort(com_distances_cluster_centers)
    ordered_pmf = pmf[indices_order]
    ordered_cluster_centers = com_distances_cluster_centers[indices_order]

    pmf_baseline = []
    for i_cluster_center, (cluster_center, _pmf) in enumerate(zip(ordered_cluster_centers, ordered_pmf)):
        if unbound_threshold < cluster_center and np.isfinite(_pmf):
            pmf_baseline.append(_pmf)
    ordered_pmf -= np.mean(pmf_baseline)

    return ordered_cluster_centers, ordered_pmf


def plot_pmf(ax, ordered_cluster_centers, ordered_pmf, xmin, xmax, ymin, ymax, bound_threshold=None, unbound_threshold=None):
    ax.hlines(0, xmin=xmin, xmax=xmax, color="black")
    if bound_threshold != None and unbound_threshold != None:
        ax.vlines(bound_threshold,   ymin=ymin, ymax=ymax, color="black", linestyles="--")
        ax.vlines(unbound_threshold, ymin=ymin, ymax=ymax, color="black", linestyles="--")

    ax.plot(ordered_cluster_centers, ordered_pmf)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def calculate_g_pmf(cluster_centers, stationary_distributions, bound_threshold, unbound_threshold):
    bound_probability = 0
    unbound_probability = 0

    for cluster_center, stationary_distribution in zip(cluster_centers, stationary_distributions):
        com_distances_cluster_centers = np.linalg.norm(cluster_center)
        
        if com_distances_cluster_centers < bound_threshold:
            bound_probability += stationary_distribution
        if  unbound_threshold < com_distances_cluster_centers:
            unbound_probability += stationary_distribution

    g_pmf = -np.log(bound_probability/unbound_probability)
    g_pmf *= scipy.constants.Boltzmann
    g_pmf *= TEMPERATURE
    g_pmf *= scipy.constants.Avogadro
    g_pmf /= scipy.constants.calorie
    g_pmf /= 1000

    return g_pmf


def calculate_volume_correction(com_vectors, unbound_threshold):
    unbound_com_vectors = com_vectors[unbound_threshold < np.linalg.norm(com_vectors, axis=2)]
    hull = ConvexHull(unbound_com_vectors)

    volume_correction = -np.log(hull.volume/1.661)
    volume_correction *= scipy.constants.Boltzmann
    volume_correction *= TEMPERATURE
    volume_correction *= scipy.constants.Avogadro
    volume_correction /= scipy.constants.calorie
    volume_correction /= 1000

    return volume_correction


pmfs = [None for _ in range(len(trajs))]

for trial in range(len(trajs)):
    print(trial)
    pmfs[trial] = {}
    for i_clusters, n_clusters in enumerate(n_clusters_candidates):
        pmfs[trial][n_clusters] = {}
        for lagtime in stationary_distributions[trial][n_clusters]:
            try:
                pmfs[trial][n_clusters][lagtime] = calculate_pmf(stationary_distributions[trial][n_clusters][lagtime])
            except:
                print("failed")
                continue

bound_thresholds = [1 for _ in range(len(trajs))]
unbound_thresholds = [2.5 for _ in range(len(trajs))]

n_clusters = [20 for _ in range(len(trajs))]
lagtimes   = [20 for _ in range(len(trajs))]

#fig_trial, ax_trial = plt.subplots(1, len(trajs))
fig_trial, ax_trial = plt.subplots(4,5)
fig_trial.set_size_inches(1.5*len(trajs),18)

fig_all, ax_all = plt.subplots(1,1)
fig_all.set_size_inches(6, 4)

ax_all.hlines(0, xmin=0, xmax=12, color="black")

for trial in range(len(trajs)):
    ordered_cluster_centers, ordered_pmf = calculate_pmf_plot(cluster_centers[trial][n_clusters[trial]], pmfs[trial][n_clusters[trial]][lagtimes[trial]], unbound_thresholds[trial])

    plot_pmf(ax_all,          ordered_cluster_centers, ordered_pmf, 0, 5, -10, 2)
    # plot_pmf(ax_trial[trial], ordered_cluster_centers, ordered_pmf, 0, 6, -10, 2, bound_threshold=bound_thresholds[trial], unbound_threshold=unbound_thresholds[trial])
    # ax_trial[trial].set_title(trial)
    row = trial // 5
    col = trial % 5

    plot_pmf(ax_trial[row, col], ordered_cluster_centers, ordered_pmf, 0, 5, -10, 2,
            bound_threshold=bound_thresholds[trial], unbound_threshold=unbound_thresholds[trial])
    ax_trial[row, col].set_title("trial " + str(trial+1))

fig_trial.tight_layout(pad=2.0, w_pad=1.5, h_pad=2.0)
fig_all.tight_layout()

fig_trial.savefig("./figures/pmf_each_trial.png", dpi=400, bbox_inches="tight")
fig_all.savefig("./figures/pmf_all_trials.png", dpi=400, bbox_inches="tight")

g_pmfs = np.zeros(len(trajs))
volume_corrections = np.zeros(len(trajs))
for trial in range(len(trajs)):
    g_pmfs[trial] = calculate_g_pmf(cluster_centers[trial][n_clusters[trial]], stationary_distributions[trial][n_clusters[trial]][lagtimes[trial]], bound_thresholds[trial], unbound_thresholds[trial])
    volume_corrections[trial] = calculate_volume_correction(pacs_com_vectors[trial], unbound_thresholds[trial])

g_pmfs = np.zeros(len(trajs))
volume_corrections = np.zeros(len(trajs))
for trial in range(len(trajs)):
    g_pmfs[trial] = calculate_g_pmf(cluster_centers[trial][n_clusters[trial]], stationary_distributions[trial][n_clusters[trial]][lagtimes[trial]], bound_thresholds[trial], unbound_thresholds[trial])
    volume_corrections[trial] = calculate_volume_correction(pacs_com_vectors[trial], unbound_thresholds[trial])

g_stds = g_pmfs + volume_corrections

for trial in range(len(trajs)):
    print(trial, g_pmfs[trial], volume_corrections[trial], g_stds[trial])
print()
print("mean", g_pmfs.mean(), volume_corrections.mean(), g_stds.mean())
print("std ", g_pmfs.std(),  volume_corrections.std(),  g_stds.std())
print(g_stds)

g_stds_ordered = []
for trial, g_std in enumerate(g_stds):
    g_stds_ordered.append((trial,  g_std))
g_stds_ordered.sort(key=lambda x: x[1], reverse=True)

print(g_stds_ordered)
print(f"[{','.join([str(val[0]) for val in g_stds_ordered])}]")

