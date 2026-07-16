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

if True:
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


if True:
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
