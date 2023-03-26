#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering library for clustintime
"""
import os
import random

import infomap
import networkx as nx  # creation, manipulation and study of the structure, dynamics and functions of complex networks
import nibabel as nib
import numpy as np
import pandas as pd
from community import community_louvain
from networkx.algorithms import community
from sklearn.cluster import AgglomerativeClustering, KMeans

import clustintime.processing as proc


def generate_maps(labels, directory, data, masker, prefix):
    """
    Creates the activity maps of each cluster

    Parameters
    ----------
    labels : numpy array
        Array of assigned clusters.
    directory : str or path
        Fullpath to the saving directory.
    data_file : str or path
        Fullpath to the data file.
    mask_file : str or path
        Fullpath to the mask file.
    prefix : str
        Prefix for the saved files.

    Returns
    -------
    None.

    """

    unique, counts = np.unique(labels, return_counts=True)  # find the labels and how many in each of them
    unique = unique[counts > 1]
    counts = counts[counts > 1]

    for map_idx in range(len(unique)):
        mean_img = np.mean(data[labels == map_idx + 1], axis=0)
        if mean_img.min() / mean_img.max() < 0.9:
            mean_img = mean_img / mean_img.max()
        mean_img_3d = masker.inverse_transform(mean_img)  # Transform the averaged image into a 3D image

        nib.save(
            mean_img_3d,
            os.path.join(directory, f"{directory}/{prefix}_cluster_{map_idx+1}.nii.gz"),
        )


def find_communities(graph):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id and return number of communities found.
    """
    infomap_network = infomap.Infomap("--two-level")

    print("Building Infomap network from a NetworkX graph...")
    for edge in graph.edges():
        infomap_network.network.addLink(*edge)

    print("Find communities with Infomap...")
    infomap_network.run()

    print(f"Found {infomap_network.numTopModules()} modules with codelength: {infomap_network.codelength}")

    communities = {}
    for node in infomap_network.iterLeafNodes():
        communities[node.physicalId] = node.moduleIndex()

    nx.set_node_attributes(graph, values=communities, name="community")
    return communities


def consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr):
    """
    This algorithm samples the data and clusters it with a defined algorithm.
    With the results of each cluster, it creates a consensus matrix.
    The consensus matrix is then clustered a hundred times.
    If the results are the same in every run, that will be the returned labels.
    If the results are not unanimous, the algorithm will return to the sampling step.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    algorithm : function
        Algorithm to employ


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    npoints = len(indexes)
    sum_connectivity_matrix = np.zeros([npoints, npoints])
    index_matrix_sum = np.zeros([npoints, npoints])

    are_clusters_stable = False
    while not are_clusters_stable:
        for _ in range(1):
            sampling = np.sort(random.sample(range(npoints), round(npoints * 0.6)))
            index_matrix = pd.DataFrame([0] * npoints)
            index_matrix[0][sampling] = 1
            index_matrix_sum = index_matrix_sum + np.dot(index_matrix, np.transpose(index_matrix))
            data_sampled = corr_map[sampling, :][:, sampling]
            if algorithm in (info_map, greedy_mod, louvain):  # pylint: disable=comparison-with-callable
                _, idx = algorithm(data_sampled, indexes, thr, nscans)
                idx = np.transpose(
                    pd.DataFrame([idx, sampling])
                )  # Create a vector that combines the previous indexes and the labels
            else:
                idx = algorithm(data_sampled, indexes, nscans, n_clusters)
                idx = np.transpose(
                    pd.DataFrame([idx, sampling])
                )  # Create a vector that combines the previous indexes and the labels
            idx = idx.set_index(1)
            idx = idx[np.logical_not(np.isnan(idx[0]))]
            labels = np.array([0] * npoints)
            labels[sampling] = idx[0]
            connectivity_matrix = proc.compute_connectivity_matrix(npoints, labels)
        sum_connectivity_matrix = sum_connectivity_matrix + connectivity_matrix
        _consensus = np.divide(sum_connectivity_matrix, index_matrix_sum)
        _consensus[_consensus < proc.find_threshold_bfs(_consensus)] = 0
        final_labels = algorithm(_consensus, indexes, 0, nscans)
        thr = proc.find_threshold_bfs(_consensus)
        _consensus[_consensus <= thr] = 0
        aux = proc.compute_connectivity_matrix(npoints, final_labels[1])

        are_clusters_stable = check_if_clusters_stable(algorithm, _consensus, indexes, thr, npoints, aux)

    return final_labels


def check_if_clusters_stable(algorithm, _consensus, indexes, thr, npoints, aux):
    for _ in range(1):
        labels = algorithm(corr_map=_consensus, indexes=indexes, thr=thr, nscans=npoints)
        connect = proc.compute_connectivity_matrix(npoints, labels[1])

        if not np.array_equal(aux, connect):
            return False

    return True


def k_means(corr_map, indexes, nscans, n_clusters, seed=0):
    """
    K-Means uses a pre-stablished number of centroids and iterations defined by the user.
    The algorithms places the centroids at random locations (real or imaginary,
    that represent the centre of the cluster) and then allocates each data point to the nearest cluster.
    Afterwards, it will optimise the position of those centroids in the number of iterations defined.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    nscans : int
        Number of scans.
    n_clusters : int
        Number of clusters.
    seed : int
        Seed for the algorithm

    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying K-Means ...")
    print(" ")

    print(" ")

    k_m = KMeans(n_clusters=n_clusters, random_state=seed)
    labels = k_m.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        i = int(i)
        final_labels[i] = labels[0][i] + 1
    print("K-Means applied!")

    return final_labels


def agglomerative_clustering(corr_map, indexes, nscans, n_clusters, affinity="euclidean", linkage="ward"):
    """
    Agglomerative Clustering recursively merges the pair of clusters that minimally increases a given linkage distance.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    nscans : int
        Number of scans.
    n_clusters : int
        Number of clusters.
    affinity : str, optional
        Metric used to compute the linkage. Can be `euclidean`, `l1`, `l2`, `manhattan`, `cosine`, or `precomputed`.
        If linkage is `ward`, only `euclidean` is accepted.
        If `precomputed`, a distance matrix (instead of a similarity matrix) is needed as input for the fit method.
        The default is `euclidean`
    linkage : str, optional:
        Linkage criterion to use. The linkage criterion determines which distance to use between sets of observation.
        The algorithm will merge the pairs of cluster that minimize this criterion.
        The options are `ward`, `complete`, `average`, `single`

    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Agglomerative Clustering ...")
    print(" ")

    print(" ")

    agglo_clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity=affinity, linkage=linkage)
    labels = agglo_clustering.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        i = int(i)
        final_labels[i] = labels[0][i] + 1
    print("Agglomerative Clustering applied!")

    return final_labels


def info_map(corr_map, indexes, thr, nscans):
    """
    InfoMap uses information theory to find communities. In particular, it employs the Huffman code to understand the
    flow of information within a graph. This code assigns a prefix to each node, then a prefix to each community.
    When a random walker enters a network, the probability that it transitions between two nodes is given by its
    Markov transition matrix. Nonetheless, once the walker find itself inside a region, it is relatively improbable
    that it transitions onto another.
    InfoMap uses a random walker and applies the aforementioned theories to to find regions and nodes belonging to them.


    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    thr : int
        Percentile threshold for the binarization.
    nscans : int
        Number of scans.

    Returns
    -------
    corr_map : matrix
        Binary correlation map.
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying InfoMap... ")
    print(" ")

    corr_map = proc.thr_index(corr_map, thr)
    corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation

    graph = nx.from_numpy_matrix(np.matrix(corr_smooth_binary))  # Again the binary
    coms = find_communities(graph)  # Clustering

    coms_labels = np.zeros(corr_map.shape[0])

    for key, value in coms.items():
        coms_labels[key] = value

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(sum(nscans))

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Infomap applied")

    return corr_smooth_binary, final_labels


def louvain(corr_map, indexes, thr, nscans):
    """
    Louvain's algorithm maximises modularity and implements an extra step to ensure community properties in the network.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    thr : int
        Percentile threshold for the binarization.
    nscans : int
        Number of scans.

    Returns
    -------
    corr_map : matrix
        Binary correlation map.
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Louvain... ")
    print(" ")

    # compute the best partition

    corr_map = proc.thr_index(corr_map, thr)

    corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation

    graph = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
    partition = community_louvain.best_partition(graph)

    coms_labels = np.zeros(corr_map.shape[0])

    for key, value in partition.items():
        coms_labels[key] = value

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Louvain applied")

    return corr_smooth_binary, final_labels


def greedy_mod(corr_map, indexes, thr, nscans):
    """
    Greedy modularity maximises modularity.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    thr : int
        Percentile threshold for the binarization.
    nscans : int
        Number of scans.

    Returns
    -------
    corr_map : matrix
        Binary correlation map.
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Greedy Modularity... ")
    print(" ")

    # compute the best partition

    corr_map = proc.thr_index(corr_map, thr)

    corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation

    graph = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
    partition = list(community.greedy_modularity_communities(graph))

    coms_labels = np.zeros(corr_map.shape[0])

    for idx, com in enumerate(partition):
        coms_labels[list(com)] = idx + 1

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Greedy Modularity applied")

    return corr_smooth_binary, final_labels
