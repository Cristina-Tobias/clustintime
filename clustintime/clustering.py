#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering library for clustintime
"""
import os

import infomap
import networkx as nx  # creation, manipulation and study of the structure, dynamics and functions of complex networks
import nibabel as nib
import numpy as np
import pandas as pd
import random

import clustintime.processing as proc
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from networkx.algorithms import community
from community import community_louvain


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

    unique, counts = np.unique(
        labels, return_counts=True
    )  # find the labels and how many in each of them
    unique = unique[counts > 1]
    counts = counts[counts > 1]

    for map_idx in range(len(unique)):
        mean_img = np.mean(data[labels == map_idx + 1], axis=0)
        if mean_img.min() / mean_img.max() < 0.9:
            mean_img = mean_img / mean_img.max()
        mean_img_3d = masker.inverse_transform(
            mean_img
        )  # Transform the averaged image into a 3D image

        nib.save(
            mean_img_3d,
            os.path.join(directory, f"{directory}/{prefix}_cluster_{map_idx+1}.nii.gz"),
        )


def findCommunities(G):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id and return number of communities found.
    """
    infomapX = infomap.Infomap("--two-level")

    print("Building Infomap network from a NetworkX graph...")
    for e in G.edges():
        infomapX.network.addLink(*e)

    print("Find communities with Infomap...")
    infomapX.run()

    print(
        "Found {} modules with codelength: {}".format(
            infomapX.numTopModules(), infomapX.codelength
        )
    )

    communities = {}
    for node in infomapX.iterLeafNodes():
        communities[node.physicalId] = node.moduleIndex()

    nx.set_node_attributes(G, values=communities, name="community")
    return communities


def consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr):
    """
    This algorithm samples the data and clusters it with a defined algorithm. With the results of each cluster, it creates a consensus matrix.
    The consensus matrix is then clustered a hundred times. If the results are the same in every run, that will be the returned labels.
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
    M_sum = np.zeros([npoints, npoints])
    I_sum = np.zeros([npoints, npoints])
    while 1:
        for i in range(100):
            sampling = np.sort(random.sample(range(npoints), round(npoints * 0.6)))
            I = pd.DataFrame([0] * npoints)
            I[0][sampling] = 1
            I_sum = I_sum + np.dot(I, np.transpose(I))
            data_sampled = corr_map[sampling, :][:, sampling]
            if algorithm == Info_Map or algorithm == Greedy_Mod or algorithm == Louvain:
                corr_mat, idx = algorithm(data_sampled, indexes, thr, nscans)
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
            M = proc.compute_connectivity_matrix(npoints, labels)
        M_sum = M_sum + M
        Consensus = np.divide(M_sum, I_sum)
        Consensus[Consensus < proc.find_threshold_bfs(Consensus)] = 0
        final_labels = algorithm(Consensus, indexes, 0, nscans)
        thr = proc.find_threshold_bfs(Consensus)
        Consensus[Consensus <= thr] = 0
        aux = proc.compute_connectivity_matrix(npoints, final_labels[1])
        boolean = True
        for i in range(100):
            labels = algorithm(corr_map=Consensus, indexes=indexes, thr=thr, nscans=npoints)
            connect = proc.compute_connectivity_matrix(npoints, labels[1])

            if np.array_equal(aux, connect) == False:
                boolean = False
                break
        if boolean:
            break
    return labels[1]


def K_Means(corr_map, indexes, nscans, n_clusters, seed=0):
    """
    K-Means uses a pre-stablished number of centroids and iterations defined by the user.
    The algorithms places the centroids at random locations (real or imaginary, that represent the centre of the cluster) and then allocates each data point to the nearest cluster.
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

    KM = KMeans(n_clusters=n_clusters, random_state=seed)
    labels = KM.fit_predict(corr_map)

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


def Agglomerative_Clustering(
    corr_map, indexes, nscans, n_clusters, affinity="euclidean", linkage="ward"
):
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
        Linkage criterion to use. The linkage criterion determines which distance to use between sets of observation. The algorithm will merge the pairs of cluster that minimize this criterion.
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

    AG = AgglomerativeClustering(n_clusters=n_clusters, affinity=affinity, linkage=linkage)
    labels = AG.fit_predict(corr_map)

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


def Info_Map(corr_map, indexes, thr, nscans):
    """
    InfoMap uses information theory to find communities. In particular, it employs the Huffman code to understand the flow of information within a graph. This code assigns a prefix to each node, then a prefix to each community.
    When a random walker enters a network, the probability that it transitions between two nodes is given by its Markov transition matrix. Nonetheless, once the walker find itself inside a region, it is relatively improbable that it transitions onto another.
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

    G = nx.from_numpy_matrix(np.matrix(corr_smooth_binary))  # Again the binary
    coms = findCommunities(G)  # Clustering

    coms_labels = np.zeros(corr_map.shape[0])

    for ii in coms:
        coms_labels[ii] = coms[ii]

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Infomap applied")

    return corr_smooth_binary, final_labels


def Louvain(corr_map, indexes, thr, nscans):
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

    G = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
    partition = community_louvain.best_partition(G)

    coms_labels = np.zeros(corr_map.shape[0])

    for ii in partition:
        coms_labels[ii] = partition[ii]

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Louvain applied")

    return corr_smooth_binary, final_labels


def Greedy_Mod(corr_map, indexes, thr, nscans):
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

    G = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
    partition = list(community.greedy_modularity_communities(G))

    coms_labels = np.zeros(corr_map.shape[0])

    for num, ii in enumerate(partition):
        coms_labels[list(ii)] = num + 1

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Greedy Modularity applied")

    return corr_smooth_binary, final_labels
