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

import clustintime.Visualization as vis
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import MeanShift
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
        mean_img = mean_img/mean_img.max()
        mean_img_3d = masker.inverse_transform(
            mean_img
        )  # Transform the averaged image into a 3D image
        
        nib.save(
            mean_img_3d, os.path.join(directory, f"{directory}/{prefix}_cluster_{map_idx+1}.nii.gz")
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


def K_Means(corr_map, indexes, nscans, n_clusters, task=[], TR=0.5, saving_dir=".", prefix="", seed = 0):
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
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results
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

    KM = KMeans(n_clusters=n_clusters, random_state = seed)
    labels = KM.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("K-Means applied!")
    vis.plot_labels(final_labels, "Assigned clusters for K_Means", task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for K_Means", task = task, TR = TR,  saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels


def Agglomerative_Clustering(corr_map, indexes, nscans, n_clusters, affinity = 'euclidean', linkage = 'ward',task=[], TR=0.5, saving_dir=".", prefix=""):
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
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Agglomerative Clustering ...")
    print(" ")

    print(" ")

    AG = AgglomerativeClustering(n_clusters = n_clusters, affinity = affinity, linkage = linkage)
    labels = AG.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("Agglomerative Clustering applied!")
    vis.plot_labels(final_labels, "Assigned clusters for Agglomerative Clustering", task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for Agglomerative Clustering", task = task, TR = TR, saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels

def dbscan(corr_map, indexes, nscans, eps, metric = 'euclidean', algorithm = 'auto',task=[], TR=0.5, saving_dir=".", prefix=""):
    """
    Density-Based Spatial Clustering of Applications with Noise. Finds core samples of high density and expands clusters from them.

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
    eps : float
        The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
    metric : str, optional
        Metric used to compute the linkage. Can be `euclidean`, `l1`, `l2`, `manhattan`, `cosine`, or `precomputed`. 
        If linkage is `ward`, only `euclidean` is accepted. 
        If `precomputed`, a distance matrix (instead of a similarity matrix) is needed as input for the fit method.
        The default is `euclidean`
    algorithm : str, optional:
        The algorithm to be used by the NearestNeighbors module to compute pointwise distances and find nearest neighbors.
        The options are `auto`, `ball_tree`, `kd_tree`, `brute`. The default is `auto`
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying DBScan ...")
    print(" ")

    print(" ")

    DBScan = DBSCAN(eps = eps, metric = metric, algorithm = algorithm).fit(corr_map)
    labels = DBScan.labels_

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("DBScan applied!")
    vis.plot_labels(final_labels, "Assigned clusters for DBScan", task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for DBScan", task = task, TR = TR, saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels

def optics(corr_map, indexes, nscans, algorithm = 'auto',task=[], TR=0.5, saving_dir=".", prefix=""):
    """
    OPTICS (Ordering Points To Identify the Clustering Structure), closely related to DBSCAN, finds core sample of high density and expands clusters from them. 
    Unlike DBSCAN, keeps cluster hierarchy for a variable neighborhood radius. Better suited for usage on large datasets than the current sklearn implementation of DBSCAN.
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

    algorithm : str, optional:
        The algorithm to be used by the NearestNeighbors module to compute pointwise distances and find nearest neighbors.
        The options are `auto`, `ball_tree`, `kd_tree`, `brute`. The default is `auto`
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying OPTICS ...")
    print(" ")

    print(" ")

    optics = OPTICS(algorithm = algorithm).fit(corr_map)
    labels = optics.labels_

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label

    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("DBScan applied!")
    vis.plot_labels(final_labels, "Assigned clusters for DBScan", task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for DBScan", task = task, TR = TR, saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels

def Info_Map(corr_map, indexes, thr, nscans, task=[], TR=0.5, saving_dir=".", prefix=""):
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
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results
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

    
    corr_map[corr_map < np.percentile(corr_map, thr)] = 0
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

    vis.plot_labels(final_labels, "Labels for the infoMap algorithm", task, TR)
    vis.plot_heatmap(final_labels, "Labels for the infoMap algorithm", task = task, TR = TR,  saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return corr_smooth_binary, final_labels

def Louvain(corr_map, indexes, thr, nscans, task=[], TR=0.5, saving_dir=".", prefix=""):
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
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results
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

    corr_map[corr_map < np.percentile(corr_map, thr)] = 0
    
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

    vis.plot_labels(final_labels, "Labels for the louvain algorithm", task, TR)
    vis.plot_heatmap(final_labels, "Labels for the louvain algorithm", task = task, TR = TR,  saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return corr_smooth_binary, final_labels

def Greedy_Mod(corr_map, indexes, thr, nscans, task=[], TR=0.5, saving_dir=".", prefix=""):
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
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results
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

    corr_map[corr_map < np.percentile(corr_map, thr)] = 0
    
    corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation

    G = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
    partition = list(community.greedy_modularity_communities(G))

    coms_labels = np.zeros(corr_map.shape[0])

    for ii in partition:
        coms_labels[ii] = partition[ii]

    labels = np.transpose(pd.DataFrame([coms_labels, indexes]))
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    for i in labels.index:
        final_labels[int(i)] = labels[0][int(i)] + 1

    print("Greedy Modularity applied")

    vis.plot_labels(final_labels, "Labels for the Greedy Modularity algorithm", task, TR)
    vis.plot_heatmap(final_labels, "Labels for the Greedy Modularity algorithm", task, TR,  saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return corr_smooth_binary, final_labels

def Affinity_Propagation(corr_map, indexes, nscans,damping=0.5,seed = 0, task=[], TR=0.5, saving_dir=".", prefix=""):
    """
    AffinityPropagation creates clusters by sending messages between pairs of samples until convergence. 
    A dataset is then described using a small number of exemplars, which are identified as those most representative of other samples. 
    The messages sent between pairs represent the suitability for one sample to be the exemplar of the other, which is updated in response to the values from other pairs. 
    This updating happens iteratively until convergence, at which point the final exemplars are chosen, and hence the final clustering is given.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    nscans : int
        Number of scans.
    seed : int, optional
        Random State for the algorithm. The default is 0
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Affinity Propagation ...")
    print(" ")

    print(" ")

    Prop = AffinityPropagation(damping = damping,affinity = 'precomputed')
    labels = Prop.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label


    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("Affinity Propagation applied!")
    vis.plot_labels(final_labels, "Assigned clusters for Affinity Propagation" , task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for Affinity Propagation ", task = task, TR = TR, saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels

def Mean_Shift(corr_map, indexes, nscans, task=[], TR=0.5, saving_dir=".", prefix=""):
    """
    Mean shift clustering aims to discover “blobs” in a smooth density of samples. 
    It is a centroid-based algorithm, which works by updating candidates for centroids to be the mean of the points within a given region. 
    These candidates are then filtered in a post-processing stage to eliminate near-duplicates to form the final set of centroids.
    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    indexes : numpy array
        Indexes of the relevant time-points.
    nscans : int
        Number of scans.
    seed : int, optional
        Random State for the algorithm. The default is 0
    task : dictionary or list, optional
        Structure containing the timings of the task. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.
    saving_dir : str or path
        Fullpath of saving directory of the results
    prefix : str
        Desired name for the results


    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying Mean Shift ...")
    print(" ")

    print(" ")

    MS = MeanShift()
    labels = MS.fit_predict(corr_map)

    labels = np.transpose(
        pd.DataFrame([labels, indexes])
    )  # Create a vector that combines the previous indexes and the labels
    labels = labels.set_index(1)

    final_labels = np.zeros(nscans)

    # assign to each timepoint their label


    for i in labels.index:
        final_labels[i] = labels[0][i] + 1
    print("Mean Shift applied!")
    vis.plot_labels(final_labels, "Assigned clusters for Mean Shift" , task, TR)
    vis.plot_heatmap(final_labels, "Assigned clusters for Mean Shift", task = task, TR = TR, saving_dir = saving_dir, prefix = prefix)
    vis.show_table(final_labels, saving_dir, prefix)
    return final_labels