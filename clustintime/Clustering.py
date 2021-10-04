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
import Processing as proc
import Visualization as vis
from nilearn.input_data import NiftiMasker
from sklearn.cluster import KMeans


def generate_maps(labels, directory, data_file, mask_file, prefix):
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
    masker = NiftiMasker(mask_file, standardize="zscore")
    data = masker.fit(data_file)
    unique, counts = np.unique(
        labels, return_counts=True
    )  # find the labels and how many in each of them
    unique = unique[counts > 1]
    counts = counts[counts > 1]

    for map_idx in range(len(unique)):
        mean_img = np.mean(data[labels == map_idx], axis=0)
        mean_img_3d = masker.inverse_transform(
            mean_img
        )  # Transform the averaged image into a 3D image
        nib.save(
            mean_img_3d, os.path.join(directory, f"{directory}/{prefix}_cluster_{map_idx}.nii.gz")
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


def K_Means(corr_map, indexes, nscans, n_clusters, task=[], TR=0.5, saving_dir=".", prefix=""):
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

    Returns
    -------
    final_labels : numpy array
        Assigned clusters to each time-point.

    """
    print(" ")
    print("Applying K-Means ...")
    print(" ")

    print(" ")

    KM = KMeans(n_clusters=n_clusters)
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
    vis.plot_heatmap(final_labels, "Assigned clusters for K_Means", task, TR)
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

    corr_map = pd.DataFrame(corr_map).loc[indexes, indexes]
    corr_map[corr_map < np.percentile(corr_map, thr)] = 0
    corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation
    corr_smooth_binary = np.matrix(corr_smooth_binary)

    G = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
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
    vis.plot_heatmap(final_labels, "Labels for the infoMap algorithm", task, TR)
    vis.show_table(final_labels, saving_dir, prefix)
    return corr_map, final_labels
