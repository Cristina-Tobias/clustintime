#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustering library for clustintime
"""
import os
import infomap
import networkx as nx  
import nibabel as nib
import numpy as np
import pandas as pd
from community import community_louvain
from networkx.algorithms import community
from sklearn.cluster import AgglomerativeClustering, KMeans
from clustintime.processing import Processing

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
    if len(labels.shape)>1:
        nr_labels = len(labels)
        nr_subjects = labels.shape[1]
        labels = np.reshape(labels, (nr_labels*nr_subjects,), order='F')
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

def create_labels(clustering_results, nscans):
    if type(nscans) is not int:
        labels = np.zeros([np.max(nscans), len(nscans)])
        subject = 0
        for idx in clustering_results.index:
            labels[int(idx) - subject*max(nscans), subject] = clustering_results[0][int(idx)] + 1
            if (int(idx)+1)%max(nscans) == 0:
                subject+=1
    else:
        labels = np.zeros(nscans)
        for idx in clustering_results.index:
            labels[int(idx)] = clustering_results[0][int(idx)] + 1
    return labels

class Clustering:
    def __init__(self, corr_map, indices, nscans):
        self.corr_map = corr_map
        self.indices = indices
        self.nscans = nscans
    def k_means(self, n_clusters, seed=None):
        """
        K-Means uses a pre-stablished number of centroids and iterations defined by the user.
        The algorithms places the centroids at random locations (real or imaginary,
        that represent the centre of the cluster) and then allocates each data point to the nearest cluster.
        Afterwards, it will optimise the position of those centroids in the number of iterations defined.
    
        Parameters
        ----------
        n_clusters : int
            Number of clusters.
        seed : int
            Seed for the algorithm
    
        Returns
        -------
        final_labels : numpy array
            Assigned clusters to each time-point.
    
        """
 
        k_m = KMeans(n_clusters=n_clusters, random_state=seed)
        labels = k_m.fit_predict(self.corr_map)
        labels = np.transpose(
            pd.DataFrame([labels, self.indices])
        )  
        labels = labels.set_index(1)  
        results = create_labels(labels, self.nscans)
    
        return results
    
    
    def agglomerative_clustering(self, n_clusters, affinity="euclidean", linkage="ward"):
    
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
    
        agglo_clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity=affinity, linkage=linkage)
        labels = agglo_clustering.fit_predict(self.corr_map)
        labels = np.transpose(
            pd.DataFrame([labels, self.indices])
        )  
        labels = labels.set_index(1)
    
        results = create_labels(labels, self.nscans)
   
        return results
    
    
    def info_map(self,thr):
        """
        InfoMap uses information theory to find communities. In particular, it employs the Huffman code to understand the
        flow of information within a graph. This code assigns a prefix to each node, then a prefix to each community.
        When a random walker enters a network, the probability that it transitions between two nodes is given by its
        Markov transition matrix. Nonetheless, once the walker find itself inside a region, it is relatively improbable
        that it transitions onto another.
        InfoMap uses a random walker and applies the aforementioned theories to to find regions and nodes belonging to them.
    
    
        Parameters
        ----------
        thr : int
            Percentile threshold for the binarization.

    
        Returns
        -------
        corr_map : matrix
            Binary correlation map.
        final_labels : numpy array
            Assigned clusters to each time-point.
    
        """

        corr_map = Processing(self.corr_map).thr_index(thr)
        corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation
    
        graph = nx.from_numpy_array(np.matrix(corr_smooth_binary))  # Again the binary
        coms = find_communities(graph)  # Clustering
    
        coms_labels = np.zeros(corr_map.shape[0])
    
        for key, value in coms.items():
            coms_labels[key] = value
    
        labels = np.transpose(pd.DataFrame([coms_labels, self.indices]))
        labels = labels.set_index(1)
    
        results = create_labels(labels, self.nscans)


    
        return corr_smooth_binary, results
    
    
    def louvain(self, thr):
    
        """
        Louvain's algorithm maximises modularity and implements an extra step to ensure community properties in the network.
    
        Parameters
        ----------

        thr : int
            Percentile threshold for the binarization.

        Returns
        -------
        corr_map : matrix
            Binary correlation map.
        final_labels : numpy array
            Assigned clusters to each time-point.
    
        """

    
        corr_map = Processing(self.corr_map).thr_index(thr)
    
        corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation
    
        graph = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
        partition = community_louvain.best_partition(graph)
    
        coms_labels = np.zeros(corr_map.shape[0])
    
        for key, value in partition.items():
            coms_labels[key] = value
    
        labels = np.transpose(pd.DataFrame([coms_labels, self.indices]))
        labels = labels.set_index(1)
    
        results = create_labels(labels, self.nscans)
    
        print("Louvain applied")
    
        return corr_smooth_binary, results
    
    
    def greedy_mod(self, thr):
    
        """
        Greedy modularity maximises modularity.
    
        Parameters
        ----------

        thr : int
            Percentile threshold for the binarization.

        Returns
        -------
        corr_map : matrix
            Binary correlation map.
        final_labels : numpy array
            Assigned clusters to each time-point.
    
        """
 
        corr_map = Processing(self.corr_map).thr_index(thr)
    
        corr_smooth_binary = corr_map != 0  # Find all the voxels with correlation
    
        graph = nx.from_numpy_matrix(corr_smooth_binary)  # Again the binary
        partition = list(community.greedy_modularity_communities(graph))
    
        coms_labels = np.zeros(corr_map.shape[0])
    
        for idx, com in enumerate(partition):
            coms_labels[list(com)] = idx + 1
    
        labels = np.transpose(pd.DataFrame([coms_labels, self.indices]))
        labels = labels.set_index(1)
    
        results = create_labels(labels, self.nscans)

    
        return corr_smooth_binary, results
    
        



