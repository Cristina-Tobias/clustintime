#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 10:16:51 2021

@author: ctobias
"""
import numpy as np
import Visualization as vis
import infomap
from sklearn.cluster import KMeans
import networkx as nx # creation, manipulation and study of the structure, dynamics and functions of complex networks
import os
import nibabel as nib
from nilearn.masking import apply_mask 
from nilearn.input_data import NiftiMasker

def save_maps(labels, directory, data_file, mask_file):
    data = apply_mask(data_file,mask_file) # apply mask to the fitted signal
    masker = NiftiMasker(mask_file, standardize='zscore')
    unique, counts = np.unique(labels, return_counts = True) # find the labels and how many in each of them
    unique = unique[counts > 1]
    counts = counts[counts > 1]
    
    for map_idx in range(len(unique)):
        mean_img = np.mean(data[labels == map_idx], axis = 0)
        mean_img_3d = masker.inverse_transform(mean_img) #Transform the averaged image into a 3D image
        nib.save(mean_img_3d, os.path.join(directory, f'{directory}/filename_{map_idx}.nii.gz'))
    os.system('module load afni/latest')
    os.system(f'3drefit -space ORIG -view orig {directory}/filename_*')

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
    infomapX.run();

    print("Found {} modules with codelength: {}".format(infomapX.numTopModules(), infomapX.codelength))

    communities = {}
    for node in infomapX.iterLeafNodes():
        communities[node.physicalId] = node.moduleIndex()

    nx.set_node_attributes(G, values=communities, name='community')
    return communities

def K_Means(corr_map, n_clusters, multi_tap = [], single_tap = []):
    print(' ')
    print('K-Means ')
    print(' ')
    print('K-Means uses a pre-stablished number of centroids and iterations defined by the user. The algorithms places the centroids at random locations (real or imaginary, that represent the centre of the cluster) and then allocates each data point to the nearest cluster. Afterwards, it will optimise the position of those centroids in the number of iterations defined.')
    print(' ')
            
    KM= KMeans(n_clusters=n_clusters)
    labels = KM.fit_predict(corr_map)
    vis.plot_labels(labels, 'Assigned clusters for K_Means',multi_tap, single_tap)
    vis.show_table(labels)
    return labels
        
def Info_Map(corr_map, thr, multi_tap = [], single_tap = []):
    print(' ')
    print('InfoMap')
    print(' ')
    print('InfoMap uses information theory to find communities. In particular, it employs the Huffman code to understand the flow of information within a graph. This code assigns a prefix to each node, then a prefix to each community.',
          'When a random walker enters a network, the probability that it transitions between two nodes is given by its Markov transition matrix. Nonetheless, once the walker find itself inside a region, it is relatively improbable that it transitions onto another.',
          'InfoMap uses a random walker and applies the aforementioned theories to to find regions and nodes belonging to them.')
    print(' ')

    corr_smooth_w5 = corr_map.copy()
    thr_3 = np.percentile(corr_smooth_w5, thr)
    corr_smooth_w5[corr_smooth_w5 < thr_3] = 0
    corr_map = corr_smooth_w5 != 0 # Find all the voxels with correlation
    
    G = nx.from_numpy_matrix(corr_map) # Again the binary
    coms = findCommunities(G) # Clustering
    
    coms_labels = np.zeros((corr_map.shape[0]))
            
    for ii in coms:
        coms_labels[ii] = coms[ii]
        
    vis.plot_labels(coms_labels, 'Labels for the infoMap algorithm', multi_tap, single_tap)
    vis.show_table(coms_labels)   
    return corr_map, coms_labels     
