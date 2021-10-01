#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 09:39:05 2021

@author: ctobias
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt # For graphs
import pandas as pd

import networkx as nx # creation, manipulation and study of the structure, dynamics and functions of complex networks
from IPython.display import HTML, display

from kmapper import KeplerMapper, Cover
from sklearn.decomposition import PCA 
from sklearn.manifold import TSNE
from umap.umap_ import UMAP
from dyneusr.mapper.utils import optimize_dbscan
from dyneusr import DyNeuGraph


def display_table(data):
    html = "<table>"
    for row in data:
        html += "<tr>"
        for field in row:
            html += "<td>%s<td>"%(field)
        html += "</tr>"
    html += "</table>"
    display(HTML(html))


def plot_labels(labels, title, multi_tap = [], single_tap = [], TR = 0.5):
    
    plt.figure(figsize= [32,16])
    plt.xlim([0, len(labels)*TR])
    plt.ylim([0,labels.max()])
    plt.xlabel('Time in seconds', fontsize = 30)
    plt.ylabel('# Cluster', fontsize = 30)
    
    for i in range(int(labels.max())):
        prueba = np.array([0]*len(labels))
        prueba[np.where(labels == i+1)] = labels[np.where(labels == i+1)] 
        plt.fill_between(np.linspace(0,len(labels)*TR,len(labels)),prueba)

        




def show_table(labels):
    Cluster,Count = np.unique(labels, return_counts = True)
    array= [Cluster, Count, Count/len(labels)]
    display_table(np.transpose([['Cluster'], ['Count'], ['Percentage']]))
    display_table(np.transpose(array))
    
    
def plot_two_maps(map_1, map_2,title_1, title_2 ,single_tap = [], multi_tap = [], contrast = 1 ):
    fig = plt.figure(figsize=(32, 9))
    gs = fig.add_gridspec(1, 2, hspace=0)
        
    (ax1, ax2) = gs.subplots()
    ax1.imshow(map_1, aspect='auto', vmin = -contrast, vmax = contrast, cmap='RdBu_r')
    # Vertical lines to delimit the instant in which the event occurs
    limit = map_1.shape[0]
    ax1.vlines(single_tap, ymin=0, ymax=limit, linewidth=.7)
    ax1.vlines(multi_tap, ymin=0, ymax=limit, linewidth=.7)
    ax1.hlines(single_tap, xmin=0, xmax=limit,
                       linewidth=.7)  # Same for horizontal
    ax1.hlines(multi_tap, xmin=0, xmax=limit, linewidth=.7)
    ax1.set_title(title_1)
    ax1.set_xlim([0,limit])
    ax1.set_ylim([limit,0])
            
        # Plot
        
    im = ax2.imshow(map_2, aspect='auto',vmin=-contrast, vmax=contrast, cmap='RdBu_r')
    limit = map_2.shape[0]
    ax2.vlines(single_tap, ymin=0, ymax=limit, linewidth=.7)
    ax2.vlines(multi_tap, ymin=0, ymax=limit, linewidth=.7)
    ax2.hlines(single_tap, xmin=0, xmax=limit,
                       linewidth=.7)  # Same for horizontal
    ax2.hlines(multi_tap, xmin=0, xmax=limit, linewidth=.7)
    ax2.set_title(title_2)
    ax2.set_xlim([0,limit])
    ax2.set_ylim([limit,0])
    fig.colorbar(im)
    
def Dyn(corr_map, labels,output_file = './dyneusr.html'):
    y = pd.get_dummies(labels.astype(str))
        
    mapper = KeplerMapper(verbose=2)
    X = corr_map
    # Configure projection
    pca = PCA(2, random_state=1)
    umap = UMAP(n_components=2, init=pca.fit_transform(X))
    
        # Construct lens and generate the shape graph
    lens = mapper.fit_transform(
           umap.fit_transform(X, y=None),
           projection=[0,1])  
    graph = mapper.map(
           lens, X=corr_map, 
           cover=Cover(30, 0.4),
           clusterer=optimize_dbscan(corr_map, k=3, p=100.0), )



        # Convert to a DyNeuGraph
    dG = DyNeuGraph(G=graph, y=y)
    
        # Define some custom_layouts
    dG.add_custom_layout(lens, name='lens')
    dG.add_custom_layout(nx.spring_layout, name='nx.spring')
    dG.add_custom_layout(nx.kamada_kawai_layout, name='nx.kamada_kawai')
    dG.add_custom_layout(nx.spectral_layout, name='nx.spectral')
    dG.add_custom_layout(nx.circular_layout, name='nx.circular')

        # Configure some projections
    pca = PCA(2, random_state=1)
    tsne = TSNE(2, init='pca', random_state=1)
    umap = UMAP(n_components=2, init=pca.fit_transform(corr_map))
        
        # Add projections as custom_layouts
    dG.add_custom_layout(pca.fit_transform(X), name='PCA')
    dG.add_custom_layout(tsne.fit_transform(X), name='TSNE')
    dG.add_custom_layout(umap.fit_transform(X, y=None), name='UMAP')
        
        # Visualize 
    dG.visualize(output_file)