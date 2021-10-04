#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Result visualization for clustintime
"""

import matplotlib.pyplot as plt  # For graphs
import networkx as nx  # creation, manipulation and study of the structure, dynamics and functions of complex networks

# Libraries
import numpy as np
import pandas as pd
import seaborn as sns
from dyneusr import DyNeuGraph
from dyneusr.mapper.utils import optimize_dbscan
from IPython.display import HTML, display
from kmapper import Cover, KeplerMapper
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap.umap_ import UMAP


def display_table(data):
    """
    HTML code for showing a table of the results

    Parameters
    ----------
    data : numpy array
    Returns
    -------
    None.

    """
    html = "<table>"
    for row in data:
        html += "<tr>"
        for field in row:
            html += "<td>%s<td>" % (field)
        html += "</tr>"
    html += "</table>"
    display(HTML(html))


def plot_labels(labels, title, task=[], TR=0.5):
    """
    Visualization of all the clusters found along time

    Parameters
    ----------
    labels : numpy array
        Array of assigned clusters.
    title : str
        Title for the plot.
    task : dictionary or list, optional
        Structure containing the times when the task is performed. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.

    Returns
    -------
    None.

    """

    colors = sns.color_palette("bright", len(task))

    plt.figure(figsize=[32, 16])
    plt.xlim([0, len(labels) * TR])
    plt.ylim([0, labels.max()])
    plt.xlabel("Time in seconds", fontsize=30)
    plt.ylabel("# Cluster", fontsize=30)

    for i in range(int(labels.max())):
        selected_labels = np.array([0] * len(labels))
        selected_labels[np.where(labels == i + 1)] = labels[np.where(labels == i + 1)]
        plt.fill_between(np.linspace(0, len(labels) * TR, len(labels)), selected_labels)

    for i in range(len(task)):
        plt.vlines(task[i], 0, labels.max(), linewidth=1.2, colors=colors[i])


def plot_heatmap(labels, title, task=[], TR=0.5):
    """
    Visualization of the clusters separately

    Parameters
    ----------
    labels : numpy array
        Array of assigned clusters.
    title : str
        Title for the plot.
    task : dictionary or list, optional
        Structure containing the times when the task is performed. The default is [].
    TR : float, optional
        TR of the data. The default is 0.5.

    Returns
    -------
    None.

    """
    
    heatmatrix = np.zeros([int(labels.max()), len(labels)])
    rownames = np.zeros([int(labels.max())]).astype(str)
    x = np.linspace(0, len(labels) * TR, len(labels)).astype(int)
    
    for i in range(int(labels.max())):
        selected_labels = np.array([0] * len(labels))
        selected_labels[np.where(labels == i + 1)] = 1
        heatmatrix[i] = selected_labels
        rownames[i] = f'Cluster {i+1}'
    heatmatrix = pd.DataFrame(heatmatrix,columns = x ,index = rownames)
    colors = sns.color_palette("bright", len(task))
    sns.heatmap(heatmatrix,  cmap = 'YlGnBu', xticklabels = 150)
    plt.xlabel('Time in seconds', fontsize = 10)
    for j in range(len(task)):
        plt.vlines(task[j]/TR, 0, labels.max() ,linewidth=1.2, colors=colors[j])
    
    
    
        


def show_table(labels, saving_dir, prefix):
    """
    Table display of the results

    Parameters
    ----------
    labels : numpy array
        Array of assigned clusters.

    Returns
    -------
    None.

    """
    Cluster, Count = np.unique(labels, return_counts=True)
    array = [Cluster, Count, Count / len(labels)]
    display_table(np.transpose([["Cluster"], ["Count"], ["Percentage"]]))
    display_table(np.transpose(array))
    table_result = pd.DataFrame({"Cluster": array[0], "Count": array[1], "Percentage": array[2]})
    table_result.to_csv(f"{saving_dir}/{prefix}_Results.csv")


def plot_two_matrixes(map_1, map_2, title_1, title_2, task=[], contrast=1):
    """
    Graphical comparison between two correlation maps

    Parameters
    ----------
    map_1 : matrix
        Correlation matrix before processing.
    map_2 : matrix
        Correlation matrix after processing.
    title_1 : str
        Title for the first matrix.
    title_2 : str
        Title for the second matrix.
    task : dictionary or list, optional
        Structure containing the times when the task is performed. The default is [].
    contrast : int, optional
        Range of values of the correlation matrixes. The default is 1.

    Returns
    -------
    None.

    """
    fig = plt.figure(figsize=(32, 9))
    gs = fig.add_gridspec(1, 2, hspace=0)

    (ax1, ax2) = gs.subplots()
    ax1.imshow(map_1, aspect="auto", vmin=-contrast, vmax=contrast, cmap="RdBu_r")
    # Vertical lines to delimit the instant in which the event occurs
    limit = map_1.shape[0]
    for i in len(task):
        ax1.vlines(task[i], ymin=0, ymax=limit, linewidth=0.7)
        ax1.hlines(task[i], xmin=0, xmax=limit, linewidth=0.7)  # Same for horizontal
    ax1.set_title(title_1)
    ax1.set_xlim([0, limit])
    ax1.set_ylim([limit, 0])

    # Plot

    im = ax2.imshow(map_2, aspect="auto", vmin=-contrast, vmax=contrast, cmap="RdBu_r")
    limit = map_2.shape[0]
    for i in len(task):
        ax2.vlines(task[i], ymin=0, ymax=limit, linewidth=0.7)
        ax2.hlines(task[i], xmin=0, xmax=limit, linewidth=0.7)  # Same for horizontal
    ax2.set_title(title_2)
    ax2.set_xlim([0, limit])
    ax2.set_ylim([limit, 0])
    fig.colorbar(im)


def Dyn(corr_map, labels, output_file="./dyneusr.html"):
    """
    DyNeuSR Visualization of the results

    Parameters
    ----------
    corr_map : matrix
        Analyzed correlation matrix .
    labels : numpy array
        Array of assigned clusters.
    output_file : str or path, optional
        Desired saving path for the generated html.
        The default is "./dyneusr.html".

    Returns
    -------
    None.

    """
    y = pd.get_dummies(labels.astype(str))

    mapper = KeplerMapper(verbose=2)
    X = corr_map
    # Configure projection
    pca = PCA(2, random_state=1)
    umap = UMAP(n_components=2, init=pca.fit_transform(X))

    # Construct lens and generate the shape graph
    lens = mapper.fit_transform(umap.fit_transform(X, y=None), projection=[0, 1])
    graph = mapper.map(
        lens,
        X=corr_map,
        cover=Cover(30, 0.4),
        clusterer=optimize_dbscan(corr_map, k=3, p=100.0),
    )

    # Convert to a DyNeuGraph
    dG = DyNeuGraph(G=graph, y=y)

    # Define some custom_layouts
    dG.add_custom_layout(lens, name="lens")
    dG.add_custom_layout(nx.spring_layout, name="nx.spring")
    dG.add_custom_layout(nx.kamada_kawai_layout, name="nx.kamada_kawai")
    dG.add_custom_layout(nx.spectral_layout, name="nx.spectral")
    dG.add_custom_layout(nx.circular_layout, name="nx.circular")

    # Configure some projections
    pca = PCA(2, random_state=1)
    tsne = TSNE(2, init="pca", random_state=1)
    umap = UMAP(n_components=2, init=pca.fit_transform(corr_map))

    # Add projections as custom_layouts
    dG.add_custom_layout(pca.fit_transform(X), name="PCA")
    dG.add_custom_layout(tsne.fit_transform(X), name="TSNE")
    dG.add_custom_layout(umap.fit_transform(X, y=None), name="UMAP")

    # Visualize
    dG.visualize(output_file)
