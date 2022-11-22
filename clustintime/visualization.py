#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Result visualization for clustintime
"""

# import networkx as nx
# creation, manipulation and study of the structure, dynamics and functions of complex networks

from matplotlib import patches
import matplotlib.pyplot as plt  # For graphs
import networkx as nx

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


class Visualization:
    def __init__(self, title, saving_dir, prefix, tasks, repetition_time, labels):

        if tasks is None:
            tasks = []

        self.title = title
        self.saving_dir = saving_dir
        self.prefix = prefix
        self.tasks = tasks
        self.repetition_time = repetition_time
        self.labels = labels

    def plot_heatmap(self):
        """
        Visualization of the clusters separately

        Parameters
        ----------
        labels : numpy array
            Array of assigned clusters.
        title : str
            Title for the plot.
        saving_dir : str or path
            Saving directory for the plot
        prefix : str
            Prefix for the image
        task : dictionary or list, optional
            Structure containing the times when the task is performed. The default is [].
        repetition_time : float, optional
            TR of the data. The default is 0.5.

        Returns
        -------
        None.

        """

        plt.figure(figsize=[8, 8])
        heatmatrix = np.zeros([int(self.labels.max()), len(self.labels)])
        rownames = np.zeros([int(self.labels.max())]).astype(str)
        x_values = np.linspace(0, len(self.labels) * self.repetition_time, len(self.labels)).astype(int)
        file1d = pd.DataFrame()
        for i in range(int(self.labels.max())):
            selected_labels = np.array([0] * len(self.labels))
            selected_labels[np.where(self.labels == i + 1)] = 1
            file1d = pd.concat([file1d, pd.DataFrame(selected_labels)], axis=1)
            heatmatrix[i] = selected_labels
            rownames[i] = f"# {i+1}"

        heatmatrix = pd.DataFrame(heatmatrix, columns=x_values, index=rownames)
        colors = sns.color_palette("Dark2", len(self.tasks) + 1)
        sns.heatmap(heatmatrix, cmap="Greys", xticklabels=150, cbar=False)
        plt.xlabel("Time in seconds", fontsize=10)
        plt.ylabel("Clusters", fontsize=10)
        legends = np.zeros([len(self.tasks)]).astype(str)
        rectangles = [patches.Rectangle((0, 0), 1, 1, facecolor=colors[0])]
        for idx, task in enumerate(self.tasks.values()):
            plt.vlines(task / self.repetition_time, 0, self.labels.max(), linewidth=1.2, colors=colors[idx], alpha=0.5)
            legends[idx] = f"task {idx}"
            rectangles.append(patches.Rectangle((0, 0), 1, 1, facecolor=colors[idx + 1]))
        plt.title(self.title)
        plt.legend(
            (rectangles),
            np.array(legends),
            bbox_to_anchor=[1, 0.5],
            loc="center left",
            handlelength=1,
            handleheight=1,
        )
        plt.savefig(f"{self.saving_dir}/{self.prefix}_heatmap.png")
        np.savetxt(f"{self.saving_dir}/{self.prefix}_heatmap.1D", file1d)

    def show_table(self):
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
        cluster, count = np.unique(self.labels, return_counts=True)
        array = [cluster, count, count / len(self.labels)]
        table_result = pd.DataFrame({"Cluster": array[0], "Count": array[1], "Percentage": array[2]})
        print(table_result)
        table_result.to_csv(f"{self.saving_dir}/{self.prefix}_Results.csv")

    def plot_two_matrices(self, map_1, map_2, title_1, title_2, contrast=1):
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
        tasks : dictionary or list, optional
            Structure containing the times when the task is performed. The default is [].
        contrast : int, optional
            Range of values of the correlation matrixes. The default is 1.
        repetition_time: float, optional
            TR of the data. The default is 0.5

        Returns
        -------
        None.

        """
        if self.tasks is None:
            self.tasks = []

        fig = plt.figure(figsize=(16, 8))
        grid_spec = fig.add_gridspec(1, 2, hspace=0)

        (ax1, ax2) = grid_spec.subplots()
        ax1.imshow(map_1, aspect="equal", vmin=-contrast, vmax=contrast, cmap="RdBu_r")
        # Vertical lines to delimit the instant in which the event occurs
        limit = map_1.shape[0]
        for task in self.tasks.values():
            ax1.vlines(task / self.repetition_time, ymin=0, ymax=limit, linewidth=0.7)
            ax1.hlines(task / self.repetition_time, xmin=0, xmax=limit, linewidth=0.7)  # Same for horizontal
        ax1.set_title(title_1)
        ax1.set_xlim([0, limit])
        ax1.set_ylim([limit, 0])

        # Plot

        image = ax2.imshow(map_2, aspect="equal", vmin=-contrast, vmax=contrast, cmap="RdBu_r")
        limit = map_2.shape[0]
        for task in self.tasks.values():
            ax2.vlines(task / self.repetition_time, ymin=0, ymax=limit, linewidth=0.7)
            ax2.hlines(task / self.repetition_time, xmin=0, xmax=limit, linewidth=0.7)  # Same for horizontal
        ax2.set_title(title_2)
        ax2.set_xlim([0, limit])
        ax2.set_ylim([limit, 0])
        fig.colorbar(image)
        plt.savefig(f"{self.saving_dir}/{self.prefix}_matrix_comparison.png")

    def generate_dyneusr_visualization(self, corr_map):
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
        y_values = pd.get_dummies(self.labels.astype(str))

        mapper = KeplerMapper(verbose=1)
        # Configure projection
        pca = PCA(2, random_state=1)
        umap = UMAP(n_components=2, init=pca.fit_transform(corr_map))

        # Construct lens and generate the shape graph
        lens = mapper.fit_transform(umap.fit_transform(corr_map, y=None), projection=[0, 1])
        graph = mapper.map(
            lens,
            X=corr_map,
            cover=Cover(20, 0.7),
            clusterer=optimize_dbscan(corr_map, k=3, p=100.0),
        )

        # Convert to a DyNeuGraph
        dyneusr_graph = DyNeuGraph(G=graph, y=y_values)

        # Define some custom_layouts
        dyneusr_graph.add_custom_layout(lens, name="lens")
        dyneusr_graph.add_custom_layout(nx.spring_layout, name="nx.spring")
        dyneusr_graph.add_custom_layout(nx.kamada_kawai_layout, name="nx.kamada_kawai")
        dyneusr_graph.add_custom_layout(nx.spectral_layout, name="nx.spectral")
        dyneusr_graph.add_custom_layout(nx.circular_layout, name="nx.circular")

        # Configure some projections
        pca = PCA(2, random_state=1)
        tsne = TSNE(2, init="pca", random_state=1)
        umap = UMAP(n_components=2, init=pca.fit_transform(corr_map))

        # Add projections as custom_layouts
        dyneusr_graph.add_custom_layout(pca.fit_transform(corr_map), name="PCA")
        dyneusr_graph.add_custom_layout(tsne.fit_transform(corr_map), name="TSNE")
        dyneusr_graph.add_custom_layout(umap.fit_transform(corr_map, y=None), name="UMAP")

        # Visualize
        output_file = f"{self.saving_dir}_{self.prefix}_dyneusr.html"
        dyneusr_graph.visualize(output_file)
