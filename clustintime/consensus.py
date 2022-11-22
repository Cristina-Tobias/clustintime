import random

import numpy as np
import pandas as pd
from clustintime.clustering import info_map, louvain, greedy_mod


def compute_connectivity_matrix(n_items, labels):
    connectivity_matrix = np.zeros([n_items, n_items])
    for j in range(n_items):
        if labels[j] > 0:
            row = np.where(labels == labels[j])
            connectivity_matrix[j, row] = 1
    return connectivity_matrix


def find_threshold_bfs(array):
    first_node = 0
    last_node = len(array) - 1
    probabilities = np.unique(array.ravel())
    low = 0
    high = len(probabilities)

    while high - low > 1:
        i = (high + low) // 2
        prob = probabilities[i]
        copied_array = np.array(array)
        copied_array[copied_array < prob] = 0.0
        if bfs(copied_array, first_node, last_node):
            low = i
        else:
            high = i

    return probabilities[low]


def bfs(graph, source, dest):
    """Perform breadth-first search starting at source. If dest is reached,
    return True, otherwise, return False."""
    # Based on http://www.ics.uci.edu/~eppstein/PADS/BFS.py
    # by D. Eppstein, July 2004.
    visited = set([source])
    nodes = np.arange(0, len(graph))
    stack = [(source, nodes[graph[source] > 0])]
    while stack:
        _, children = stack[0]
        for child in children:
            if child == dest:
                return True
            if child not in visited:
                visited.add(child)
                stack.append((child, nodes[graph[child] > 0]))
        stack.pop(0)
    return False


class Consensus:
    def __init__(self, algorithm, thr, n_clusters, n_scans):
        self.algorithm = algorithm
        self.threshold = thr
        self.n_clusters = n_clusters
        self.n_scans = n_scans

    def find_clusters_with_consensus(self, corr_map, indexes):
        npoints = len(indexes)
        sum_connectivity_matrix = np.zeros([npoints, npoints])
        index_matrix_sum = np.zeros([npoints, npoints])

        are_clusters_stable = False
        while not are_clusters_stable:
            for _ in range(100):
                sampling = np.sort(random.sample(range(npoints), round(npoints * 0.6)))
                index_matrix = pd.DataFrame([0] * npoints)
                index_matrix[0][sampling] = 1
                index_matrix_sum = index_matrix_sum + np.dot(index_matrix, np.transpose(index_matrix))
                data_sampled = corr_map[sampling, :][:, sampling]

                idx = self.get_indexes(data_sampled, indexes)

                idx = np.transpose(pd.DataFrame([idx, sampling]))
                idx = idx.set_index(1)
                idx = idx[np.logical_not(np.isnan(idx[0]))]
                labels = np.array([0] * npoints)
                labels[sampling] = idx[0]
                connectivity_matrix = compute_connectivity_matrix(npoints, labels)
            sum_connectivity_matrix = sum_connectivity_matrix + connectivity_matrix
            _consensus = np.divide(sum_connectivity_matrix, index_matrix_sum)
            _consensus[_consensus < find_threshold_bfs(_consensus)] = 0

            final_labels = self.get_labels(_consensus, indexes)

            thr = find_threshold_bfs(_consensus)
            _consensus[_consensus <= thr] = 0
            if self.algorithm in (info_map, greedy_mod, louvain):
                final_labels = final_labels[1]
            whole_connectivity_matrix = compute_connectivity_matrix(npoints, final_labels)

            are_clusters_stable = self.check_if_clusters_stable(
                _consensus, indexes, npoints, whole_connectivity_matrix
            )

        return final_labels

    def get_indexes(self, data_sampled, indexes):
        indexes = self.get_labels(data_sampled, indexes)
        if self.algorithm in (info_map, greedy_mod, louvain):
            return indexes[1]
        return indexes

    def get_labels(self, data_sampled, indexes):
        if self.algorithm in (info_map, greedy_mod, louvain):  # pylint: disable=comparison-with-callable
            return self.algorithm(data_sampled, indexes, self.threshold, self.n_scans)
        else:
            return self.algorithm(data_sampled, indexes, self.n_scans, self.n_clusters)

    def check_if_clusters_stable(self, _consensus, indexes, npoints, whole_connectivity_matrix):
        for _ in range(100):
            labels = self.get_labels(_consensus, indexes)
            if self.algorithm in (info_map, greedy_mod, louvain):
                labels = labels[1]
            connect = compute_connectivity_matrix(npoints, labels)

            if not np.array_equal(whole_connectivity_matrix, connect):
                return False

        return True
