import random

import numpy as np
import pandas as pd
from clustintime.clustering import info_map, louvain, greedy_mod

import clustintime.processing as proc


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

                idx = np.transpose(
                    pd.DataFrame([idx, sampling])
                )
                idx = idx.set_index(1)
                idx = idx[np.logical_not(np.isnan(idx[0]))]
                labels = np.array([0] * npoints)
                labels[sampling] = idx[0]
                connectivity_matrix = proc.compute_connectivity_matrix(npoints, labels)
            sum_connectivity_matrix = sum_connectivity_matrix + connectivity_matrix
            _consensus = np.divide(sum_connectivity_matrix, index_matrix_sum)
            _consensus[_consensus < proc.find_threshold_bfs(_consensus)] = 0

            final_labels = self.get_labels(_consensus, indexes)

            thr = proc.find_threshold_bfs(_consensus)
            _consensus[_consensus <= thr] = 0
            whole_connectivity_matrix = proc.compute_connectivity_matrix(npoints, final_labels[1])

            are_clusters_stable = self.check_if_clusters_stable(_consensus, indexes, npoints, whole_connectivity_matrix)

        return final_labels[1]

    def get_indexes(self, data_sampled, indexes):
        labels = self.get_labels(data_sampled, indexes)
        if self.algorithm in (info_map, greedy_mod, louvain):
            return labels[1]

        return labels

    def get_labels(self, data_sampled, indexes):
        if self.algorithm in (info_map, greedy_mod, louvain):  # pylint: disable=comparison-with-callable
            return self.algorithm(data_sampled, indexes, self.threshold, self.n_scans)
        else:
            return self.algorithm(data_sampled, indexes, self.n_scans, self.n_clusters)

    def check_if_clusters_stable(self, _consensus, indexes, npoints, whole_connectivity_matrix):
        for _ in range(100):
            labels = self.algorithm(corr_map=_consensus, indexes=indexes, thr=self.threshold, nscans=npoints)
            connect = proc.compute_connectivity_matrix(npoints, labels[1])

            if not np.array_equal(whole_connectivity_matrix, connect):
                return False

        return True