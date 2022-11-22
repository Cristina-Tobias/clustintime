#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Processing library for clustintime
"""

import matplotlib.pyplot as plt  # For graphs

# Libraries
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

import clustintime.visualization as visualization


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


def rss_peaks(corr_map, near):
    """
    Calculates the RSS of the correlation maps and returns the indexes of the time-points with the highest scores
    and the time-points nearby.

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    near : int
        Nearby points to the RSS peaks to be considered.

    Returns
    -------
    new_peaks : numpy array
        Array containing the indexes of the most relevant time-points.

    """
    rss_values = [np.sum(np.square(i)) ** 0.5 for i in corr_map]
    peaks = find_peaks(rss_values)[0]
    new_peaks = [0] * 2 * near * len(peaks)

    for peak in peaks:
        if peak <= near:
            new_peaks[0:near] = range(near)
        elif peak >= len(rss_values) - near:
            new_peaks[(len(rss_values) - near) : peak] = range((len(rss_values) - near), peak)
        else:
            new_peaks[(peak - near) : (peak + near)] = range((peak - near), (peak + near))
    new_peaks = np.array(new_peaks)
    new_peaks = new_peaks[new_peaks != 0]
    new_peaks = np.array(list(set(new_peaks)))
    plt.figure(figsize=[16, 8])
    plt.plot(rss_values)

    filtered_plot = np.array([np.mean(rss_values)] * len(rss_values))
    filtered_plot[[new_peaks]] = np.array(rss_values)[new_peaks]

    plt.plot(filtered_plot)
    plt.title("RSS of original data vs filtered RSS")
    plt.legend(["Original RSS", "Filtered RSS"])
    plt.xlabel("TR", fontsize=10)
    plt.ylabel("RSS", fontsize=10)

    return new_peaks


def thr_index(corr_map, thr):
    """
    Removes time-points that have a correlation under a certain threshold

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    thr : int
        Threshold percentile.

    Returns
    -------
    thresh_corr_map: matrix
        Correlation map containing only the relevant time-points.
    keep : numpy array
        Array containing the indexes of the relevant time-points.

    """

    if corr_map.max().max() == 0:
        corr_map[abs(corr_map) < np.percentile(abs(corr_map), thr)] = 0
    else:
        corr_map[corr_map < np.percentile(corr_map, thr)] = 0
    return corr_map


def correlation_with_window(data, window_length):
    """
    Calculates the correlation using a sliding window

    Parameters
    ----------
    data : matrix
        fMRI data.
    window_length : int
        size of the sliding window.

    Returns
    -------
    corr_map_window : matrix
        Correlation map of the data with the selected window.

    """
    concat_data = np.zeros((data.shape[0], data.shape[1] * (window_length + 1)))

    for timepoint in range(data.shape[0]):
        temp = data[timepoint, :]  # Select a whole row

        if timepoint < (data.shape[0] - window_length):
            # The first rows (those whose index is between nrow - window lenght) will follow this sequence:
            for window_idx in range(window_length):  # the loop runs for the whole window length
                temp = np.concatenate(
                    (temp, data[timepoint + window_idx + 1, :])
                )  # The data is concatenated for all the rows in the window
        else:
            # The last rows will be concatenated (since there are less rows than the specified length once the loop
            # finishes, you can exit it)
            for window_idx in range(window_length):
                if (timepoint + window_idx + 1) < (data.shape[0]):
                    temp = np.concatenate((temp, data[timepoint + window_idx + 1, :]))
                else:
                    break
        concat_data[timepoint, range(len(temp))] = temp
    corr_map_window = np.corrcoef(concat_data, rowvar=True)

    return corr_map_window


def preprocess(
    corr_map, analysis, saving_dir=".", prefix="", near=1, thr=95, contrast=1, task=None, repetition_time=0.5
):
    """
    Main workflow for the processing algorithms

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    analysis : str
        Desired type of processing, the options are `double`, `thr`, `RSS`.
    saving_dir : str or path
        Saving directory for figures
    prefix : str
        Prefix for figure names
    near : int, optional
        Desired size of RSS peaks for the option `RSS`. The default is 1.
    thr : int, optional
        Threshold percentile for the option `thr`. The default is 95.
    contrast : int, optional
        Range of values of the correlation maps. The default is 1.
    task : dictionary or list, optional
        Structure containing the timings when the task is performed. The default is [].
    repetition_time: float, optional
        TR of the data. The default is 0.5

    Returns
    -------
    corr_map : matrix
        Correlation map of the resulting processing.
    indexes : numpy array
        Array specifying the relevant time-points.

    """

    if task is None:
        task = []

    indexes = range(corr_map.shape[0])
    if analysis == "thr":
        aux = corr_map.copy()
        aux = thr_index(aux, thr)
        visualization.plot_two_matrixes(
            corr_map,
            aux,
            "Original matrix",
            "Filtered matrix",
            tasks=task,
            contrast=contrast,
            repetition_time=repetition_time,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_thr_{thr}",
        )
        corr_map = thr_index(corr_map, thr)
    elif analysis == "RSS":
        indexes = rss_peaks(corr_map, near)
        visualization.plot_two_matrixes(
            corr_map,
            pd.DataFrame(corr_map).loc[indexes, indexes],
            "Original matrix",
            "Filtered matrix",
            tasks=task,
            contrast=contrast,
            repetition_time=repetition_time,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_RSS_{near}",
        )
        corr_map = pd.DataFrame(corr_map).loc[indexes, indexes]
    elif analysis == "double":
        visualization.plot_two_matrixes(
            corr_map,
            np.nan_to_num(np.corrcoef(corr_map)),
            "Original matrix",
            "Double correlation matrix",
            tasks=task,
            contrast=contrast,
            repetition_time=repetition_time,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_double",
        )
        corr_map = np.nan_to_num(np.corrcoef(corr_map))

    return corr_map, indexes
