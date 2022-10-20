#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Processing library for clustintime
"""

import matplotlib.pyplot as plt  # For graphs

# Libraries
import numpy as np
import pandas as pd
import Visualization as vis
from scipy.signal import find_peaks


def RSS_peaks(corr_map, near):
    """
    Calculates the RSS of the correlation maps and returns the indexes of the time-points with the highest
    scores and the time-points nearby.

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
    RSS_1 = [np.sum(np.square(i)) ** 0.5 for i in corr_map]
    peaks = find_peaks(RSS_1)[0]
    new_peaks = [0] * 2 * near * len(peaks)

    for i in range(len(peaks)):
        if peaks[i] == 0:
            new_peaks[0 : 2 * near] = range(2 * near)
        elif peaks[i] == len(RSS_1):
            new_peaks[(len(RSS_1) - near) : (len(RSS_1) + near)] = range(
                (len(RSS_1) - near), (len(RSS_1) + near)
            )
        else:
            new_peaks[(peaks[i] - near) : (peaks[i] + near)] = range(
                (peaks[i] - near), (peaks[i] + near)
            )
    new_peaks = np.array(new_peaks)
    new_peaks = new_peaks[new_peaks != 0]
    new_peaks = np.array(list(set(new_peaks)))
    plt.plot(RSS_1)

    filtered_plot = np.array([np.mean(RSS_1)] * len(RSS_1))
    filtered_plot[new_peaks] = np.array(RSS_1)[new_peaks]

    plt.plot(filtered_plot)
    plt.title("RSS of original data vs filtered RSS")
    plt.legend(["Original RSS"], "Filtered RSS")

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
    triangle = np.triu(corr_map, 10)  # remove everything below the 20th diagonal

    thr_triangle = np.percentile(abs(triangle), thr)  # threshold the matrix
    aux = corr_map.copy()

    aux[abs(aux) < thr_triangle] = 0  # set to 0 everything below the thr

    zeroes = np.array([list(x).count(0) for x in aux])
    keep = np.array(np.where(zeroes < np.percentile(zeroes, thr)))[0]
    thresh_corr_map = pd.DataFrame(corr_map).loc[keep, keep]
    return np.matrix(thresh_corr_map), keep


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
            # The last rows will be concatenated (since there are less rows than the specified length once the
            # loop finishes, you can exit it)
            for window_idx in range(window_length):
                if (timepoint + window_idx + 1) < (data.shape[0]):
                    temp = np.concatenate((temp, data[timepoint + window_idx + 1, :]))
                else:
                    break
        concat_data[timepoint, range(len(temp))] = temp
    corr_map_window = np.corrcoef(concat_data, rowvar=True)

    return corr_map_window


def preprocess(
    corr_map, analysis, data=None, window_size=1, near=1, thr=95, contrast=1, task=[], TR=0.5
):
    """
    Main workflow for the processing algorithms

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    analysis : str
        Desired type of processing, the options are `double`, `thr`, `RSS`, `window`.
    data : matrix, optional
        fMRI data for the ẁindow` option. The default is None.
    window_size : int, optional
        Window size for the `window` option. The default is 1.
    near : int, optional
        Desired size of RSS peaks for the option `RSS`. The default is 1.
    thr : int, optional
        Threshold percentile for the option `thr`. The default is 95.
    contrast : int, optional
        Range of values of the correlation maps. The default is 1.
    task : dictionary or list, optional
        Structure containing the timings when the task is performed. The default is [].
    TR: float, optional
        TR of the data. The default is 0.5

    Returns
    -------
    corr_map : matrix
        Correlation map of the resulting processing.
    indexes : numpy array
        Array specifying the relevant time-points.

    """
    indexes = range(corr_map.shape[0])
    if analysis == "thr":
        aux, indexes = thr_index(corr_map, thr)
        vis.plot_two_matrixes(
            corr_map, aux, "Original matrix", "Filtered matrix", task, contrast, TR
        )
        corr_map, indexes = thr_index(corr_map, thr)
    elif analysis == "RSS":
        indexes = RSS_peaks(corr_map, near)
        vis.plot_two_matrixes(
            corr_map,
            pd.DataFrame(corr_map).loc[indexes, indexes],
            "Original matrix",
            "Filtered matrix",
            task,
            contrast,
            TR,
        )
        corr_map = pd.DataFrame(corr_map).loc[indexes, indexes]
    elif analysis == "double":
        vis.plot_two_matrixes(
            corr_map,
            np.nan_to_num(np.corrcoef(corr_map)),
            "Original matrix",
            "Double correlation matrix",
            task,
            contrast,
            TR,
        )
        corr_map = np.nan_to_num(np.corrcoef(corr_map))
    elif analysis == "window":
        vis.plot_two_matrixes(
            corr_map,
            correlation_with_window(data, window_size),
            "Original matrix",
            f"Correlation with sliding window size {window_size} ",
            task,
            contrast,
            TR,
        )
        corr_map = correlation_with_window(data, window_size)
    return corr_map, indexes
