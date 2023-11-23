#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clustintime is a python3 library meant to apply clustering algorithms over spatio-temporal fMRI data

It supports ``.nii.gz`` files as well as ``.txt`` files.

It requires python 3.6 or above, as well as the modules:
    - `numpy`
    - `matplotlib`
    - `pandas`
"""


import sys

# Libraries
import numpy as np
import pandas as pd
from nilearn.input_data import NiftiMasker

from clustintime.clustintime.clustintime.clustering import Clustering, generate_maps
from clustintime.clustintime.clustintime.cli.run_clustintime import _get_parser
from clustintime.clustintime.clustintime.consensus import Consensus
from clustintime.clustintime.clustintime.visualization import Visualization
from clustintime.clustintime.clustintime.processing import Processing


def load_data(data_paths, mask_paths):
    """
    Load and mask data with atlas using NiftiLabelsMasker.
    """
    # Initialize masker object
    masker = NiftiMasker(mask_img=mask_paths)

    # If n_pathss is 1 mask and return data
    if not isinstance(data_paths, list):
        data_masked = masker.fit_transform(data_paths)
        nscans = [data_masked.shape[0]]
    else:
        # If n_pathss is > 1, mask each paths in data list separately and
        # concatenate the masked data.
        for path_idx, path in enumerate(data_paths):
            if path_idx == 0:
                data_masked = masker.fit_transform(path)
                nscans = data_masked.shape[0]
            else:
                next_subject = masker.fit_transform(path)
                data_masked = np.concatenate((data_masked, next_subject), axis=0)
                nscans = np.append(nscans, next_subject.shape[0])

    return data_masked, masker, nscans

def implement_algorithm(
    algorithm,
    consensus,
    thr,
    n_clusters,
    nscans,
    corr_map,
    indices,
    seed,
    affinity,
    linkage,
    visualization_parameters,
    contrast,
):
    clustering_parameters = Clustering(corr_map, indices, nscans)
    if algorithm == "infomap":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering_parameters, algorithm, thr).find_clusters_with_consensus()
        else:
            corr_map, labels = clustering_parameters.info_map(thr=thr)
    elif algorithm == "KMeans":
        if consensus:
            labels = Consensus(clustering_parameters, algorithm, [n_clusters, seed]).find_clusters_with_consensus()
        else:
            labels = clustering_parameters.k_means(
                n_clusters=n_clusters, seed=seed
            )
    elif algorithm == "Agglomerative":
        if consensus:
            labels = Consensus(clustering_parameters, algorithm, 
                               [n_clusters, affinity, linkage]).find_clusters_with_consensus()
        else:
            labels = clustering_parameters.agglomerative_clustering(
                n_clusters=n_clusters,
                affinity=affinity,
                linkage=linkage,
            )
    elif algorithm == "Louvain":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering_parameters, algorithm, thr).find_clusters_with_consensus()
        else:
            corr_map, labels = clustering_parameters.louvain(thr)
    elif algorithm == "Greedy":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering_parameters, algorithm, thr).find_clusters_with_consensus()
        else:
            corr_map, labels = clustering_parameters.greedy_mod(thr)

    if algorithm in ["infomap", "Louvain", "Greedy"]:
        visualization_parameters.plot_two_matrices(
            corr_map_2, corr_map, "Original correlation map", "Binary correlation map", contrast=contrast
        )
    return labels

def preprocess(corr_map, analysis, near, thr):
    indices = range(corr_map.shape[0])
    if analysis == "thr":
        parameter = thr
        new_corr_map = Processing(corr_map).thr_index(thr)
    elif analysis == "RSS":
        indices = Processing(corr_map).rss_peaks(near)
        parameter = near
        new_corr_map = pd.DataFrame(corr_map).loc[indices, indices]
    return new_corr_map, corr_map, indices, parameter

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
    self.corr_map_window : matrix
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
            # The last rows will be concatenated (since there are less rows than the specified length
            # once the loop finishes, you can exit it)

            for window_idx in range(window_length):
                if (timepoint + window_idx + 1) < (data.shape[0]):
                    temp = np.concatenate((temp, data[timepoint + window_idx + 1, :]))
                else:
                    break
        concat_data[timepoint, range(len(temp))] = temp
    corr_map_window = np.corrcoef(concat_data, rowvar=True)

    return corr_map_window        

def clustintime(
    data_paths,
    mask_path,
    component="whole",
    timings_file=None,
    correlation="standard",
    process_type=None,
    window_size=1,
    near=1,
    thr=95,
    contrast=1,
    repetition_time=0.5,
    affinity="euclidean",
    linkage="ward",
    algorithm="infomap",
    consensus=False,
    n_clusters=7,
    save_maps=True,
    saving_dir=".",
    prefix="",
    seed=0,
    generate_dyneusr_graph=False,
    title="",
):
    """
    Run main workflow of clustintime.

    Estimates the functional connectivity of the data, processes the data and performs the clustering analysis,
    It returns summary on screen as well as the graphs corresponding to the data, the processing and the results.


    Parameters
    ----------
    data_paths : str or path
        Fullpath to the data to be analyzed.
    mask_path : str or path
        Fullpath to the corresponding mask.
    component : str, optional
        Desired component of the signal to analyze, the options are `whole`, `positive`, `negative`.
        The default is `whole`.
    timings_file : str or path, optional
        path to `.txt` files containing timings of the analyzed task.
        The default is `None`
    correlation : str, optional
        Desired type of correlation, the options are `standard`, `window`
        The default is `standard`
    process_type : str, optional
        Desired type of processing, the options are `None`, `double`, `thr`, `RSS`, `window`.
        The default is `None`.
    window_size : int, optional
        Window size for the `window` correlation option.
        The default is 1.
    near : int, optional
        Nearby time-points to select when performing `RSS` processing.
        The default is 1.
    thr : int, optional
        Threshold percentile for the `thr` processing.
        The default is 95.
    contrast : float, optional
        Range of values for the correlation matrixes.
        The default is 1.
    repetition_time : float, optional
        Repetition time of the data.
        The default is 0.5.
    algorithm : str, optional
        Desired clustering algorithm for the analysis, the options are `infomap`, `Agglomerative`, `Louvain`, `Greedy` and `KMeans`.
        The default is "infomap".
    consensus : bool, optional
        Boolean that indicates whether to use consensus clustering in the algorithm or not.
        The default is False.
    n_clusters : int, optional
        Desired number of groups for the K Means algorithm.
        The default is 7.
    save_maps : bool, optional
        Boolean that indicates whether the results must be saved or not.
        The default is True.
    saving_dir : str or path, optional
        Fullpath to the saving path.
        The default is ".".
    prefix: str, optional
        Prefix for the saved outcomes
    generate_dyneusr_graph : bool, optional
        Generate a DyNeuSR graph. Default is `False`
    title: str, optional
        Title for the graphs

    Returns
    -------
    None.

    """

    data, masker, nscans = load_data(data_paths, mask_path)


    if component == "negative":
        data[data > 0] = 0
    elif component == "positive":
        data[data < 0] = 0

    # Create data
    if timings_file is not None:
        # Load timings
        # 1D files are a txt file with the times in which the events occur. They are divided by the repetition_time.
        task = {}
        if isinstance(timings_file, str):
            timings_file = [timings_file]
        for idx, _dir in enumerate(timings_file):
            task[idx] = np.loadtxt(_dir)
    else:
        task = None

    if correlation == "standard":
        corr_map = np.nan_to_num(np.corrcoef(data))
    else:
        corr_map = np.nan_to_num(correlation_with_window(data, window_size))

    indexes = range(corr_map.shape[0])

    if process_type is not None:
        new_corr_map, corr_map, indexes, parameter = preprocess(
            corr_map=corr_map, analysis=process_type, near=near, thr=thr
        )
        
        Visualization(
            tasks=task,
            contrast=contrast,
            repetition_time=repetition_time,
            saving_dir=saving_dir,
            title=None,
            labels=None,
            prefix=f"{prefix}_orig_{process_type}_{parameter}",
        ).plot_two_matrixes(corr_map, new_corr_map, "Original matrix", "Filtered matrix", contrast)
        corr_map = new_corr_map

    visualization_parameters = Visualization(
        tasks=task,
        saving_dir=saving_dir,
        repetition_time=repetition_time,
        prefix=f"{prefix}_orig_binary",
        labels=None,
        title=title,
    )

    labels = implement_algorithm(
        algorithm,
        consensus,
        thr,
        n_clusters,
        nscans,
        corr_map,
        indexes,
        seed,
        affinity,
        linkage,
        visualization_parameters,
        contrast,
    )

    Visualization(
        tasks=task, saving_dir=saving_dir, repetition_time=repetition_time, prefix=prefix, labels=labels, title=title
    ).plot_heatmap(nscans)
    Visualization(
        tasks=task, saving_dir=saving_dir, repetition_time=repetition_time, prefix=prefix, labels=labels, title=title
    ).show_table()
    if save_maps:
        generate_maps(labels, saving_dir, data, masker, prefix)

    if generate_dyneusr_graph:
        Visualization(
            tasks=task,
            saving_dir=saving_dir,
            repetition_time=repetition_time,
            prefix=prefix,
            labels=labels,
            title=title,
        ).generate_dyneusr_visualization(corr_map)


def _main(argv=None):
    print(sys.argv)
    options = _get_parser().parse_args(argv)
    clustintime(**vars(options))
    if __name__ == "__main__":
        _main(sys.argv[1:])
