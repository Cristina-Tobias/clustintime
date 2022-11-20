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


import os
import sys

# Libraries
import numpy as np
from nilearn.input_data import NiftiMasker
from nilearn.masking import apply_mask
from scipy.signal import find_peaks

import clustintime.clustering as clus
from clustintime.consensus import Consensus
import clustintime.processing as proc
import clustintime.visualization as vis
from clustintime.cli.run_clime import _get_parser


def load_data(data, mask_file):
    """
    Load and mask data with atlas using NiftiLabelsMasker.
    """
    # Initialize masker object
    masker = NiftiMasker(mask_img=mask_file)

    # If n_echos is 1 (single echo), mask and return data
    if len(data) == 1:
        # If data is a list, keep only first element
        if isinstance(data, list):
            data = data[0]
            data_masked = masker.fit_transform(data)
    else:
        # If n_echos is > 1 (multi-echo), mask each echo in data list separately and
        # concatenate the masked data.
        # If n_echos and len(data) are equal, read data.
        for echo_idx, echo in enumerate(data):
            if echo_idx == 0:
                data_masked = masker.fit_transform(echo)
            else:
                data_masked = np.concatenate((data_masked, masker.fit_transform(echo)), axis=0)
                #  If n_echos is different from len(data), raise error.

    return data_masked, masker


def clustintime(
    data_file,
    mask_file,
    component="whole",
    timings_file=None,
    correlation="standard",
    processing=None,
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
    fir=False,
    title="",
):
    """
    Run main workflow of clustintime.

    Estimates the functional connectivity of the data, processes the data and performs the clustering analysis,
    It returns summary on screen as well as the graphs corresponding to the data, the processing and the results.


    Parameters
    ----------
    data_file : str or path
        Fullpath to the data to be analyzed.
    mask_file : str or path
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
    processing : str, optional
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
        TR of the data.
        The default is 0.5.
    algorithm : str, optional
        Desired clustering algorithm for the analysis, the options are `infomap` and `KMeans`.
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
        Fullpath to the saving directory.
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

    # data_prueba = load_data(data_file, mask_file)
    # data = data_prueba[0]
    masker = NiftiMasker(mask_img=mask_file)
    masker.fit(data_file)

    # print("Applying mask ...")
    # print(" ")

    data = apply_mask(data_file, mask_file)  # apply mask to the fitted signal
    print("Mask applied!")

    if component == "negative":
        data[data > 0] = 0
    elif component == "positive":
        data[data < 0] = 0

    # Create data
    if timings_file is not None:
        # Load timings
        # 1D files are a txt file with the times in which the events occur. They are divided by TR
        task = {}
        if isinstance(timings_file, str):
            timings_file = [timings_file]
        for idx, _dir in enumerate(timings_file):
            task[idx] = np.loadtxt(_dir)

    else:
        task = []
    if correlation == "standard":
        corr_map = np.nan_to_num(np.corrcoef(data))
    else:
        corr_map = np.nan_to_num(proc.correlation_with_window(data, window_size))

    nscans = corr_map.shape[0]
    indexes = range(corr_map.shape[0])

    if processing is not None:
        corr_map, indexes = proc.preprocess(
            corr_map=corr_map,
            analysis=processing,
            near=near,
            thr=thr,
            contrast=contrast,
            task=task,
            repetition_time=repetition_time,
        )

    if algorithm == "infomap":
        corr_map_2 = corr_map.copy()
        if consensus:
            algorithm = clus.info_map
            labels = Consensus(algorithm, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            corr_map, labels = clus.info_map(corr_map, indexes, thr, nscans)

        vis.plot_two_matrixes(
            corr_map_2,
            corr_map,
            "Original correlation map",
            "Binary correlation map",
            tasks=task,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_binary",
            repetition_time=repetition_time,
            contrast=contrast,
        )
    elif algorithm == "KMeans":
        if consensus:
            algorithm = clus.k_means
            labels = clus.consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr)
        else:
            labels = clus.k_means(corr_map=corr_map, indexes=indexes, nscans=nscans, n_clusters=n_clusters, seed=seed)
    elif algorithm == "Agglomerative":
        if consensus:
            algorithm = clus.agglomerative_clustering
            labels = clus.consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr)
        else:
            labels = clus.agglomerative_clustering(
                corr_map=corr_map,
                indexes=indexes,
                nscans=nscans,
                n_clusters=n_clusters,
                affinity=affinity,
                linkage=linkage,
            )
    elif algorithm == "Louvain":
        corr_map_2 = corr_map.copy()
        if consensus:
            algorithm = clus.louvain
            labels = clus.consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr)
        else:
            corr_map, labels = clus.louvain(corr_map, indexes, thr, nscans=nscans)
        vis.plot_two_matrixes(
            corr_map_2,
            corr_map,
            "Original correlation map",
            "Binary correlation map",
            tasks=task,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_binary",
            repetition_time=repetition_time,
            contrast=contrast,
        )
    elif algorithm == "Greedy":
        corr_map_2 = corr_map.copy()
        if consensus:
            algorithm = clus.greedy_mod
            labels = clus.consensus(corr_map, indexes, nscans, n_clusters, algorithm, thr)
        else:
            corr_map, labels = clus.greedy_mod(corr_map, indexes, thr, nscans=nscans)
        vis.plot_two_matrixes(
            corr_map_2,
            corr_map,
            "Original correlation map",
            "Binary correlation map",
            tasks=task,
            saving_dir=saving_dir,
            prefix=f"{prefix}_orig_binary",
            repetition_time=repetition_time,
            contrast=contrast,
        )

    vis.plot_heatmap(labels, title, tasks=task, repetition_time=repetition_time, saving_dir=saving_dir, prefix=prefix)
    vis.show_table(labels, saving_dir, prefix)
    if save_maps:
        clus.generate_maps(labels, saving_dir, data, masker, prefix)

    if generate_dyneusr_graph:
        vis.generate_dyneusr_visualization(corr_map, labels, output_file=f"{saving_dir}/dyneusr_{prefix}.html")
    if fir:
        if os.path.exists(f"{saving_dir}/fir") == 0:
            os.mkdir(f"{saving_dir}/fir")
        for i in range(int(max(labels))):
            all_time_points = np.where(labels == i + 1)[0]
            difference = np.diff(all_time_points)
            select = find_peaks(difference)[0]
            fir_timepoints = np.insert(all_time_points[select + 1], 0, all_time_points[0])
            vis.plot_heatmap(
                labels,
                f"FIR onsets for cluster {i+1}",
                f"{saving_dir}/fir",
                f"{prefix}_fir_{i+1}",
                tasks=[fir_timepoints * repetition_time],
                repetition_time=repetition_time,
            )
            np.savetxt(f"{saving_dir}/fir/{prefix}_FIR_Cluster_{i+1}.1D", fir_timepoints * repetition_time)


def _main(argv=None):
    print(sys.argv)
    options = _get_parser().parse_args(argv)
    clustintime(**vars(options))
    if __name__ == "__main__":
        _main(sys.argv[1:])
