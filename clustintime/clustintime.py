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
from nilearn.input_data import NiftiMasker
from nilearn.masking import apply_mask


import clustintime.clustering as clustering
from clustintime.consensus import Consensus
import clustintime.processing as processing
from clustintime.visualization import Visualization
from clustintime.cli.run_clustintime import _get_parser


def load_data(data, mask_directory):
    """
    Load and mask data with atlas using NiftiLabelsMasker.
    """
    # Initialize masker object
    masker = NiftiMasker(mask_img=mask_directory)

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

def implement_algorithm(algorithm, consensus, thr, n_clusters, nscans, corr_map, indexes, seed, affinity, linkage, visualization_parameters, contrast):
    if algorithm == "infomap":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering.info_map, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            corr_map, labels = clustering.info_map(corr_map, indexes, thr, nscans=nscans)
    elif algorithm == "KMeans":
        if consensus:
            labels = Consensus(clustering.k_means, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            labels = clustering.k_means(corr_map=corr_map, indexes=indexes, nscans=nscans, n_clusters=n_clusters, seed=seed)
    elif algorithm == "Agglomerative":
        if consensus:
            labels = Consensus(clustering.agglomerative_clustering, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            labels = clustering.agglomerative_clustering(corr_map=corr_map, indexes=indexes, nscans=nscans, n_clusters=n_clusters, affinity=affinity,linkage=linkage)
    elif algorithm == "Louvain":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering.louvain, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            corr_map, labels = clustering.louvain(corr_map, indexes, thr, nscans=nscans)
    elif algorithm == "Greedy":
        corr_map_2 = corr_map.copy()
        if consensus:
            labels = Consensus(clustering.greedy_mod, thr, n_clusters, nscans).find_clusters_with_consensus(corr_map, indexes)
        else:
            corr_map, labels = clustering.greedy_mod(corr_map, indexes, thr, nscans=nscans)

    if algorithm  == "infomap" or "Louvain" or "Greedy":
        visualization_parameters.plot_two_matrixes(
            corr_map_2,
            corr_map,
            "Original correlation map",
            "Binary correlation map",
            contrast=contrast
        )
    return labels
    
    

def clustintime(
    data_directory,
    mask_directory,
    component="whole",
    timings_file=None,
    correlation="standard",
    proc=None,
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
    data_directory : str or path
        Fullpath to the data to be analyzed.
    mask_directory : str or path
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
  
    # data_prueba = load_data(data_directory, mask_directory)
    # data = data_prueba[0]
    masker = NiftiMasker(mask_img=mask_directory)
    masker.fit(data_directory)

    # print("Applying mask ...")
    # print(" ")

    data = apply_mask(data_directory, mask_directory)  # apply mask to the fitted signal
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

    if correlation == "standard":
        corr_map = np.nan_to_num(np.corrcoef(data))
    else:
        corr_map = np.nan_to_num(processing.correlation_with_window(data, window_size))

    nscans = corr_map.shape[0]
    indexes = range(corr_map.shape[0])

    if proc is not None:
        new_corr_map, corr_map, indexes, parameter = processing.preprocess(
            corr_map=corr_map,
            analysis=processing,
            near=near,
            thr=thr
        )
        Visualization(tasks=task, contrast=contrast,
                  repetition_time=repetition_time,
                  saving_dir=saving_dir, title=None,
                  labels=None, prefix=f"{prefix}_orig_{proc}_{parameter}").plot_two_matrixes(corr_map,
                                                                                                 new_corr_map,
                                                                                                 "Original matrix",
                                                                                                 "Filtered matrix",
                                                                                                 contrast)
        corr_map = new_corr_map
        
    visualization_parameters = Visualization(tasks=task, saving_dir=saving_dir, 
                                             repetition_time=repetition_time, 
                                             prefix=f"{prefix}_orig_binary", 
                                             labels=None, title=title)

    labels = implement_algorithm(algorithm, consensus, thr, n_clusters, nscans, corr_map, indexes, seed, affinity, linkage, visualization_parameters, contrast)
    
    Visualization(tasks=task, saving_dir=saving_dir, 
                                             repetition_time=repetition_time, 
                                             prefix=prefix, 
                                             labels=labels, title=title).plot_heatmap()
    Visualization(tasks=task, saving_dir=saving_dir, 
                                             repetition_time=repetition_time, 
                                             prefix=prefix, 
                                             labels=labels, title=title).show_table()
    if save_maps:
        clustering.generate_maps(labels, saving_dir, data, masker, prefix)

    if generate_dyneusr_graph:
        Visualization(tasks=task, saving_dir=saving_dir, 
                                             repetition_time=repetition_time, 
                                             prefix=prefix, 
                                             labels=labels, title=title).generate_dyneusr_visualization(corr_map)



def _main(argv=None):
    print(sys.argv)
    options = _get_parser().parse_args(argv)
    clustintime(**vars(options))
    if __name__ == "__main__":
        _main(sys.argv[1:])
