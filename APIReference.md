# APIReference

## Overview

Clustintime is a Python library designed for applying clustering algorithms to spatio-temporal fMRI data. It supports both .nii.gz and .txt files and requires Python 3.6 or above. The required modules include:

- numpy
- matplotlib
- pandas

## Functions and Classes

**load_data** Function

`````
def load_data(data_paths, mask_paths):
    """
    Load and mask data with atlas using NiftiLabelsMasker.
    """
`````
This function loads and masks data with an atlas using NiftiLabelsMasker.

Arguments:
- data_paths: Str or path - Full path to the data to be analyzed.
- mask_paths: Str or path - Full path to the corresponding mask.

Returns:
- data_masked: Numpy array - Masked data.
- masker: Masker - Masker object.
- nscans: List - Number of scans.

**implement_algorithm** Function

`````
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
    """
    Perform clustering analysis using the specified algorithm.
    """
`````

This function performs clustering analysis using the specified algorithm, considering options such as consensus clustering, threshold, number of clusters, etc.

Parameters:
- algorithm (str): Desired clustering algorithm for the analysis.
- consensus (bool): Boolean that indicates whether to use consensus clustering in the algorithm.
- thr (int): Threshold percentile for processing.
- n_clusters (int): Desired number of clusters for KMeans algorithm.
- nscans (list): Number of scans.
- corr_map (numpy array): Correlation map of the data.
- indices (range): Range of indices.
- seed (int): Seed for the algorithm.
- affinity (str): Affinity parameter for agglomerative clustering.
- linkage (str): Linkage criterion for agglomerative clustering.
- visualization_parameters (Visualization): Visualization parameters.
- contrast (float): Range of values for correlation matrices.

Returns:
- labels (numpy array): Assigned clusters.

**preprocess** Function

`````
def preprocess(corr_map, analysis, near, thr):
    """
    Preprocess the correlation map based on the specified analysis.
    """
`````
This function preprocesses the correlation map based on the specified analysis, such as thresholding or selecting nearby time-points.

Parameters:
- corr_map (numpy array): Correlation map of the data.
- analysis (str): Type of analysis (e.g., "thr", "RSS").
- near (int): Nearby time-points for RSS processing.
- thr (int): Threshold percentile for "thr" processing.

Returns:
- new_corr_map (numpy array): Preprocessed correlation map.
- corr_map (numpy array): Original correlation map.
- indices (range): Range of indices.
- parameter (int): Parameter used for analysis.

**correlation_with_window** Function

`````
def correlation_with_window(data, window_length):
    """
    Calculates the correlation using a sliding window.
    """
`````

This function calculates the correlation using a sliding window on fMRI data.

Parameters:
- data (matrix): fMRI data.
- window_length (int): Size of the sliding window.

Returns:
- corr_map_window (matrix): Correlation map of the data with the selected window.

**clustintime** Function
`````
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
    Run the main workflow of clustintime.
    """
`````
This function runs the main workflow of clustintime, estimating functional connectivity, processing data, and performing clustering analysis.

Parameters:
- data_paths (str or path): Fullpath to the data to be analyzed.
- mask_path (str or path): Fullpath to the corresponding mask.
- component (str, optional): Desired component of the signal to analyze.
- timings_file (str or path, optional): Path to .txt files containing timings of the analyzed task.
- correlation (str, optional): Desired type of correlation.
- process_type (str, optional): Desired type of processing.
- window_size (int, optional): Window size for the correlation option.
- near (int, optional): Nearby time-points for processing.
- thr (int, optional): Threshold percentile for processing.
- contrast (float, optional): Range of values for the correlation matrix.
- repetition_time (float, optional): Repetition time of the data.
- algorithm (str, optional): Desired clustering algorithm for the analysis.
- consensus (bool, optional): Boolean indicating whether to use consensus clustering.
- n_clusters (int, optional): Desired number of groups for the K Means algorithm.
- save_maps (bool, optional): Boolean indicating whether to save the results.
- saving_dir (str or path, optional): Fullpath to the saving path.
- prefix (str, optional): Prefix for the saved outcomes.
- seed (int, optional): Seed for the algorithm.
- generate_dyneusr_graph (bool, optional): Generate a DyNeuSR graph.
- title (str, optional): Title for the graphs.

## Command-Line Interface (CLI)
The script also includes a command-line interface (CLI) provided by the _main function. This CLI accepts various command-line arguments for running clustintime with different configurations. Run the script from the command line with appropriate arguments.

## Note
This documentation provides an overview of the functions and classes in the clustintime script. For more detailed information on each function's parameters and behavior, refer to the code and comments within the script.

- Provide example usage for each API endpoint.