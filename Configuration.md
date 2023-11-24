# Configuration

There are several parameters and settings that can be customized. These settings control the behavior of the clustintime library and the clustering analysis.

Let's go through some of the main customizable configurations:

## Main Function Parameters (clustintime function):
- data_paths: Full path to the data to be analyzed.
- mask_path: Full path to the corresponding mask.
- component: Desired component of the signal to analyze (whole, positive, negative).
- timings_file: Path to .txt files containing timings of the analyzed task.
- correlation: Desired type of correlation (standard, window).
- process_type: Desired type of processing (None, double, thr, RSS, window).
- window_size: Window size for the window correlation option.
- near: Nearby time-points to select when performing RSS processing.
- thr: Threshold percentile for the thr processing.
- contrast: Range of values for the correlation matrices.
- repetition_time: Repetition time of the data.
- algorithm: Desired clustering algorithm for the analysis (infomap, Agglomerative, Louvain, Greedy, KMeans).
- consensus: Boolean indicating whether to use consensus clustering in the algorithm or not.
- n_clusters: Desired number of groups for the K Means algorithm.
- save_maps: Boolean indicating whether the results must be saved or not.
- saving_dir: Full path to the saving path.
- prefix: Prefix for the saved outcomes.
- seed: Seed for reproducibility.
- generate_dyneusr_graph: Generate a DyNeuSR graph.
- title: Title for the graphs.

## Command Line Interface (CLI) Parameters:
The command-line interface allows users to specify these parameters when running the script from the command line.

For example:
`````
python script.py --data_paths /path/to/data --mask_path /path/to/mask --algorithm infomap --n_clusters 5
`````
This command would run the script with the specified data and mask paths, using the Infomap clustering algorithm with 5 clusters.

## Other Configurations:
Within the code, there are several constants and default values that can be adjusted based on specific needs. For example, default values for certain parameters, file extensions (e.g., .nii.gz, .txt), and required Python version (python 3.6 or above) are set as constants in the script.

Overall, users can customize various aspects of the analysis by adjusting these parameters based on their specific data and research requirements.

## Example 1: Basic Clustering
Suppose you have fMRI data in the NIfTI format and want to perform basic clustering using the default settings:
`````
clustintime(
    data_paths='/path/to/your/data.nii.gz',
    mask_path='/path/to/your/mask.nii.gz',
    algorithm='infomap',
    n_clusters=7,
    save_maps=True,
    saving_dir='/path/to/save/results/',
    prefix='basic_clustering'
)
`````
In this example, the script will load the data, apply the Infomap clustering algorithm with 7 clusters, and save the results in the specified directory with a prefix.

## Example 2: Windowed Correlation
If you want to use a sliding window for correlation calculation:
`````
clustintime(
    data_paths='/path/to/your/data.nii.gz',
    mask_path='/path/to/your/mask.nii.gz',
    correlation='window',
    window_size=3,
    save_maps=True,
    saving_dir='/path/to/save/results/',
    prefix='windowed_correlation'
)
`````
This example uses a sliding window of size 3 for correlation calculation instead of the standard correlation.

## Example 3: Thresholding
For thresholding-based processing:
`````
clustintime(
    data_paths='/path/to/your/data.nii.gz',
    mask_path='/path/to/your/mask.nii.gz',
    process_type='thr',
    thr=90,
    save_maps=True,
    saving_dir='/path/to/save/results/',
    prefix='thresholding'
)
`````
This example applies thresholding at the 90th percentile.

## Example 4: KMeans Clustering with Consensus
Using KMeans clustering with consensus:
`````
clustintime(
    data_paths='/path/to/your/data.nii.gz',
    mask_path='/path/to/your/mask.nii.gz',
    algorithm='KMeans',
    consensus=True,
    n_clusters=5,
    save_maps=True,
    saving_dir='/path/to/save/results/',
    prefix='kmeans_with_consensus'
)
`````
This example applies KMeans clustering with consensus to improve stability.

## Example 5: DyNeuSR Graph Generation
Generating a DyNeuSR graph:
`````
clustintime(
    data_paths='/path/to/your/data.nii.gz',
    mask_path='/path/to/your/mask.nii.gz',
    generate_dyneusr_graph=True,
    save_maps=True,
    saving_dir='/path/to/save/results/',
    prefix='dyneusr_graph'
)
`````
This example generates a DyNeuSR graph in addition to the clustering results.

These examples demonstrate how you can customize the script for different clustering scenarios and data processing options. Adjust the parameters based on your specific analysis requirements.