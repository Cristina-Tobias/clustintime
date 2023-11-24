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