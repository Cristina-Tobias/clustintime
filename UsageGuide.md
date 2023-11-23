# Usage Guide

## Core Components and Features:

1. Data Loading: The script can handle loading data from .nii.gz files as well as .txt files. It uses the Nilearn library's NiftiMasker for masking and loading data.

2. Clustering Algorithms: The toolbox supports multiple clustering algorithms, including Infomap, KMeans, Agglomerative, Louvain, and Greedy. Consensus clustering is an option for achieving robust results.

3. Data Preprocessing: Various preprocessing steps are available, such as thresholding, RSS (Region Sum of Squares), and windowed correlation. The correlation_with_window function calculates correlation using a sliding window.

4. Visualization: The toolbox provides visualization capabilities, including plotting correlation matrices, heatmaps, and tables. Options for saving the results and generating DyNeuSR graphs.

5. Command-Line Interface (CLI): The script includes a command-line interface for easy execution with command-line arguments.

6. Dependencies: The toolbox relies on external libraries, including NumPy, Pandas, Nilearn, and Matplotlib.

Include code snippets and examples for common use cases.

## How to Use:

Users can run the script from the command line, providing various options and parameters for data analysis.

The clustintime function is the main workflow, where users can specify data paths, clustering parameters, and other analysis options.

## Example Usage:

`````clustintime(
    data_paths="path/to/fmri_data.nii.gz",
    mask_path="path/to/mask.nii.gz",
    component="whole",
    correlation="standard",
    process_type=None,
    algorithm="infomap",
    consensus=False,
    n_clusters=7,
    save_maps=True,
    saving_dir="output/",
    prefix="result",
    seed=0,
    generate_dyneusr_graph=False,
    title="Clustintime Analysis",
)
`````

## How to Run:
The script can be executed from the command line with appropriate arguments.

Command-line arguments are parsed using the _get_parser function.

## Configuration Options and Parameters:

Refer to the Input section for more detail.

## Important Considerations and Best Practices 

1. Data Quality and Preprocessing:
Ensure that your input fMRI data and mask are of high quality and properly preprocessed.
Check for artifacts, motion effects, and other sources of noise.

2. Parameter Tuning:
Experiment with different parameter settings, especially for clustering algorithms.
Sensitivity to parameter choices may vary based on the characteristics of your data.

3. Consensus Clustering:
Understand the implications of using consensus clustering.
It can improve robustness but may also lead to increased computational costs.

4. Visualization:
Use visualization tools to interpret and validate your results.
Visualize correlation matrices, clustering results, and any additional relevant information.

5. Documentation:
Document your analysis workflow and parameter choices.
Include details on data sources, preprocessing steps, and any modifications made to the toolbox.

6. Reproducibility:
Set a random seed for reproducibility, especially if your analysis involves randomness (e.g., KMeans clustering).
Document the software versions and library dependencies.

7. Validation:
If possible, validate your clustering results using ground truth or external metrics.
Consider using known anatomical or functional regions for validation.

8. Parameter Sensitivity Analysis:
Conduct sensitivity analyses to evaluate how changes in parameters affect your results.
Assess the stability and reliability of the identified clusters.

9. Parallel Processing:
Depending on the size of your dataset, consider leveraging parallel processing to speed up computations.

10. Data Size:
Be mindful of the size of your data, as large datasets may require significant computational resources.

11. Task Timing:
If applicable, provide task timing files to enhance the specificity of your analysis.

12. Collaboration and Sharing:
If sharing your results or collaborating with others, provide clear documentation and share the version of the toolbox and libraries used.

13. Ethical Considerations:
If your analysis involves human subjects, ensure compliance with ethical standards and obtain necessary approvals.

14. Read the Documentation:
Thoroughly read the documentation of the "clustintime" toolbox to understand all available features and options.

15. Stay Informed:
Keep abreast of updates and releases for the "clustintime" toolbox or any similar tool.
Participate in relevant forums or communities for discussions and support.

