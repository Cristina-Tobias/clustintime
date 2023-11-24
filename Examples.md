# Examples

Let's explore a variety of examples showcasing different use cases of the "Clustintime" toolbox:

## Example 1: Resting-State Network Identification
**Scenario:** Researchers want to identify resting-state networks in a group of healthy subjects.
`````
clustintime(
    data_paths='/path/to/resting_state_data.nii.gz',
    mask_path='/path/to/brain_mask.nii.gz',
    correlation='standard',
    algorithm='infomap',
    n_clusters=10,
    save_maps=True,
    saving_dir='/path/to/results/',
    prefix='resting_state_networks'
)
`````
**Benefit:** The toolbox applies the Infomap algorithm to identify clusters of brain regions with similar resting-state connectivity patterns, helping to delineate distinct resting-state networks.

## Example 2: Task-Related Connectivity Analysis
**Scenario:** Investigating task-related functional connectivity during a working memory task.
`````
clustintime(
    data_paths='/path/to/task_data.nii.gz',
    mask_path='/path/to/brain_mask.nii.gz',
    correlation='window',
    window_size=4,
    algorithm='KMeans',
    n_clusters=8,
    save_maps=True,
    saving_dir='/path/to/results/',
    prefix='task_connectivity'
)
`````
**Benefit:** Using a sliding window correlation approach, the toolbox applies KMeans clustering to identify dynamic changes in functional connectivity patterns related to the working memory task.

## Example 3: Comparative Study in Psychiatric Research
**Scenario:** Comparing functional connectivity between healthy controls and individuals with schizophrenia.
`````
clustintime(
    data_paths=['/path/to/controls_data.nii.gz', '/path/to/schizophrenia_data.nii.gz'],
    mask_path='/path/to/brain_mask.nii.gz',
    algorithm='Agglomerative',
    n_clusters=5,
    save_maps=True,
    saving_dir='/path/to/results/',
    prefix='psychiatric_comparison'
)
`````
**Benefit:** The toolbox applies agglomerative clustering to identify differences in functional connectivity patterns between the control and schizophrenia groups, potentially revealing biomarkers.

## Example 4: Longitudinal Study of Aging
**Scenario:** Analyzing longitudinal fMRI data to study age-related changes in brain connectivity.
`````
clustintime(
    data_paths='/path/to/longitudinal_data.nii.gz',
    mask_path='/path/to/brain_mask.nii.gz',
    algorithm='Louvain',
    n_clusters=6,
    save_maps=True,
    saving_dir='/path/to/results/',
    prefix='aging_study'
)
`````
**Benefit:** Louvain clustering is applied to track changes in functional connectivity patterns over time, providing insights into age-related alterations in brain networks.

## Example 5: Personalized Medicine - Individual Differences
**Scenario:** Investigating individual differences in brain connectivity for personalized medicine.
`````
clustintime(
    data_paths='/path/to/individual_data.nii.gz',
    mask_path='/path/to/brain_mask.nii.gz',
    algorithm='Greedy',
    consensus=True,
    n_clusters=3,
    save_maps=True,
    saving_dir='/path/to/results/',
    prefix='individual_differences'
)
`````
**Benefit:** The toolbox applies the Greedy algorithm with consensus clustering to identify stable individual differences in functional connectivity, potentially contributing to personalized treatment strategies.

These examples demonstrate the versatility of the "Clustintime" toolbox in various neuroimaging scenarios, from basic resting-state network identification to complex comparative studies and personalized medicine applications. Researchers can customize the toolbox for their specific use cases by adjusting parameters and algorithms based on the nature of their fMRI data and research questions.

# Integration into Existing Projects
Integrating the "Clustintime" toolbox into an existing project involves incorporating the necessary functions and scripts into your codebase and then calling the relevant functions with your specific data and parameters. 

Below is a step-by-step guide on how you can integrate the toolbox into an existing Python project:

## Step 1: Download the Toolbox
Download the "Clustintime" toolbox or clone the repository from the source. You can typically find this on a platform like GitHub.
`````
git clone https://github.com/username/clustintime.git
`````
## Step 2: Structure Your Project
Organize your existing project structure and create a directory for the "Clustintime" toolbox.
For example:
`````
your_project/
|-- clustintime/
|   |-- __init__.py
|   |-- clustering.py
|   |-- cli/
|   |   |-- __init__.py
|   |   |-- run_clustintime.py
|   |-- consensus.py
|   |-- processing.py
|   |-- visualization.py
|-- your_code.py
|-- your_data/
|-- ...
`````
## Step 3: Import Necessary Functions
In your Python script or module (your_code.py), import the necessary functions from the "Clustintime" toolbox.
For example:
`````
from clustintime.clustering import Clustering, generate_maps
from clustintime.consensus import Consensus
from clustintime.visualization import Visualization
from clustintime.processing import Processing
from clustintime.cli.run_clustintime import _get_parser, clustintime
`````
## Step 4: Configure Parameters
Configure the parameters based on your project requirements.
For instance:
`````
data_paths = '/path/to/your/data.nii.gz'
mask_path = '/path/to/your/mask.nii.gz'
algorithm = 'infomap'
n_clusters = 7
save_maps = True
saving_dir = '/path/to/save/results/'
prefix = 'your_project'
`````
## Step 5: Call the Clustintime Function
Call the clustintime function with your configured parameters:
`````
clustintime(
    data_paths=data_paths,
    mask_path=mask_path,
    algorithm=algorithm,
    n_clusters=n_clusters,
    save_maps=save_maps,
    saving_dir=saving_dir,
    prefix=prefix
)
`````
## Step 6: Execute Your Code
Execute your Python script or run your module, and the "Clustintime" toolbox functions will be integrated into your existing project. Adjust parameters, algorithms, and other settings as needed for your specific analysis.

## Example Integration:
`````
from clustintime.clustering import Clustering, generate_maps
from clustintime.consensus import Consensus
from clustintime.visualization import Visualization
from clustintime.processing import Processing
from clustintime.cli.run_clustintime import _get_parser, clustintime

# Configure Parameters
data_paths = '/path/to/your/data.nii.gz'
mask_path = '/path/to/your/mask.nii.gz'
algorithm = 'infomap'
n_clusters = 7
save_maps = True
saving_dir = '/path/to/save/results/'
prefix = 'your_project'

# Call Clustintime Function
clustintime(
    data_paths=data_paths,
    mask_path=mask_path,
    algorithm=algorithm,
    n_clusters=n_clusters,
    save_maps=save_maps,
    saving_dir=saving_dir,
    prefix=prefix
)
`````
This example assumes that you have the "Clustintime" toolbox in a directory within your project. Adjust the import statements and paths accordingly based on your specific project structure.

Now, you have successfully integrated the "Clustintime" toolbox into your existing project, and you can leverage its functionalities for analyzing fMRI data within the context of your research or application.

# Real World Scenarios
The "Clustintime" toolbox, designed for applying clustering algorithms to spatio-temporal fMRI data, can be beneficial in various real-world scenarios in the field of neuroimaging research.

Here are some examples:

## Brain Connectivity Studies:
- Scenario: Researchers often investigate functional connectivity patterns in the brain to understand how different regions communicate during specific tasks or at rest.
- Benefit: Clustering algorithms can help identify groups of brain regions with similar connectivity patterns, providing insights into functional networks and their dynamics.
## Task-Related Functional MRI (tfMRI) Analysis:
- Scenario: Studying brain activity during specific tasks or stimuli is common in cognitive neuroscience.
- Benefit: The toolbox allows researchers to analyze task-related fMRI data, identify functional clusters related to specific stimuli, and explore how these clusters evolve over time.
## Identification of Resting-State Networks:
- Scenario: Resting-state fMRI is used to study intrinsic brain activity in the absence of explicit tasks.
- Benefit: Clustering can help identify resting-state networks, revealing patterns of synchronized activity between different brain regions during resting conditions.
## Clinical Applications in Psychiatry:
- Scenario: Studying brain connectivity can be crucial in psychiatric research to understand abnormalities in patients with conditions like schizophrenia or depression.
- Benefit: Clustering can aid in identifying unique connectivity patterns associated with specific psychiatric disorders, potentially contributing to diagnostic and treatment strategies.
## Exploration of Dynamic Connectivity:
- Scenario: Brain connectivity is not static; it changes over time. Understanding dynamic connectivity patterns is essential.
- Benefit: The toolbox's ability to analyze dynamic connectivity, especially with methods like sliding window correlation, enables researchers to capture time-varying patterns in brain networks.
## Comparative Studies Across Populations:
- Scenario: Researchers may want to compare connectivity patterns between different populations, such as healthy controls and individuals with a neurological disorder.
- Benefit: Clustering helps identify differences in functional organization, potentially leading to biomarker discovery or a deeper understanding of neurological conditions.
## Longitudinal Studies:
- Scenario: Investigating changes in brain connectivity over time in longitudinal studies is common.
- Benefit: The toolbox allows for the analysis of longitudinal data, helping researchers track changes in functional connectivity patterns within individuals over extended periods.
## Exploration of Individual Differences:
- Scenario: Understanding how individuals differ in their brain connectivity can be valuable for personalized medicine.
- Benefit: Clustering can reveal individual variability in functional organization, aiding in the identification of personalized biomarkers or treatment approaches.

In summary, the "Clustintime" toolbox is beneficial in various neuroimaging scenarios, providing researchers with tools to analyze and interpret complex spatio-temporal fMRI data, ultimately advancing our understanding of brain function in health and disease.