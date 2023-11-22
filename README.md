# Clustintime, a computational and visualization tool for time clustering of fMRI data.

Welcome! This toolbox contains all the functions required to run the whole project pipeline from one command.


## Motivation

In certain experimental contexts, particularly in naturalistic paradigms, clinical experiments (such as epilepsy), 
and during resting states, the precise timing of neuronal activity remains elusive. 
To address this challenge, researchers often turn to clustering algorithms in the time domain. 
This approach enables the identification of spatial patterns in neural activity by grouping neurons or brain 
regions together based on the similarity of their temporal dynamics.

The concept of similarity in time domain is pivotal to this methodology. It involves assessing the degree of
resemblance or corresponding between temporal patterns of different signals or data points.
When we talk about measuring similarity in the time domain, we are investigating how closely the temporal profiles
of 2 or more signals align over a specified period.

This analytical strategy becomes invaluable in scenarios where the temporal sequence of neuronal events is not 
predetermined. By employing clustering algorithms sensitive to temporal dynamics, researchers can uncover 
meaningful patterns in neural activity, shedding light on how different neurons or brain regions organize their 
responses over time. This approach is particularly relevant for understanding the intricacies of neural processing 
in real-world, uncontrolled conditions or clinical settings, where precise timing information may be challenging 
to obtain.

This toolbox aims to fill the gap in time domain analysis. ????

## Introduction

Spatial clustering is conventionally employed to identify patterns of synchronous BOLD fMRI signals across regions
(Thirion, B. 2014)

There is increasing interest in applying clustering algorithms in the temporal domain (Jo, Y. 2021)

### Key points of time domain clustering algorithms

1. Flexibility. They adapt to the inherent variability in timing, allowing the identification of patterns without
strict temporal contraints.
2. Pattern discovery that may have been overlooked in time-locked analysis.
3. Applicability across experimental designs, from exploratory research of spontaneous neural dynamics to controlled
tasks.

## Instructions

The toolbox is composed of 4 libraries:

1. func_for_clustering.py : This file contains the main pipeline
2. Visualization.py : This file has the necessary functions to visualise the results
3. Clustering.py : This file contains the algorithms employed in the pipeline and is programmed to show the results of the selected algorithm
4. Processing.py : This file has smoothing functions as well as filtering algorithms to clean the data


### func_for_clustering

The required arguments for this functions are:

1. The data directory
2. The mask directory

Optional arguments are:

1. Processing: this argument can take the values of 'double', 'RSS', 'thr' or 'window'
2. finger_tapping: this argument can take the values of 'No' or 'Yes'. If 'Yes', the directory of the timings and TR ought to be specified
3. window_size: necessary argument if Processing takes the value of 'window'
4. near: necessary argument if Processing takes the value of 'RSS'
5. thr: necessary argument if Processing takes the value of 'thr'
6. contrast: contrast for the correlation maps
7. dir_path: directory for the finger_tapping timings
8. repetition_time: Repetition Time (TR) of the data, only necessary in finger_tapping events
9. Algorithm: the possible values are 'infomap' or 'KMeans' 
10. thr_infomap: necessary argument if the selected algorithm is infomap
11. n_clusters: necessary argument if the selected algorithm is KMeans
12. save_maps: can take values of 'No' or 'Yes'
13. saving_dir: the results and plot will be saved onto the current directory if no other are specified

By default this library executes the pipeline with no processing and applies the infomap algorithm with a threshold of the 90th percentile

### Visualization

This library is composed of several functions for visualization.

1. display_table: displays an html table 
2. plot_labels: plots the time_series of the found clusters
3. RSS_peaks: plots the RSS of the data and the selected peaks. It return the indexes of the selected peaks
4. show_table: employs display_table to show a table with the clusters, the number of time-points per cluster and the relative importance of each cluster
5. plot_two_maps: plots the comparison of the original correlation map vs the one after preprocessing
6. Dyn: saves a DyNeuSR html of the results

### Clustering

This library contains the necessary functions to carry out the clustering operations:

1. save_maps: saves the results of the clustering into an specified directory
2. findCommunities: internal algorithm of infomap
3. K_Means: KMeans algorithm with display of results
4. Info_Map: InfoMap algorithm with display of results

### Processing

This library is employed by the main pipeline to apply filters in the data:

1. thr_index: employs a threshold based on percentile and removes those points below the threshold
2. correlation_with_window: calculates the correlation of the data using a sliding window of the specified length
3. preprocess: depending on the settings employed in the main pipeline, it will return a processed correlation map 


