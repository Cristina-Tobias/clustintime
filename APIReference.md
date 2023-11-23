# APIReference

This section provides a documentation of all the classes, functions and methods.

## clustintime Module Documentation

**_get_parser Function**

`````def _get_parser():
    """
    Parse command line inputs for clustintime

    Returns
    -------
    parser.parse_args() : argparse dict
    """
`````

This function defines a command-line argument parser for the clustintime module. It is used to parse user inputs when running clustintime from the command line.

## Arguments:

- i, --input-file: The name or full path to the file containing the fMRI data (required).
- m, --mask-file: The name or full path to the file containing the mask for the data (required).
- com, --component: Desired component of the signal to analyze, with options 'whole', 'positive', 'negative' (optional, default is 'whole').
- tf, --timings-file: The name or full path to the file(s) containing the onset of the tasks (optional, default is None).
- cor, --correlation: Desired type of correlation, with options 'standard', 'window' (optional, default is 'standard').
- p, --processing: The name of the desired type of processing (optional, default is None).
- ws, --window-size: Window size for processing when '-p' is 'window' (optional, default is 1).
- n, --near: Number of points to use with RSS when '-p' is 'RSS' (optional, default is 1).
- thr, --threshold: Threshold percentile when '-p' is 'thr' (optional, default is 95).
- c, --contrast: Range of values for the correlation maps (optional, default is 1).
- tr, --tr: Repetition time for the data (optional, default is 0.5).
- alg, --algorithm: Algorithm to be employed for clustering (optional, default is 'infomap').
- aff, --affinity: Affinity to use (optional, default is 'euclidean').
- li, --linkage: Linkage criterion for the 'Agglomerative Clustering' algorithm (optional, default is 'ward').
- nc, --n-clusters: Number of clusters for some sklearn algorithms (optional, default is 7).
- sm, --save-maps: Save the generated maps (optional, default is True).
- con, --consensus: Use consensus in the clustering algorithm (optional, default is False).
- sd, --saving-dir: The name or full path to the saving directory (optional, default is '.').
- pre, --prefix: Prefix for the saved data (optional, default is '.').
- s, --seed: Seed for the KMeans algorithm (optional, default is 0).
- dyn, --dyneusr: Generate and save DyneuSR map (optional, default is False).
- t, --title: Title for the figures (optional, default is '').

**__main__ Block**

`````
if __name__ == "__main__":
    raise RuntimeError(
        "clustintime/cli/run_clustintime.py should not be run directly; \n"
        "Please `pip install` clustintime and use the"
        "`clustintime` command"
    )
`````

This block of code ensures that the module is not run directly but is intended to be installed as part of the clustintime package. Users are encouraged to install clustintime and then use the clustintime command.

- Detailed documentation of all classes, functions, and methods.
- Include information on parameters, return values, and exceptions.
- Provide example usage for each API endpoint.