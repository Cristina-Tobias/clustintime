#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Processing library for clustintime
"""

import matplotlib.pyplot as plt  # For graphs

# Libraries
import numpy as np
from scipy.signal import find_peaks

class Processing:
    """ Main workflow for the processing algorithms

    Parameters
    ----------
    corr_map : matrix
        Correlation map of the data.
    analysis : str
        Desired type of processing, the options are `double`, `thr`, `RSS`.
    saving_dir : str or path
        Saving directory for figures
    prefix : str
        Prefix for figure names
    near : int, optional
        Desired size of RSS peaks for the option `RSS`. The default is 1.
    thr : int, optional
        Threshold percentile for the option `thr`. The default is 95.
    contrast : int, optional
        Range of values of the correlation maps. The default is 1.
    task : dictionary or list, optional
        Structure containing the timings when the task is performed. The default is [].
    repetition_time: float, optional
        Repetition time (TR) of the data. The default is 0.5
    Returns
    -------
    corr_map : matrix
        Correlation map of the resulting processing.
    indexes : numpy array
        Array specifying the relevant time-points.

    """
    def __init__(self, corr_map):
        self.corr_map = corr_map

        

    def rss_peaks(self, near):
        """
        Calculates the RSS of the correlation maps and returns the indexes of the time-points with the highest scores
        and the time-points nearby.
        Parameters
        ----------
        self.corr_map : matrix
            Correlation map of the data.
        near : int
            nearby points to the RSS peaks to be considered.
    
        Returns
        -------
        new_peaks : numpy array
            Array containing the indexes of the most relevant time-points.
    
        """
        rss_values = [np.sum(np.square(i)) ** 0.5 for i in self.corr_map]
        peaks = find_peaks(rss_values)[0]
        new_peaks = [0] * 2 * near * len(peaks)
    
        for peak in peaks:
            if peak <= near:
                new_peaks[0:near] = range(near)
            elif peak >= len(rss_values) - near:
                new_peaks[(len(rss_values) - near) : peak] = range((len(rss_values) - near), peak)
            else:
                new_peaks[(peak - near) : (peak + near)] = range((peak - near), (peak + near))
        new_peaks = np.array(new_peaks)
        new_peaks = new_peaks[new_peaks != 0]
        new_peaks = np.array(list(set(new_peaks)))
        plt.figure(figsize=[16, 8])
        plt.plot(rss_values)
        filtered_plot = np.array([np.mean(rss_values)] * len(rss_values))
        filtered_plot[[new_peaks]] = np.array(rss_values)[new_peaks]
        plt.plot(filtered_plot)
        plt.title("RSS of original data vs filtered RSS")
        plt.legend(["Original RSS", "Filtered RSS"])
        plt.xlabel("TR", fontsize=10)
        plt.ylabel("RSS", fontsize=10)
    
        return new_peaks
    
    
    def thr_index(self, thr):
        """
        Removes time-points that have a correlation under a certain threshold
    
        Parameters
        ----------
        self.corr_map : matrix
            Correlation map of the data.
        thr : int
            threshold percentile.
    
        Returns
        -------
        thresh_self.corr_map: matrix
            Correlation map containing only the relevant time-points.
        keep : numpy array
            Array containing the indexes of the relevant time-points.
    
        """
    
        if self.corr_map.max().max() == 0:
            self.corr_map[abs(self.corr_map) < np.percentile(abs(self.corr_map), thr)] = 0
        else:
            self.corr_map[self.corr_map < np.percentile(self.corr_map, thr)] = 0
        return self.corr_map
    
    

    
    
