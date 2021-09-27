#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 09:50:31 2021

@author: ctobias
"""

# Libraries
import numpy as np
import pandas as pd
import Visualization as vis



def thr_index(data, thr):
    triangle=np.triu(data, 10) #remove everything below the 20th diagonal

    thr = np.percentile(abs(triangle), thr) # threshold the matrix

    triangle_thr = triangle.copy()
    triangle_thr[abs(triangle_thr) < thr] = 0 #set to 0 everything below the thr

    zeroes = np.array([list(x).count(0) for x in triangle_thr])
    keep = np.array(np.where(zeroes < np.percentile(zeroes,thr)))[0]
    prueba = pd.DataFrame(data).loc[keep,keep]
    return prueba

def correlation_with_window(data, window_length):
    concat_data = np.zeros((data.shape[0],data.shape[1]*(window_length + 1)))
    
    for timepoint in range(data.shape[0]):
        temp = data[timepoint, :] # Select a whole row
        
        if timepoint < (data.shape[0] - window_length): 
            # The first rows (those whose index is between nrow - window lenght) will follow this sequence:
            for window_idx in range(window_length): # the loop runs for the whole window length
                temp = np.concatenate((temp, data[timepoint + window_idx + 1, :])) # The data is concatenated for all the rows in the window
        else:
            # The last rows will be concatenated (since there are less rows than the specified length once the loop finishes, you can exit it)
            for window_idx in range(window_length):
                if (timepoint + window_idx + 1) < (data.shape[0]):
                    temp = np.concatenate((temp, data[timepoint + window_idx + 1, :]))
                else:
                    break
        concat_data[timepoint, range(len(temp)) ] = temp
    corr_map_window = np.corrcoef(concat_data, rowvar = True)

    return corr_map_window

def preprocess(corr_map,analysis, data, window_size = 1, near = 1, thr = 0.95, contrast = 1, single_tap = [], multi_tap = []):

    if analysis == 'thr':
        vis.plot_two_maps(corr_map, thr_index(corr_map, thr),'Original matrix', 'Filtered matrix',[], [], contrast)
        corr_map = thr_index(corr_map, thr)
    elif analysis == 'RSS':
        indexes = vis.RSS_peaks(corr_map, near)
        vis.plot_two_maps(corr_map, pd.DataFrame(corr_map).loc[indexes,indexes],'Original matrix', 'Filtered matrix',single_tap, multi_tap, 0.25)
        corr_map = pd.DataFrame(corr_map).loc[indexes,indexes]            
    elif analysis == 'double':
        vis.plot_two_maps(corr_map, np.nan_to_num(np.corrcoef(corr_map)),'Original matrix', 'Double correlation matrix', single_tap, multi_tap, 0.75)
        corr_map = np.nan_to_num(np.corrcoef(corr_map))
    elif analysis == 'window':
        vis.plot_two_maps(corr_map, correlation_with_window(data, window_size),'Original matrix', f'Correlation with sliding window size {window_size} ', single_tap, multi_tap, 0.25)
        corr_map = correlation_with_window(data, window_size)
    return corr_map