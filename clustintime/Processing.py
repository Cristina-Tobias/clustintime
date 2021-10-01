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
from scipy.signal import find_peaks
import matplotlib.pyplot as plt # For graphs

def RSS_peaks(data, near):
    RSS_1 = [np.sum(np.square(i))**0.5 for i in data]
    peaks = find_peaks(RSS_1)[0]
    new_peaks = [0]*2*near*len(peaks)

    for i in range(len(peaks)):
        if peaks[i] == 0:
            new_peaks[0:2*near]=range(2*near)
        elif peaks[i] == len(RSS_1):
             new_peaks[(len(RSS_1)-near):(len(RSS_1)+near)] = range((len(RSS_1)-near),(len(RSS_1)+near))
        else:
            new_peaks[(peaks[i]-near):(peaks[i]+near)] = range((peaks[i]-near),(peaks[i]+near))
    new_peaks = np.array(new_peaks)
    new_peaks = new_peaks[new_peaks != 0]
    new_peaks = np.array(list(set(new_peaks)))   
    plt.plot(RSS_1)
    plt.plot(np.array(RSS_1)[new_peaks])
    plt.title('RSS of original data vs filtered RSS')
    plt.legend(['Original RSS'], 'Filtered RSS')
    
    return new_peaks

def thr_index(data, thr):
    triangle=np.triu(data, 10) #remove everything below the 20th diagonal

    thr_triangle = np.percentile(abs(triangle), thr) # threshold the matrix
    aux = data.copy()

    aux[abs(aux) < thr_triangle] = 0 #set to 0 everything below the thr

    zeroes = np.array([list(x).count(0) for x in aux])
    keep = np.array(np.where(zeroes < np.percentile(zeroes,thr)))[0]
    prueba = pd.DataFrame(data).loc[keep,keep]
    return np.matrix(prueba), keep

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

def preprocess(corr_map,analysis, data = None, window_size = 1, near = 1, thr = 95, contrast = 1, single_tap = [], multi_tap = []):
    indexes = range(corr_map.shape[0])
    if analysis == 'thr':
        aux, indexes = thr_index(corr_map, thr)
        vis.plot_two_maps(corr_map,aux ,'Original matrix', 'Filtered matrix',[], [], contrast)
        corr_map, indexes = thr_index(corr_map, thr)
    elif analysis == 'RSS':
        indexes = RSS_peaks(corr_map, near)
        vis.plot_two_maps(corr_map, pd.DataFrame(corr_map).loc[indexes,indexes],'Original matrix', 'Filtered matrix',single_tap, multi_tap, 0.25)
        corr_map = pd.DataFrame(corr_map).loc[indexes,indexes]            
    elif analysis == 'double':
        vis.plot_two_maps(corr_map, np.nan_to_num(np.corrcoef(corr_map)),'Original matrix', 'Double correlation matrix', single_tap, multi_tap, 0.75)
        corr_map = np.nan_to_num(np.corrcoef(corr_map))
    elif analysis == 'window':
        vis.plot_two_maps(corr_map, correlation_with_window(data, window_size),'Original matrix', f'Correlation with sliding window size {window_size} ', single_tap, multi_tap, 0.25)
        corr_map = correlation_with_window(data, window_size)
    return corr_map, indexes