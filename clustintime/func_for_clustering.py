#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 12:56:04 2021

@author: ctobias
"""
# Libraries
import numpy as np
import os # access system's identity
import Processing as proc
import Clustering as clus
import Visualization as vis

from nilearn.masking import apply_mask 
from nilearn.input_data import NiftiMasker


def pipeline(data_directory, mask_directory, processing = None, finger_tapping = 'No',  window_size = 1, near = 1, thr = 95, contrast = 1, dir_path = '.', TR = 0.5, algorithm = 'infomap', thr_infomap = 90, n_clusters = 7, save_maps = 'No', saving_dir = '.'):
    masker = NiftiMasker(mask_img = mask_directory)
    masker.fit(data_directory)
    
    print('Applying mask ...')
    print(' ')
    
    data = apply_mask(data_directory,mask_directory) # apply mask to the fitted signal    
    
    print('Mask applied!')
        # Create data
    if finger_tapping == 'Yes':
        # Load timings
        # 1D files are a txt file with the times in which the events occur. They are divided by TR
        single_tap = np.loadtxt(os.path.join(dir_path, 'timings_SINGLETAP.1D')) / TR
        multi_tap = np.loadtxt(os.path.join(dir_path, 'timings_MULTITAP.1D')) / TR
    else:
        single_tap = []
        multi_tap = []

    corr_map = np.nan_to_num(np.corrcoef(data))
    nscans = corr_map.shape[0]
    indexes = range(corr_map.shape[0])
    if processing != None:
        corr_map, indexes = proc.preprocess(corr_map = corr_map, analysis = processing, data = data, window_size = window_size, near = near, thr = thr, contrast = contrast, single_tap = single_tap, multi_tap = multi_tap)    
    if algorithm == 'infomap':
       corr_map, labels = clus.Info_Map(data, indexes,thr_infomap, nscans = nscans ,multi_tap = multi_tap, single_tap = single_tap)
    elif algorithm == 'KMeans':
       labels = clus.K_Means(corr_map=corr_map,indexes = indexes, nscans = nscans ,n_clusters=n_clusters, multi_tap = multi_tap, single_tap = single_tap)
    if save_maps == 'Yes':
        clus.save_maps(labels, saving_dir, data_directory, mask_directory)
    vis.Dyn(corr_map, labels,output_file = f'{saving_dir}/dyneusr_{algorithm}.html')
 
    

 
        
    
