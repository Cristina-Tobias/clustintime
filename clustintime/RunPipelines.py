#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 15:28:12 2021

@author: ctobias
"""

import clustintime as clime

data_directory = (
    "/bcbl/home/public/Epilepsy7T/study475/S41/S41vol.pySPFM.LARS.STC.betafitts_99.nii.gz"
)
mask_directory = "/bcbl/home/public/Epilepsy7T/study475/S41/S41vol.mask.nii.gz"
single_tap = "/bcbl/home/public/Epilepsy7T/study475/timings_SINGLETAP.1D"
multi_tap = "/bcbl/home/public/Epilepsy7T/study475/timings_MULTITAP.1D"
clime.clustintime(
    data_directory,
    mask_directory,
    single_tap,
    multi_tap,
    processing=None,
    timings=True,
    window_size=1,
    near=1,
    thr=95,
    contrast=1,
    TR=0.5,
    algorithm="infomap",
    thr_infomap=90,
    n_clusters=7,
    save_maps=False,
    saving_dir=".",
    prefix="Infomap",
)
