#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 17:00:49 2023

@author: crito25
"""

from clustintime.clustintime.clustintime import implement_algorithm
from clustintime.clustintime.visualization import Visualization
import pandas as pd
import numpy as np


def test_k_means():
    corr_map = pd.read_csv('./test/data/correlation.csv', index_col=0)
    expected_labels = np.load('./test/data/labels_kmeans.npy').reshape([284,])
    visualization_parameters = Visualization(
        tasks=None,
        saving_dir='.',
        repetition_time=0.5,
        prefix=".",
        labels=None,
        title=None,
    )
    indices= range(corr_map.shape[0])
    labels = implement_algorithm(
        algorithm='KMeans',
        consensus=False,
        thr=95,
        n_clusters=7,
        nscans=284,
        corr_map=corr_map,
        indices=indices,
        seed=0,
        affinity="euclidean",
        linkage="ward",
        visualization_parameters=visualization_parameters,
        contrast=1,
    )
    assert np.array_equal(labels,expected_labels)

def test_infomap():
    corr_map = pd.read_csv('./test/data/correlation.csv', index_col=0)
    expected_labels = np.load('./test/data/labels_infomap.npy').reshape([284,])
    visualization_parameters = Visualization(
        tasks=None,
        saving_dir='.',
        repetition_time=0.5,
        prefix=".",
        labels=None,
        title=None,
    )
    indices= range(corr_map.shape[0])
    labels = implement_algorithm(
        algorithm='infomap',
        consensus=False,
        thr=95,
        n_clusters=7,
        nscans=284,
        corr_map=corr_map,
        indices=indices,
        seed=0,
        affinity="euclidean",
        linkage="ward",
        visualization_parameters=visualization_parameters,
        contrast=1,
    )
    assert np.array_equal(labels,expected_labels)
    
def test_agglomerative():
    corr_map = pd.read_csv('./test/data/correlation.csv', index_col=0)
    expected_labels = np.load('./test/data/labels_agglomerative.npy').reshape([284,])
    visualization_parameters = Visualization(
        tasks=None,
        saving_dir='.',
        repetition_time=0.5,
        prefix=".",
        labels=None,
        title=None,
    )
    indices= range(corr_map.shape[0])
    labels = implement_algorithm(
        algorithm='Agglomerative',
        consensus=False,
        thr=95,
        n_clusters=7,
        nscans=284,
        corr_map=corr_map,
        indices=indices,
        seed=0,
        affinity="euclidean",
        linkage="ward",
        visualization_parameters=visualization_parameters,
        contrast=1,
    )
    assert np.array_equal(labels,expected_labels)

def test_louvain():
    corr_map = pd.read_csv('./test/data/correlation.csv', index_col=0)
    expected_labels = np.load('./test/data/labels_louvain.npy').reshape([284,])
    visualization_parameters = Visualization(
        tasks=None,
        saving_dir='.',
        repetition_time=0.5,
        prefix=".",
        labels=None,
        title=None,
    )
    indices= range(corr_map.shape[0])
    labels = implement_algorithm(
        algorithm='Louvain',
        consensus=False,
        thr=95,
        n_clusters=7,
        nscans=284,
        corr_map=corr_map,
        indices=indices,
        seed=0,
        affinity="euclidean",
        linkage="ward",
        visualization_parameters=visualization_parameters,
        contrast=1,
    )
    assert np.array_equal(labels,expected_labels)
def test_greedy():
    corr_map = pd.read_csv('./test/data/correlation.csv', index_col=0)
    expected_labels = np.load('./test/data/labels_greedy.npy').reshape([284,])
    visualization_parameters = Visualization(
        tasks=None,
        saving_dir='.',
        repetition_time=0.5,
        prefix=".",
        labels=None,
        title=None,
    )
    indices= range(corr_map.shape[0])
    labels = implement_algorithm(
        algorithm='Greedy',
        consensus=False,
        thr=95,
        n_clusters=7,
        nscans=284,
        corr_map=corr_map,
        indices=indices,
        seed=0,
        affinity="euclidean",
        linkage="ward",
        visualization_parameters=visualization_parameters,
        contrast=1,
    )
    assert np.array_equal(labels,expected_labels)

