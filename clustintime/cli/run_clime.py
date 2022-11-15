#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Parser for clustintime """

import argparse


def _get_parser():
    '''
    Parse command line inputs for clustintime

    Returns
    -------
    parser.parse_args() : argparse dict
    
    

    '''
    
    parser = argparse.ArgumentParser()
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('Required Argument:')
    
    required.add_argument('-i', '--input-file',
                          dest = 'data_file',
                          type = str,
                          help = 'The name or fullpath to the file containing the fMRI data',
                          required = True)
    
    required.add_argument('-m', '--mask-file',
                          dest = 'mask_file',
                          type = str,
                          help = 'The name or fullpath to the file containing the mask' 
                                  'for the data',
                          required = True)
    optional.add_argument('-com', '--component',
                          dest = 'component',
                          type = str,
                          help = 'Desired component of the signal to analyze, the options are `whole`,' 
                                 '`positive`, `negative`.',
                          required = False,
                          default = 'whole')
    
    optional.add_argument('-tf', '--timings-file',
                          dest = 'timings_file',
                          type = str,
                          help = 'The name or fullpath to the file(s) containing the onset' 
                                  'of the tasks',
                          required = False,
                          nargs = '*',
                          default = None)
    optional.add_argument('-cor', '--correlation',
                          dest = 'correlation',
                          type = str,
                          help = 'Desired type of correlation, the options are `standard`, `window`'
                                 'Default is `standard`',
                          required = False,
                          default = 'standard')  
    
    optional.add_argument('-p', '--processing',
                          dest = 'processing',
                          type = str,
                          help = 'The name of the desired type of processing.'
                                  'Default is None',
                          required = False,
                          default = None)     
    
    optional.add_argument('-ws', '--window-size',
                          dest = 'window_size',
                          type = int,
                          help = 'If \'-p\' is used and takes the value \'window\','
                                  ' it will become the window size of the processing option',
                          required = False,
                          default = 1)
    
    optional.add_argument('-n', '--near',
                          dest = 'near',
                          type = int,
                          help = 'If \'-p\' is used and takes the value \'RSS\','
                                  ' it will be used to indicate the number of points to use with the RSS',
                          required = False,
                          default = 1)
    
    optional.add_argument('-thr', '--threshold',
                          dest = 'thr',
                          type = int,
                          help = 'If \'-p\' is used and takes the value \'thr\', '
                                  'It will be used as a threshold percentile',
                          required = False,
                          default = 95)
    
    optional.add_argument('-c', '--contrast',
                          dest = 'contrast',
                          type = float,
                          help = 'Range of values for the correlation maps',
                          required = False,
                          default = 1)
    
    optional.add_argument('-tr', '--tr',
                          dest = 'TR',
                          type = float,
                          help = 'The name or fullpath to the file containing the mask' 
                                  'for the data',
                          required = False,
                          default = 0.5)
    
    optional.add_argument('-alg', '--algorithm',
                          dest = 'algorithm',
                          type = str,
                          help = 'Algorithm to be employed for clustering',
                          required = False,
                          default = 'infomap')
    optional.add_argument('-aff', '--affinity',
                          dest = 'affinity',
                          type = str,
                          help = 'Affinity to use. To see the available options see documentation of each algorithm'
                                 'Default is `euclidean`',
                          required = False,
                          default = 'euclidean') 
    
    optional.add_argument('-li', '--linkage',
                          dest = 'linkage',
                          type = str,
                          help = 'Linkage criterion for the \'Agglomerative Clustering\' algorithm. To see the available options see documentation of each algorithm'
                                 'Default is `ward`',
                          required = False,
                          default = 'ward') 
    
    optional.add_argument('-nc', '--n-clusters',
                          dest = 'n_clusters',
                          type = int,
                          help = 'Number of clusters for some sklearn algorithms',
                          required = False,
                          default = 7)
   
    
    optional.add_argument('-sm', '--save-maps',
                          dest = 'save_maps',
                          action = 'store_true',
                          help = 'Save the generated maps. ' 
                                  'Default is not to save',
                          required = False,
                          default = True)
    optional.add_argument('-con', '--consensus',
                          dest = 'consensus',
                          action = 'store_true',
                          help = 'Use consensus in the clustering algorithm.'
                                  'Default is not to use',
                          required = False,
                          default = '.')
    optional.add_argument('-sd', '--saving-dir',
                          dest = 'saving_dir',
                          type = str,
                          help = 'The name or fullpath to saving directory.'
                                  'Default is to save in the current directory',
                          required = False,
                          default = '.')
    
    optional.add_argument('-pre', '--prefix',
                          dest = 'prefix',
                          type = str,
                          help = 'Prefix for the saved data',
                          required = False,
                          default = '.')
    
    optional.add_argument('-s', '--seed',
                          dest = 'seed',
                          type = int,
                          help = 'Seed for the KMeans algorithm',
                          required = False,
                          default = 0)  
      
    optional.add_argument('-dyn', '--dyneusr',
                          dest = 'dyn',
                          action = 'store_true',
                          help = 'Generate and save DyneuSR map. ' 
                                  'Default is not to save',
                          required = False,
                          default = False)
    
    optional.add_argument('-f', '--fir',
                          dest = 'fir',
                          action = 'store_true',
                          help = 'Save the onset for FIR analysis in each cluster'
                                  'Default is not to save',
                          required = False,
                          default = '.')
    
    optional.add_argument('-t', '--title',
                          dest = 'Title',
                          type = str,
                          help = 'Title for the figures',
                          required = False,
                          default = '')
    parser._action_groups.append(optional)
    return parser

if __name__ == '__main__':
    
    raise RuntimeError('clustintime/cli/run_clime.py should not be run directly; \n'
                       'Please `pip install` clustintime and use the' 
                       '`clustintime` command')
   