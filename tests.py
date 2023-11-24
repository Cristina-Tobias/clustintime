from clustintime.pipeline import load_data, preprocess
import numpy as np
import pandas as pd

def test_load_data():
    data_file = ".//test//data//100206_stabsel_spike_MvMEPFM_G50.THR_95.STC.DR2.nii.gz"
    mask_file = ".//test//data//mask.nii.gz"
    data_masked, masker, nscans = load_data(data_paths=data_file, mask_paths=mask_file)


    assert data_masked.shape == (284,1936)
    assert str(type(masker)) == "<class 'nilearn.maskers.nifti_masker.NiftiMasker'>"
    assert nscans == [284]

def test_preprocess_thr():
    corr_map=pd.read_csv('.//test//data//correlation2.csv').to_numpy()
    near=1
    thr=95
    
    def all_asserts():
        assert parameter== parameter_check
        if process_type == 'thr':
            #range object
            assert indexes == indexes_check
        assert np.sum(new_corr_map - new_corr_check) < 1e-10
        assert np.sum(corr_map_after - corr_map_check) < 1e-10
        print('all clear')
        
    process_type= "thr"

    new_corr_map, corr_map_after, indexes, parameter = preprocess(
        corr_map=corr_map, analysis=process_type, near=near, thr=thr
    )
    
    #thr
    indexes_check = range(284)
    parameter_check= 95
    new_corr_check= np.load(".//test//data//new_corr_map.npy")
    corr_map_check= np.load(".//test//data//corr_map_after_process.npy")
    
    all_asserts()
    
def test_preprocess_RSS():
    corr_map=pd.read_csv('.//test//data//correlation2.csv').to_numpy()
    near=1
    thr=95
    
    def all_asserts():
        assert parameter== parameter_check
            #array
        assert np.all(indexes==indexes_check)
        assert np.sum(new_corr_map - new_corr_check) < 1e-10
        assert np.sum(corr_map_after - corr_map_check) < 1e-10
    
    process_type= "RSS"
    new_corr_map, corr_map_after, indexes, parameter = preprocess(
        corr_map = corr_map, analysis=process_type, near=near, thr=thr
    )
    #RSS
    new_corr_map = new_corr_map.to_numpy()
    parameter_check = 1 
    indexes_check = np.load(".//test//data//indices_RSS.npy")
    new_corr_check = pd.read_csv(".//test//data//new_corr_map_RSS.csv").to_numpy()
    corr_map_check = np.load(".//test//data//corr_map_after_RSS.npy")
    
    all_asserts()
    