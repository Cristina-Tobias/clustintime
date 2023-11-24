from clustintime.clustintime.clustintime import load_data
from clustintime.clustintime.clustintime import correlation_with_window
import pandas as pd

def test_load_data():
    data_file = "./test/data/100206_stabsel_spike_MvMEPFM_G50.THR_95.STC.DR2.nii.gz"
    mask_file = "./test/data/mask.nii.gz"
    data_masked, masker, nscans = load_data(data_paths=data_file, mask_paths=mask_file)


    assert data_masked.shape == (284,1936)
    assert str(type(masker)) == "<class 'nilearn.maskers.nifti_masker.NiftiMasker'>"
    assert nscans == [284]

def test_correlation_window():
    data_file = "./test/data/100206_stabsel_spike_MvMEPFM_G50.THR_95.STC.DR2.nii.gz"
    mask_file = "./test/data/mask.nii.gz"    
    data_masked, masker, nscans = load_data(data_paths=data_file, mask_paths=mask_file)
    obtained_matrix_correlation = correlation_with_window(data_masked, 5)
    obtained_matrix_correlation = pd.DataFrame(obtained_matrix_correlation)
    assert obtained_matrix_correlation.shape == (284,284)
    