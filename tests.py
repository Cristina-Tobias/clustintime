from clustintime.pipeline import load_data
import numpy as np

def test_load_data():
    data_file = "C://Users//ferzt//Documents//Brainhack//clustintime//test//data//100206_stabsel_spike_MvMEPFM_G50.THR_95.STC.DR2.nii.gz"
    mask_file = "C://Users//ferzt//Documents//Brainhack//clustintime//test//data//mask.nii.gz"
    data_masked, masker, nscans = load_data(data_paths=data_file, mask_paths=mask_file)


    assert data_masked.shape == (284,1936)
    assert str(type(masker)) == "<class 'nilearn.maskers.nifti_masker.NiftiMasker'>"
    assert nscans == [284]

def test_stupid():
    prueba = np.ones([3,3])
    assert prueba.shape == (3,3)
    
# print('Hello')
test_stupid()
