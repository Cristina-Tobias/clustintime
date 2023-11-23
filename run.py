from clustintime import pipeline
import matplotlib.pyplot as plt


data_file = "C://Users//ferzt//Documents//Brainhack//clustintime//test//data//100206_stabsel_spike_MvMEPFM_G50.THR_95.STC.DR2.nii.gz"
mask_file = "C://Users//ferzt//Documents//Brainhack//clustintime//test//data//mask.nii.gz"

pipeline.clustintime(data_file,mask_file)

plt.show()