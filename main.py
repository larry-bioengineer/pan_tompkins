# Larry To
# Created on: 11/18/2019
# Purpose: Re-create Pan-Tompkins algorithm for R-peak detections

# Import Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
import statistics as st
from scipy import signal

# Import Data
# df = pd.read_excel(r"C:\Users\larry\Dropbox\Learning Python\pan_tompkins\data.xlsx")
df = pd.read_excel(r"/Users/lokyiuto/Dropbox/Learning Python/pan_tompkins/data.xlsx", header=None)
# ecg = df.values.tolist()


ecg = df[0].values.tolist()
ecg = np.asarray(ecg)

fs = 1000.0
time = np.arange(0, len(ecg)/fs, 1/fs)


# plt.plot(time, ecg, label='unfiltered')
# plt.show()


# Remove mean
ecg = ecg - st.mean(ecg)

# Bandpass Filter
fc = np.array([5, 15])
fc = fc * 2 / fs


b, a = signal.butter(4, fc, 'bandpass')
ecg_f = signal.filtfilt(b, a, ecg)

plt.plot(time, ecg_f)
plt.show()
