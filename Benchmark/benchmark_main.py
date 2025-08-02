#%%
from Accuracy import *
import pandas as pd
import numpy as np


#%%
# make fake dataset
other_model = pd.read_csv("fake_dataset/fake_other_models.csv", index_col=0)
other_model.head()
ccc_evts = ["ACE_BDKRB2", "ADAM10_APP", "ADAM10_AXL", "ADAM10_CD44", "ADAM10_NOTCH3", "ADGRE5_CD55", "ADIPOQ_ADIPOR2", "ADM_CALCRL", "ADM_RAMP1", "ADM_RAMP2", "ADM_RAMP3", "ANGPTL4_ITGAV_ITGB3"]

#for _ in range(10):


#