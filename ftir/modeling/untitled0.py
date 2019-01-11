#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 11:16:33 2018

@author: jyoung
"""

# Imports and pull in some data
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from ftir.io.utils import create_df_from_multiple_files
from ftir.modeling.buffer_subtraction import find_buffer_subtraction_constant, buffer_subtract
from ftir.modeling.peak_fitting import gaussian_least_squares, secondary_structure, sd_baseline_correction
from ftir.modeling.peak_definitions import yang_d20_2015, yang_h20_2015
#%%

location = '/home/jyoung/Desktop/windows-share/Current/Analytical_Labs/FTIR001/0748/20181217_CPW_BioATR/ANAL/2DER'

d2o_samples = ['Sample in D2O_0 2Der.txt', 'OVA in D2O_0 2Der.txt', 'OVA in D2O r2_0 2Der.txt']
h2o_samples = ['Protein in H2O_0 2Der.txt']

d2o_filenames = [location + '/' + i for i in d2o_samples]
h2o_filenames = [location + '/' + i for i in h2o_samples]

d2o_files = d2o_filenames
h2o_files = h2o_filenames

d2o_df, a = create_df_from_multiple_files(d2o_files)
h2o_df, a = create_df_from_multiple_files(h2o_files)
#%%
def rubberband(y, x, ascending=False):
    v = ConvexHull(np.column_stack([x, y])).vertices
    if ascending:
        v = np.roll(v, -v.argmin())
        v = v[:v.argmax()]
    else: 
        v = np.roll(v, -v.argmax())
        v = v[:v.argmin()]

    # Create baseline using linear interpolation between vertices
    return np.interp(x, x[v], y[v])

#%%

corr = rubberband(h2o_df['freq'], h2o_df['buffer'])