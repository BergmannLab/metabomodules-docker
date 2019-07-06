# Copyright 2019 Bergmann's Lab UNIL <mattia.tomasoni@unil.ch>
#
#    This file is part of METABOMODULES Tool.
#
#    METABOMODULES Tool is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    METABOMODULES Tool is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#
###############################################################################
# This function generates pseudospectra in input for metabomatching based on
# PCs - Principal Components.
###############################################################################

'''
@author: Mattia Tomasoni
'''

import csv
import numpy as np
from sklearn import decomposition

def import_data(input_file):
    input_reader = csv.reader(open(input_file, "r"), delimiter=",")
    features = np.array(list(input_reader)).astype("double")
    # extract first column (PPM labels)
    shifts = features[0:1,1:].transpose()
    features = np.delete(features, (0), axis=0) # delete PPM labels
    features = np.delete(features, (0), axis=1) # delete participant id labels
    return [shifts, features]

def calculate_pca_loadings(features, numb_pc_to_analyze):
    pca = decomposition.PCA(n_components=numb_pc_to_analyze, svd_solver='full')
    pca.fit(features).transform(features)
    loadings = pca.components_.transpose()
    return loadings

def save_pseudospectra(output_file, shifts, loadings, numb_pc_to_analyze):
    header = "shift"
    for pc_number in range(1, numb_pc_to_analyze+1):
        header = header + "\tPC/" + str(pc_number).zfill(3)
    with open(output_file, "w") as pseudospectrum:
        pseudospectrum.write(header+"\n")
    with open(output_file, "a") as pseudospectrum:
        data = np.hstack([shifts, loadings])
        np.savetxt(pseudospectrum, data, delimiter="\t")


def generate_pca_pseudospectra(input_file, output_file, numb_pc_to_analyze):
    [shifts, features] = import_data(input_file);
    loadings = calculate_pca_loadings(features, numb_pc_to_analyze);
    save_pseudospectra(output_file, shifts, loadings, numb_pc_to_analyze)    
