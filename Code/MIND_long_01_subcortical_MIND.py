# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 09:44:18 2025

@author: natalia

Script to compute subcortical MIND networks

Copyright (C) 2026 University of Seville

Written by Natalia García San Martín (ngarcia1@us.es)

This file is part of Hierarchy Longitudinal MIND Gradients Psychosis toolkit.

Hierarchy Longitudinal MIND Gradients Psychosis toolkit is free software: 
you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or (at your option) any later version.

Hierarchy Longitudinal MIND Gradients Psychosis toolkit is distributed in the hope that 
it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Hierarchy Longitudinal MIND Gradients Psychosis toolkit. If not, see 
<https://www.gnu.org/licenses/>.

"""
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from scipy.spatial import cKDTree as KDTree
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats


"""
Utils
- functions in this section were taken from https://github.com/isebenius/MIND/blob/master/MIND_helpers.py
"""
def get_KDTree(x): #Inspired by https://gist.github.com/atabakd/ed0f7581f8510c8587bc2f41a094b518

    # Check the dimensions are consistent
    x = np.atleast_2d(x)
    
    # Build a KD tree representation of the samples
    xtree = KDTree(x)
    
    return xtree


def get_KL(x, y, xtree, ytree): #Inspired by https://gist.github.com/atabakd/ed0f7581f8510c8587bc2f41a094b518

    
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)

    n,d = x.shape
    m,dy = y.shape
    
    #Check dimensions
    assert(d == dy)

    # Get the first two nearest neighbours for x, since the closest one is the
    # sample itself.
    r = xtree.query(x, k=2, eps=.01, p=2)[0][:,1]
    s = ytree.query(x, k=1, eps=.01, p=2)[0]
    
    rs_ratio = r/s

    #Remove points with zero, nan, or infinity. This happens when two regions have a vertex with the exact same value – an occurence that basically onnly happens for the single feature MSNs
    #and has to do with FreeSurfer occasionally outputting the exact same value for different vertices.
    rs_ratio = rs_ratio[np.isfinite(rs_ratio)]
    rs_ratio = rs_ratio[rs_ratio!=0.0]
    
    # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
    # on the first term of the right hand side.

    kl = -np.log(rs_ratio).sum() * d / n + np.log(m / (n - 1.))
    kl = np.maximum(kl, 0)
    
    return kl


def calculate_mind_network(data_df, feature_cols, region_list, resample=False, n_samples = 4000):

    MIND = pd.DataFrame(np.zeros((len(region_list), len(region_list))), \
                        index = region_list, columns = region_list)

    #Get only desired regions
    data_df = data_df.loc[data_df['Label'].isin(region_list)]
    
    #Resample dataset if resample has been set to True and if it is UNIVARIATE ONLY. This should only be done if you are using a single feature which contains repeated values.
    if (len(feature_cols) == 1) and resample==True:   
        n_samples = n_samples
        resampled_dataset = pd.DataFrame(np.zeros((n_samples, len(region_list))), columns = region_list)

        for name, data in data_df.groupby('Label'):
            resampled_dataset[name] = stats.gaussian_kde(data[feature_cols[0]]).resample(n_samples)[0]

        resampled_dataset = resampled_dataset.melt(var_name = 'Label', value_name = feature_cols[0])
        data_df = resampled_dataset

    if (len(feature_cols) > 1) and resample==True:   
        raise Exception("Resampling the data is only supported if you are using a single feature -- this is because higher order density estimation can be unreliable and very computationally expensive.")

    #Check that there aren't many repeated values
    percent_unique_vals = len(data_df[feature_cols].drop_duplicates())/len(data_df[feature_cols])
    
    if percent_unique_vals < 0.8:
        raise Exception("There are many repeated values in the data, which compromises the validity of MIND calculation. Please minimize the number of repeated values in the data and try again. If you are using only one feature, try rerunning with resample=True.")

    grouped_data = data_df.groupby('Label')

    KDtrees = defaultdict(object)

    for i, (name_x, dat_x) in enumerate(grouped_data):
        tree = get_KDTree(dat_x[feature_cols])
        KDtrees[name_x] = tree
    
    used_pairs = []
    
    for i, (name_x, dat_x) in enumerate(grouped_data):
        
        for name_y, dat_y in grouped_data:
            if name_x == name_y:
                continue

            if set([name_x,name_y]) in used_pairs:
                continue

            dat_x = dat_x[feature_cols]
            dat_y = dat_y[feature_cols]

            KLa = get_KL(dat_x, dat_y, KDtrees[name_x], KDtrees[name_y])
            KLb = get_KL(dat_y, dat_x, KDtrees[name_y], KDtrees[name_x])

            kl = KLa + KLb

            MIND.at[name_x,name_y] = 1/(1+kl)
            MIND.at[name_y,name_x] = 1/(1+kl)

            used_pairs.append(set([name_x,name_y]))

    MIND = MIND[region_list].T[region_list].T
    
    return MIND



"""
scMIND
"""
sc_regions_dict = {
    'lh_thalamus': 10,
    'lh_caudate': 11,
    'lh_putamen': 12,
    'lh_pallidum': 13,
    'lh_hippocampus': 17,
    'lh_amygdala': 18,
    'lh_accumbens': 26,
    'rh_thalamus': 49,
    'rh_caudate': 50,
    'rh_putamen': 51,
    'rh_pallidum': 52,
    'rh_hippocampus': 53,
    'rh_amygdala': 54,
    'rh_accumbens': 58
}  


def number2region(index_value):
    parts = index_value.split('_')
    number = int(parts[1])
    sc_region = next((k for k, v in sc_regions_dict.items() if v == number), None)
    new_index_value = f'{sc_region}_{parts[2]}'
    return new_index_value


def prep_sc_data(sc_thickness_path, sc_logjacs_path): 
    sc_thickness = pd.read_csv(sc_thickness_path, index_col=0).T
    sc_thickness.columns = ['Feature_thickness']
    
    region_index_thickness = sc_thickness.index.map(number2region)
    sc_thickness.index = region_index_thickness
    
    sc_logjacs = pd.read_csv(sc_logjacs_path, index_col=0).T
    sc_logjacs.columns = ['Feature_area']
    
    region_index_logjacs = sc_logjacs.index.map(number2region)
    sc_logjacs.index = region_index_logjacs
    
    assert region_index_thickness.equals(region_index_logjacs)
    
    sc_data = pd.concat([sc_thickness, sc_logjacs], axis=1)
    
    sc_regions = sc_data.index.str.extract(r'(.+?)_\d+')[0]
    sc_data['Label'] = sc_regions.values
    cols = ['Label'] + [col for col in sc_data if col != 'Label']
    sc_data = sc_data[cols]
    
    sc_regions = np.array(list(sc_regions_dict.keys()))
    sc_features = ['Feature_thickness', 'Feature_area']
    
    return sc_data, sc_regions, sc_features


def compute_scMIND(logjacs_path, thickness_path):
    """
    Compute a subject's scMIND network.

    Parameters
    ----------
    logjacs_path : str
        Path to the subject's 'subj_LogJacs.csv' 
    thickness_path : str
        Path to the subject's 'subj_thick.csv' 

    Returns
    -------
    pandas.DataFrame
        Symmetric N×N similarity matrix (regions × regions) with region labels
        as both index and columns. 

    Example
    -------
    logjacs_path = '/home/lcs64/.../subj_LogJacs.csv'
    thickness_path = '/home/lcs64/.../subj_thick.csv'
    scMIND = compute_scMIND(logjacs_path, thickness_path)
    """
    
    sc_vertex_data, sc_regions, sc_features = prep_sc_data(thickness_path, logjacs_path)

    columns = ['Label'] + sc_features

    sc_vertex_data = sc_vertex_data[columns]

    for x in sc_features:
        sc_vertex_data[x] = (sc_vertex_data[x] - sc_vertex_data[x].mean())/sc_vertex_data[x].std()
    
    data_df = sc_vertex_data.reset_index(drop=True)
    
    scMIND = calculate_mind_network(data_df, sc_features, sc_regions)
    return scMIND



# NGS
saving_path = "C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/ENIGMA_Shape/"
sc_logJacs = pd.read_csv(saving_path+"groupfile_LogJacs.csv", index_col=0)

if not os.path.exists(saving_path + "sub-1006_ses-001.long.sub-1006_LogJacs.csv"):
    
    sc_thickness = pd.read_csv(saving_path+"groupfile_thick.csv", index_col=0)

    
    for patient in sc_logJacs.index:
        print(patient)
        file_logJacs = saving_path + f"{patient}_LogJacs.csv"
        file_thick = saving_path + f"{patient}_thick.csv"
        
        sc_logJacs.loc[[patient]].to_csv(file_logJacs)
        sc_thickness.loc[[patient]].to_csv(file_thick)



for patient in sc_logJacs.index:
    logJacs_path = saving_path + patient + "_LogJacs.csv"
    thickness_path = saving_path + patient + "_thick.csv"
    if not any(pd.isna(sc_logJacs.loc[patient])):
        scMIND = compute_scMIND(logJacs_path, thickness_path)
        scMIND.to_csv('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/MIND/MIND_networks_PAFIP/subcortical/'+patient+'.csv')
        print(patient)