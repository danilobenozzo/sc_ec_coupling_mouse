#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 11:23:22 2022

@author: benozzo
"""

import numpy as np
import sklearn.metrics as metrics
import pdb

def get_idx_network_subSeparated():
    
    idx_network = np.ones([37,1], dtype=int)
    
    idx_lcn = [0,1,2,3]
    idx_network[idx_lcn] = [[1] for _ in idx_lcn] #LCN 
    
    idx_vis = [5]
    idx_network[idx_vis] = [[-1] for _ in idx_vis] #VIS

    idx_aud = [6]
    idx_network[idx_aud] = [[-1] for _ in idx_aud] #AUD

    idx_dmn_mid = [7,8,9,10,11,15,16,17]
    idx_network[idx_dmn_mid] = [[4] for _ in idx_dmn_mid] #DMN midline
    
    idx_dmn_pl = [4,18,19,22] 
    idx_network[idx_dmn_pl] = [[5] for _ in idx_dmn_pl] #DMN posterolateral
    
    idx_sal = [12,13,14]
    idx_network[idx_sal] = [[6] for _ in idx_sal] #SAL
    
    idx_dmn_sub = [29]
    idx_network[idx_dmn_sub] = [[7] for _ in idx_dmn_sub] #DMN SUB

    idx_lcn_sub = [30]
    idx_network[idx_lcn_sub] = [[8] for _ in idx_lcn_sub] #LCN SUB

    idx_hippo = [20,21,24,25,26,27]
    idx_network[idx_hippo] = [[9] for _ in idx_hippo] #hippo
    
    idx_bf = [31,32,33]
    idx_network[idx_bf] = [[10] for _ in idx_bf] #basal forebrain
    
    idx_th = [34,35,36]
    idx_network[idx_th] = [[11] for _ in idx_th] #thal/hypo

    idx_th = [23]
    idx_network[idx_th] = [[-1] for _ in idx_th] #PIR

    idx_th = [28]
    idx_network[idx_th] = [[-1] for _ in idx_th] #CTXsp
    
    idx_network = np.vstack([idx_network, idx_network])
    label = ['LCN', 'DMNmid', 'DMNpl', 'SAL', 'DMNsub', 'LCNsub', 'HPF', 'basal for', 'thal/hy'] #['LCN', 'DMNmidline', 'DMNpostl', 'SAL', 'DMNsub', 'LCNsub', 'HPF', 'basal for', 'thal/hy']
    return np.squeeze(idx_network), label

def get_idx_network_iso_small(): #only LSC, DMN midline, posterolateral and SAL : negative idx if to discard
    
    idx_network = np.zeros([23,1], dtype=int)
    
    idx_lcn = [0,1,2,3]
    idx_network[idx_lcn] = [[1] for _ in idx_lcn] #LCN
    
    idx_vis = [5]
    idx_network[idx_vis] = [[-1] for _ in idx_vis] #VIS

    idx_aud = [6]
    idx_network[idx_aud] = [[-1] for _ in idx_aud] #AUD

    idx_dmn_mid = [7,8,9,10,11,15,16,17]
    idx_network[idx_dmn_mid] = [[4] for _ in idx_dmn_mid] #DMN midline
    
    idx_dmn_pl = [4,18,19,22] 
    idx_network[idx_dmn_pl] = [[5] for _ in idx_dmn_pl] #DMN posterolateral
    
    idx_sal = [12,13,14]
    idx_network[idx_sal] = [[6] for _ in idx_sal] #SAL

    idx_hippo = [20,21]
    idx_network[idx_hippo] = [[-1] for _ in idx_hippo] #hippo
    
    idx_network = np.vstack([idx_network, idx_network])
    label = ['LCN', 'DMNmid', 'DMNpl', 'SAL'] #['LCN', 'DMNmidline', 'DMNpostl', 'SAL']
    return np.squeeze(idx_network), label

def get_idx_network_iso():
    
    idx_network = np.zeros([23,1], dtype=int)
    
    idx_lcn = [0,1,2,3]
    idx_network[idx_lcn] = [[1] for _ in idx_lcn] #LCN
    
    idx_vis = [5]
    idx_network[idx_vis] = [[2] for _ in idx_vis] #VIS

    idx_aud = [6]
    idx_network[idx_aud] = [[3] for _ in idx_aud] #AUD

    idx_dmn_mid = [7,8,9,10,11,15,16,17]
    idx_network[idx_dmn_mid] = [[4] for _ in idx_dmn_mid] #DMN midline
    
    idx_dmn_pl = [4,18,19,22] 
    idx_network[idx_dmn_pl] = [[5] for _ in idx_dmn_pl] #DMN posterolateral
    
    idx_sal = [12,13,14]
    idx_network[idx_sal] = [[6] for _ in idx_sal] #SAL

    idx_hippo = [20,21]
    idx_network[idx_hippo] = [[7] for _ in idx_hippo] #hippo
    
    idx_network = np.vstack([idx_network, idx_network])
    label = ['LCN', 'VIS', 'AUD', 'DMNmidline', 'DMNpostl', 'SAL', 'HPF']
    return np.squeeze(idx_network), label

def get_idx_network():
    
    idx_network = np.zeros([37,1], dtype=int)
    
    idx_lcn = [0,1,2,3,30]
    idx_network[idx_lcn] = [[1] for _ in idx_lcn] #LCN + sub
    
    idx_vis = [5]
    idx_network[idx_vis] = [[2] for _ in idx_vis] #VIS

    idx_aud = [6]
    idx_network[idx_aud] = [[3] for _ in idx_aud] #AUD

    idx_dmn_mid = [7,8,9,10,11,15,16,17]
    idx_network[idx_dmn_mid] = [[4] for _ in idx_dmn_mid] #DMN midline
    
    idx_dmn_pl = [4,18,19,22,29] 
    idx_network[idx_dmn_pl] = [[5] for _ in idx_dmn_pl] #DMN posterolateral + sub
    
    idx_sal = [12,13,14]
    idx_network[idx_sal] = [[6] for _ in idx_sal] #SAL
    
#    idx_dmn_sub = [29]
#    idx_network[idx_dmn_sub] = [[7] for _ in idx_dmn_sub] #DMN SUB

#    idx_lcn_sub = [30]
#    idx_network[idx_lcn_sub] = [[8] for _ in idx_lcn_sub] #LCN SUB

    idx_hippo = [20,21,24,25,26,27]
    idx_network[idx_hippo] = [[7] for _ in idx_hippo] #hippo
    
    idx_bf = [31,32,33]
    idx_network[idx_bf] = [[8] for _ in idx_bf] #basal forebrain
    
    idx_th = [34,35,36]
    idx_network[idx_th] = [[9] for _ in idx_th] #thal/hypo
    
    idx_network = np.vstack([idx_network, idx_network])
    label = ['unknown', 'LCN', 'VIS', 'AUD', 'DMNmidline', 'DMNpostl', 'SAL', 'HPF', 'basal for', 'thal/hy']
    return np.squeeze(idx_network), label
    
def get_distance_matrix(V, metric='euclidean'):
    #for each sbj compute the distance matrix
    
    [n_roi, n_feat, n_sbj] = V.shape
    D = np.zeros([n_roi, n_roi, n_sbj])
    
    for i_sbj in range(n_sbj):
        D[:,:,i_sbj] = metrics.pairwise_distances(V[:,:,i_sbj], metric=metric)
    
    return D

def check_label_order(label, n_cl=2):
    
    freq, _ = np.histogram(label, bins=n_cl, range=[0, n_cl])
    new_label = np.ones(label.shape, dtype=int)
    idx_sort = np.argsort(freq)
    for i, label_i in enumerate(idx_sort):
        new_label[label == label_i] = int(n_cl-i-1)
    
    #freq, _ = np.histogram(new_label, bins=n_cl, range=[0, n_cl])
    #assert(np.unique(np.diff(np.argsort(freq)))==-1)
    return new_label