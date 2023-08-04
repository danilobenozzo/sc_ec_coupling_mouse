#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 19:21:51 2023

@author: benozzo
"""
import scipy.stats as stats
import numpy as np
from function_library import spearman_rank_permutation

def symm_level(ec, sparsified=False):
    #compute the symmetry level by the spear corr row/col
    #sparsified=True only where ec links are non-zero
    #import scipy.stats as stats
    n_roi = ec.shape[0]
    n_test_roi = np.zeros(n_roi)
    p = np.zeros(n_roi)
    s = np.zeros(n_roi)    
    #pdb.set_trace()
    for i_roi in range(n_roi):
        ec_i = ec[i_roi]
        ec_j = ec[:, i_roi]
        idx = np.ones(n_roi, dtype=bool)
        if sparsified:
            idx = np.logical_and(np.abs(ec_j) > 0, np.abs(ec_i) > 0)
        [s[i_roi], p[i_roi]] = stats.spearmanr(np.abs(ec_j[idx]), np.abs(ec_i[idx]))
        n_test_roi[i_roi] = np.sum(idx)
    
    return s, p, n_test_roi
  
def compute_within_ratio(i_roi, pos_roi):
    #compute the ratio between #roi in pos_roi that belong to the same network of i_roi / tot number of rois in that nertork
    from library_networks import get_idx_network_subSeparated
    idx_network, label = get_idx_network_subSeparated()
    return np.sum(idx_network[pos_roi] == idx_network[i_roi]) / float(np.sum(idx_network == idx_network[i_roi]))

def compute_in_network_ratio_scec(i_roi, pos_roi_sc, pos_roi_ec):
    
    from library_networks import get_idx_network_subSeparated
    idx_network, label = get_idx_network_subSeparated()
    
    pos_overlap = []
    pos_overlap += [pos_i for pos_i in pos_roi_sc if pos_i in pos_roi_ec]
    pos_out = []
    pos_out += [pos_i for pos_i in pos_roi_sc if not(pos_i in pos_roi_ec)] #roi sc not included in ec
    #pos_out += [pos_i for pos_i in pos_roi_ec if not(pos_i in pos_roi_sc)] #roi ec not included in sc
    
    #pdb.set_trace()
    #compute the ratio between sc link within network links not overlapping with ec, and sc link within network overlapping
    return np.sum(idx_network[pos_out] == idx_network[i_roi]) / float(np.sum(idx_network[pos_overlap] == idx_network[i_roi]) + np.sum(idx_network[pos_out] == idx_network[i_roi])) 

def compute_out_network_ratio_scec(i_roi, pos_roi_sc, pos_roi_ec):
    
    from library_networks import get_idx_network_subSeparated
    idx_network, label = get_idx_network_subSeparated()
    
    pos_overlap = []
    pos_overlap += [pos_i for pos_i in pos_roi_sc if pos_i in pos_roi_ec]
    pos_out = []
    pos_out += [pos_i for pos_i in pos_roi_sc if not(pos_i in pos_roi_ec)] #roi sc not included in ec
    #pos_out += [pos_i for pos_i in pos_roi_ec if not(pos_i in pos_roi_sc)] #roi ec not included in sc
    
    #pdb.set_trace()
    #compute the ratio between sc link across network links not overlapping with ec, and sc link across network overlapping
    return np.sum(idx_network[pos_out] != idx_network[i_roi]) / float(np.sum(idx_network[pos_overlap] != idx_network[i_roi]) + np.sum(idx_network[pos_out] != idx_network[i_roi]))     
    
def structure_functional_coupling(sc, ec, incoming=True, sparsified=False, n_perm=500, n_roi_min=10):
    #compute the strucutre functional coupling as in baun2020development
    #incoming=True row-wise correlation
    #sparsified=True only where both sc and ec links are non-zero
    #import scipy.stats as stats
    n_roi = sc.shape[0]
    n_test_roi = np.zeros(n_roi)
    p = np.zeros(n_roi)
    s = np.zeros(n_roi)
    if not incoming:
        sc = sc.T
        ec = ec.T
    
    #pdb.set_trace()
    for i_roi in range(n_roi):
        sc_i = sc[i_roi]
        ec_i = ec[i_roi]
        idx = np.ones(n_roi, dtype=bool)
        if sparsified:
            idx = np.logical_and(np.logical_not(np.isinf(sc_i)), np.abs(ec_i) > 0)
        #pdb.set_trace()
        [s[i_roi], p[i_roi]] = spearman_rank_permutation(sc_i[idx], np.abs(ec_i[idx]), n_perm=2) #stats.spearmanr(sc_i[idx], np.abs(ec_i[idx]))
        n_test_roi[i_roi] = np.sum(idx)
    
    n_min = n_roi_min #np.int(np.min(n_test_roi))
    s_small = np.zeros(n_roi)
    p_small = np.zeros(n_roi)
    count_overlap = np.zeros(n_roi)
    within_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] and the size of the network to which i_roi belongs to 
    in_network_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] (with max SC/EC) that do not overlap and are in SC (or EC) within network
    out_network_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] (with max SC/EC) that do not overlap and are in SC (or EC) across network
    for i_roi in range(n_roi):
        sc_i = sc[i_roi]
        ec_i = ec[i_roi]
        idx = np.ones(n_roi, dtype=bool)
        if sparsified:
            idx = np.logical_and(np.logical_not(np.isinf(sc_i)), np.abs(ec_i) > 0)
        #pdb.set_trace()
        idx_pos = np.where(idx)[0]
        idx_pos_sort = np.argsort(sc_i[idx_pos])
        idx_pos_sort_ec = np.argsort(np.abs(ec_i[idx_pos])) #check overlap between sc and abs(ec)
        count_overlap_i = []
        count_overlap_i += [np.sum(idx_pos_sort_ec[-n_min:] == idx_i) for idx_i in idx_pos_sort[-n_min:]]
        count_overlap[i_roi] = np.sum(count_overlap_i)
                
        in_network_ratio[i_roi] = compute_in_network_ratio_scec(i_roi, idx_pos[idx_pos_sort][-n_min:], idx_pos[idx_pos_sort_ec][-n_min:])
        out_network_ratio[i_roi] = compute_out_network_ratio_scec(i_roi, idx_pos[idx_pos_sort][-n_min:], idx_pos[idx_pos_sort_ec][-n_min:])
        within_ratio[i_roi] = compute_within_ratio(i_roi, idx_pos[idx_pos_sort][-n_min:]) #change idx_pos_sort or idx_pos_sort_ec
        ############## end check
        sc_i = sc_i[idx_pos[idx_pos_sort]][-n_min:] #change idx_pos_sort or idx_pos_sort_ec
        ec_i = ec_i[idx_pos[idx_pos_sort]][-n_min:] #change idx_pos_sort or idx_pos_sort_ec
        [s_small[i_roi], p_small[i_roi]] = spearman_rank_permutation(sc_i, np.abs(ec_i), n_perm=n_perm) #stats.spearmanr(sc_i, np.abs(ec_i))        
    
    #pdb.set_trace()
    return s, p, n_test_roi, s_small, p_small, in_network_ratio, out_network_ratio

def compute_in_network_ratio_ecsc(i_roi, pos_roi_sc, pos_roi_ec):
    
    from library_networks import get_idx_network_subSeparated
    idx_network, label = get_idx_network_subSeparated()
    
    pos_overlap = []
    pos_overlap += [pos_i for pos_i in pos_roi_sc if pos_i in pos_roi_ec]
    pos_out = []
    pos_out += [pos_i for pos_i in pos_roi_ec if not(pos_i in pos_roi_sc)] #roi ec not included in sc
    
    #pdb.set_trace()
    #compute the ratio between ec link within network links not overlapping with sc, and ec link within network overlapping
    return np.sum(idx_network[pos_out] == idx_network[i_roi]) / float(np.sum(idx_network[pos_overlap] == idx_network[i_roi]) + np.sum(idx_network[pos_out] == idx_network[i_roi])) 

def compute_out_network_ratio_ecsc(i_roi, pos_roi_sc, pos_roi_ec):
    
    from library_networks import get_idx_network_subSeparated
    idx_network, label = get_idx_network_subSeparated()
    
    pos_overlap = []
    pos_overlap += [pos_i for pos_i in pos_roi_sc if pos_i in pos_roi_ec]
    pos_out = []
    pos_out += [pos_i for pos_i in pos_roi_ec if not(pos_i in pos_roi_sc)] #roi ec not included in sc
    
    #pdb.set_trace()
    #compute the ratio between ec link across network links not overlapping with sc, and ec link across network overlapping
    return np.sum(idx_network[pos_out] != idx_network[i_roi]) / float(np.sum(idx_network[pos_overlap] != idx_network[i_roi]) + np.sum(idx_network[pos_out] != idx_network[i_roi]))     
    
def functional_structure_coupling(sc, ec, incoming=True, sparsified=False, n_perm=500, n_roi_min=10):
    #compute the strucutre functional coupling as in baun2020development
    #incoming=True row-wise correlation
    #sparsified=True only where both sc and ec links are non-zero
    #import scipy.stats as stats
    n_roi = sc.shape[0]
    n_test_roi = np.zeros(n_roi)
    p = np.zeros(n_roi)
    s = np.zeros(n_roi)
    if not incoming:
        sc = sc.T
        ec = ec.T
    
    #pdb.set_trace()
    for i_roi in range(n_roi):
        sc_i = sc[i_roi]
        ec_i = ec[i_roi]
        idx = np.ones(n_roi, dtype=bool)
        if sparsified:
            idx = np.logical_and(np.logical_not(np.isinf(sc_i)), np.abs(ec_i) > 0)
        #pdb.set_trace()
        [s[i_roi], p[i_roi]] = spearman_rank_permutation(sc_i[idx], np.abs(ec_i[idx]), n_perm=2) #stats.spearmanr(sc_i[idx], np.abs(ec_i[idx]))
        n_test_roi[i_roi] = np.sum(idx)
    
    n_min = n_roi_min #np.int(np.min(n_test_roi))
    s_small = np.zeros(n_roi)
    p_small = np.zeros(n_roi)
    count_overlap = np.zeros(n_roi)
    within_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] and the size of the network to which i_roi belongs to 
    in_network_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] (with max SC/EC) that do not overlap and are in SC (or EC) within network
    out_network_ratio = np.zeros(n_roi) #count the ratio between rois in [-n_min:] (with max SC/EC) that do not overlap and are in SC (or EC) across network
    for i_roi in range(n_roi):
        sc_i = sc[i_roi]
        ec_i = ec[i_roi]
        idx = np.ones(n_roi, dtype=bool)
        if sparsified:
            idx = np.logical_and(np.logical_not(np.isinf(sc_i)), np.abs(ec_i) > 0)
        #pdb.set_trace()
        idx_pos = np.where(idx)[0]
        idx_pos_sort = np.argsort(sc_i[idx_pos])
        idx_pos_sort_ec = np.argsort(np.abs(ec_i[idx_pos])) #check overlap between sc and abs(ec)
        count_overlap_i = []
        count_overlap_i += [np.sum(idx_pos_sort_ec[-n_min:] == idx_i) for idx_i in idx_pos_sort[-n_min:]]
        count_overlap[i_roi] = np.sum(count_overlap_i)
                
        in_network_ratio[i_roi] = compute_in_network_ratio_ecsc(i_roi, idx_pos[idx_pos_sort][-n_min:], idx_pos[idx_pos_sort_ec][-n_min:])
        out_network_ratio[i_roi] = compute_out_network_ratio_ecsc(i_roi, idx_pos[idx_pos_sort][-n_min:], idx_pos[idx_pos_sort_ec][-n_min:])
        within_ratio[i_roi] = compute_within_ratio(i_roi, idx_pos[idx_pos_sort_ec][-n_min:]) #change idx_pos_sort or idx_pos_sort_ec
        ############## end check
        sc_i = sc_i[idx_pos[idx_pos_sort_ec]][-n_min:] #change idx_pos_sort or idx_pos_sort_ec
        ec_i = ec_i[idx_pos[idx_pos_sort_ec]][-n_min:] #change idx_pos_sort or idx_pos_sort_ec
        [s_small[i_roi], p_small[i_roi]] = spearman_rank_permutation(sc_i, np.abs(ec_i), n_perm=n_perm) #stats.spearmanr(sc_i, np.abs(ec_i))        
    
    #pdb.set_trace()
    return s, p, n_test_roi, s_small, p_small, in_network_ratio, out_network_ratio