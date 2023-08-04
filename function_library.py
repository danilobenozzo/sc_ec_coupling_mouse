#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 16:39:48 2021

@author: benozzo
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import stats
from joblib import Parallel, delayed
import pdb

def compute_stat_test(data, group):
    "one way anova and others simple test"
    import pandas
    from statsmodels.stats.libqsturng import psturng
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats import multicomp
    
    #hval, pval_test = stats.kruskal(data[0], data[1], data[2])
    #print "pval kruskal test %f" % (pval_test)
    #hval, pval_test = stats.f_oneway(data[0], data[1], data[2])
    #print "pval anova test %f" % (pval_test)
        
    data = np.hstack(data) #np.vstack(data)
    d = {'data':data, 'group':group}
    df = pandas.DataFrame(data=d)
    model = ols('data ~ C(group)', df).fit() #C(go)*C(corr)
    model.summary()
    #Create the ANOVA table
    res = sm.stats.anova_lm(model, typ=1)
    print res 
    
    #Post hoc testing
    mc = multicomp.MultiComparison(df['data'], df['group'])
    mc_results = mc.tukeyhsd()
    print(mc_results)
    
    p_values =  psturng(np.abs(mc_results.meandiffs / mc_results.std_pairs), len(mc_results.groupsunique), mc_results.df_total)
    print 'p value: ', p_values
    
def get_idx_network():

    idx_network = np.array([7, 8, 9, 10, 11, 15, 16, 17, 44, 45, 46, 47, 48, 52, 53, 54, 4, 18, 19, 22, 41, 55, 56, 59, 29, 66,
                              0, 1, 2, 3, 30, 37, 38, 39, 40, 67,
                              12, 13, 14, 49, 50, 51,
                              5, 42,
                              6, 43,
                              20, 21, 24, 25, 26, 27, 57, 58, 61, 62, 63, 64,
                              23, 60,
                              28, 65,
                              31, 32, 33, 68, 69, 70,
                              34, 35, 71, 72,
                              36, 73], dtype=int) #DMNm, DMNpl, DMNsub, LCN, SAL , VIS, AUD, HPC, PIR, CTXsp, BF, THAL, HY [second attempt functinal network sept 24th]

    separation = np.array([8, 12, 13, 18, 21, 22, 23, 29, 30, 31, 34, 36, 37], dtype=float)*2 -1
    xtick_label = ['DMNmid', 'DMNpl', 'DMNsub', 'LCN', 'SAL', 'VIS', 'AUD', 'HPC', 'PIR', 'CTXsp', 'BF', 'THAL', 'HY']
    
    return idx_network, separation, xtick_label

def func_net_ordering(mtx):

    idx_network, separation, xtick_label = get_idx_network()
    mtx = mtx[idx_network, :][:, idx_network]
    return mtx, separation, xtick_label


def reducing_toNetwork(mtx):
    #reduce full matrix mtx IN THE HEMISPHERE ORDER to a network based smaller matrix
    
    mtx, separation, xtick_label = func_net_ordering(mtx)
    
    n = len(separation)
    reduced_mtx = np.zeros([n, n])
    separation = np.hstack([0, np.array(separation+1, dtype=int)])
    for i, idx_i in enumerate(separation[:-1]):
        for j, idx_j in enumerate(separation[:-1]):
            #print '(%d, %d) - (%d, %d)' % (idx_i, separation[i+1], idx_j, separation[j+1])
            tmp_mtx = mtx[idx_i : separation[i+1], idx_j : separation[j+1]]
            if i == j:
                #tmp_mtx is square => remove diag
                tmp_mtx = tmp_mtx[np.logical_not(np.identity(len(tmp_mtx)))]
            reduced_mtx[i, j] = np.mean(tmp_mtx) #np.mean(np.sum(tmp_mtx, 1)) #TO SET
    return reduced_mtx, xtick_label

#try to write a new version    
def get_asym_mtx(mtx):
    #compute asymetry matrix by the difference between upper and  lower triang part
    #interpretation: global exchange if interaction between (i,j)
    n = mtx.shape[0]      
    mtx_asym = np.zeros([n, n])
    idx_x, idx_y = np.triu_indices(n, 1)
    tmp = mtx[idx_x, idx_y] - mtx[idx_y, idx_x]
    mtx_asym[idx_y, idx_x] = tmp
       
    return mtx_asym, tmp

#################################################
#def get_asym_mtx(mtx, hem_sep=True):
#    #compute asymetry matrix by the difference between upper and  lower triang part
#    #interpretation: global exchange if interaction between (i,j)
#    n = mtx.shape[0]
#    if not hem_sep:
#        idx_network, separation, xtick_label = get_idx_network()
#        mtx = mtx[idx_network, :][:, idx_network]
#        
#    mtx_asym = np.zeros([n, n])
#    idx_x, idx_y = np.triu_indices(n, 1)
#    tmp = mtx[idx_x, idx_y] - mtx[idx_y, idx_x]
#    mtx_asym[idx_y, idx_x] = tmp
#    
#    idx_back = None
#    if not hem_sep:
#        idx_back = np.zeros(idx_network.shape, dtype=int)
#        for i, idx_i in enumerate(idx_network):
#            idx_back[idx_i] = i
#        mtx_asym = mtx_asym[idx_back, :][:, idx_back]
#
#    
#    return mtx_asym, tmp, idx_back

##############################################
def plot_asym_mtx(vet_asym, nS, filename_tosave, pwd_tosave, hem_sep, th_asym=1.5, barcolor=True):
    
    from fig1 import im_matrix, plot_reduced_toNetwork
    #mtx_asym_mean less variable entries across mice
    mtx_asym_mean = np.zeros([nS, nS])
    vet_asym_mean = np.mean(vet_asym, 0)
    idx_y, idx_x = np.triu_indices(nS, 1) #trick to coherently select upper and lower diag
    mtx_asym_mean[idx_x, idx_y] = vet_asym_mean
    #if not hemSep:
    #    mtx_asym_mean = mtx_asym_mean[idx_back, :][:, idx_back]
    if nS == 74:
        im_matrix(mtx_asym_mean, filename_tosave = filename_tosave, pwd_tosave = pwd_tosave, hem_sep=hem_sep, barcolor=barcolor)
    elif nS == 13:
        _, _, xtick_label = get_idx_network()
        plot_reduced_toNetwork(mtx_asym_mean, xtick_label, filename_tosave = filename_tosave, pwd_tosave=pwd_tosave, barcolor=barcolor)
    else:
        print "Error in plot_asym_mtx"
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(vet_asym, aspect='auto')
    plt.subplot(2,1,2)
    plt.plot(np.abs(np.std(vet_asym, 0) / np.mean(vet_asym, 0)))
    plt.xlim([0, vet_asym.shape[1]])
    plt.ylim([0, th_asym])
    plt.show()
    
    #th_asym = 1.5 #np.percentile(np.std(vet_asym, 0), 20)
    vet_asym_mean = np.mean(vet_asym, 0)
    vet_asym_mean[np.abs(np.std(vet_asym, 0) / np.mean(vet_asym, 0)) > th_asym] = 0
    mtx_asym_mean = np.zeros([nS, nS])
    mtx_asym_mean[idx_x, idx_y] = vet_asym_mean
    #if not hemSep:
    #    mtx_asym_mean = mtx_asym_mean[idx_back, :][:, idx_back]
    if nS == 74:
        im_matrix(mtx_asym_mean, filename_tosave = '%s_thCV' % filename_tosave, pwd_tosave = pwd_tosave, hem_sep=hem_sep, barcolor=barcolor)
    elif nS == 13:
        plot_reduced_toNetwork(mtx_asym_mean, xtick_label, filename_tosave = '%s_thCV' % filename_tosave, pwd_tosave=pwd_tosave, barcolor=barcolor)
    else:
        print "Error in plot_asym_mtx"
    
def do_permutation(x, y, seed):
    np.random.seed(seed)
    r, _ = stats.spearmanr(x[np.random.permutation(x.shape[0])], y)
    return r
    
def spearman_rank_permutation(x, y, n_perm=1000):
    
    r, p = stats.spearmanr(x, y)
    #r_perm = np.zeros(n_perm)
    
    #pdb.set_trace()
    seeds = np.arange(n_perm)
    r_perm = Parallel(n_jobs=-1)(delayed(do_permutation)(x, y, i_seed) for i_seed in seeds)
    r_perm = np.array(r_perm, dtype=float)
    
    #for i_perm in range(n_perm):
    #    r_perm[i_perm], _ = stats.spearmanr(x[np.random.permutation(x.shape[0])], y)
    p_perm = (np.sum(np.abs(r_perm) >= np.abs(r))+1) / float(n_perm+1)
    return r, p_perm

def mean_perNetwork_acrossSbj(X, pwd_tosave='/home/benozzo/Desktop/', filename='measure_perNetwork'):
    #X: [n_roi n_sbj]
    from library_networks import get_idx_network as get_network_index, get_idx_network_iso, get_idx_network_iso_small, get_idx_network_subSeparated
    from plotting_functions import figure_layout, set_size
    iso_network = [False, True]
    X_network = []
    
    for is_iso in iso_network:
        if is_iso:
            idx_network, label_network = get_idx_network_iso_small()
            idx_iso = np.hstack([np.arange(23), np.arange(23)+37]) 
            w = 1.8 #3.1
        else:
            idx_network, label_network = get_idx_network_subSeparated() 
            idx_iso = np.arange(74)          
            w = 1.8 #6.2
        
        h = 4 #3.1
        y_in = 2.5#important for figure_layout
        x_in = w / h * y_in 
        fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
        fig, ax = plt.subplots(2, 1, figsize=(x_in, y_in))
    
        #pdb.set_trace()
        n_cl = np.sum(np.unique(idx_network) > -1) #network subdivision
        cm_figure = cm.get_cmap('tab10', n_cl)
        flierprops = dict(marker='+', markerfacecolor='k', markersize=markersize, linestyle='none', markeredgecolor='k')
        medianprops = dict(linestyle='-', linewidth=linewidth, color='k')       
        #fig, ax = plt.subplots(2, 1)
        
        for i_dim, X_i in enumerate(X):
            X_i = X_i[idx_iso, :]
            X_network_i = []                                 
            
            for i_cl, cl_i in enumerate(np.unique(idx_network)): #for each network
                if ~ cl_i < 0:
                    X_network_i += [np.mean(X_i[idx_network == cl_i, :], 0)]# mean across roi for each sbj #[np.hstack(X_i[idx_network == cl_i, :])] 
            
            median_in = []
            median_in += [np.median(Xi) for Xi in X_network_i]
            idx_sort = np.argsort(median_in) #order median
            
            X_network.append(X_network_i)
            label_orded = []
            for i_idx, idx_i in enumerate(idx_sort): #(range(len(idx_sort))): #keep the same order 
                #bp_i = ax[i_dim].bar(i_idx, np.mean(X_network_i[idx_i]) - .5, width=.3, color=cm_figure.colors[idx_i], edgecolor='k', linewidth=linewidth_main, yerr=np.std(X_network_i[idx_i])) 
                bp_i = ax[i_dim].boxplot(X_network_i[idx_i], positions=[i_idx], patch_artist=True, flierprops=flierprops)
                bp_i['boxes'][0].set_facecolor(cm_figure.colors[idx_i])
                label_orded += [label_network[idx_i]]
                #ax[i_dim].set_xticklabels(label_orded[i_idx], rotation=15)
            
            #ax[i_dim].plot([-.5, n_cl-.5], [0, 0], '--k', linewidth=linewidth) #added
            ax[i_dim].set_xticks(range(n_cl))
            ax[i_dim].set_xticklabels(label_orded, rotation='vertical', fontsize=fontsize) #rotation 15
            ax[i_dim].set_xlim([-1, n_cl])
            #ax[i_dim].set_yticks(np.arange(-.5, .51, .25)) #added
            #ax[i_dim].set_yticklabels(np.arange(0, 1.1, .25)) #adde
            ax[i_dim].set_ylim([0, 1.])#([-.0005, .003])#([-2., 2.])
            ax[i_dim].grid(axis='y')
    
        ax[0].set_ylabel('incoming coupling', fontsize=fontsize)
        ax[1].set_ylabel('outgoing coupling', fontsize=fontsize)
        ax[0].tick_params(labelsize=fontsize)
        ax[1].tick_params(labelsize=fontsize)
        ax[0] = set_size(w, h, ax[0]) #one axis is enough
        fig.tight_layout()
        if is_iso:
            fig.savefig('%s%s_iso.svg' % (pwd_tosave, filename), format='svg')
        else:
            fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
        plt.close()
    
    return X_network