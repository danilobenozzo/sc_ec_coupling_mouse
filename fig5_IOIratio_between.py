#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:43:34 2022

@author: benozzo
"""
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import stats
import pickle
from library_networks import get_idx_network as get_network_index, get_idx_network_iso, get_idx_network_iso_small, get_idx_network_subSeparated
from plotting_functions import figure_layout, set_size
import numpy as np 

def plot_ratio_between(pwd_data):
    
    idx_network, label_network = get_idx_network_iso_small()
    idx_iso = np.hstack([np.arange(23), np.arange(23)+37])
    w = 3.1

    label_figure = ['SC', 'EC']
    label_path = ['sc_ec', 'ec_sc']
    
    for i_label, label_i in enumerate(label_path):
        h = 3.1
        y_in = 2.5#important for figure_layout
        x_in = w / h * y_in 
        fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
        fig, ax = plt.subplots(2, 1, figsize=(x_in, y_in))
    
        #pdb.set_trace()
        flierprops = dict(marker='+', markerfacecolor='k', markersize=markersize, linestyle='none', markeredgecolor='k')
        n_cl = np.sum(np.unique(idx_network) > -1) #network subdivision
        cm_figure = cm.get_cmap('tab10', n_cl)
        
        n_sim = [62131, 92131, 102131, 82131]
        alpha = [1, .3, .3, .3]
        delta_x = np.linspace(-.3, .3, 4)
        for i_plot, n_sim_i in enumerate(n_sim):
            
            dataset = pickle.load(open('%s%s/%s/output_onlyRate_across_%s.pickle' % (pwd_data, n_sim_i, label_i, n_sim_i)))
            X_network_i = dataset['X_network'][2:]                            
            
            for i_dim, X_i in enumerate(X_network_i): #for incoming / outoing
                
                median_in = []
                median_in += [np.median(Xi) for Xi in X_i]
                print '----------------'
                print 'sim = ', n_sim_i
                print 'median ratio per netw: ', median_in
                idx_sort = np.argsort(median_in) #order median
                
                label_orded = []
                for i_idx, idx_i in enumerate(range(len(idx_sort))): #keep the same order #enumerate(idx_sort): 
                    bp_i = ax[i_dim].bar(i_idx + delta_x[i_plot], np.mean(X_i[idx_i]) - .5, width=.2, color=cm_figure.colors[idx_i], edgecolor='k', linewidth=linewidth_main, yerr=np.std(X_i[idx_i]), alpha=alpha[i_plot]) 
                    #bp_i = ax[i_dim].boxplot(X_i[idx_i], positions=[i_idx + delta_x[i_plot]], patch_artist=True, flierprops=flierprops)
                    #bp_i['boxes'][0].set_facecolor(cm_figure.colors[idx_i])
                    label_orded += [label_network[idx_i]]
                    #ax[i_dim].set_xticklabels(label_orded[i_idx], rotation=15)
                
                ax[i_dim].plot([-.5, n_cl-.5], [0, 0], '--k', linewidth=linewidth) #added
                ax[i_dim].set_xticks(range(n_cl))
                ax[i_dim].set_xticklabels(label_orded, rotation=15, fontsize=fontsize)
                ax[i_dim].set_xlim([-1, n_cl])
                ax[i_dim].set_yticks(np.arange(-.5, .51, .25)) #added
                ax[i_dim].set_yticklabels(np.arange(0, 1.1, .25)) #adde
                ax[i_dim].set_ylim([-.5, .5]) #added to have equal ylim 8oct22
                #ax[i_dim].set_ylim([0., 1.])
                #ax[i_dim].grid(axis='y')
    
        ax[0].set_ylabel('in only %s ratio' % label_figure[i_label], fontsize=fontsize)
        ax[1].set_ylabel('out only %s ratio' % label_figure[i_label], fontsize=fontsize)
        ax[0].tick_params(labelsize=fontsize)
        ax[1].tick_params(labelsize=fontsize)
        ax[0] = set_size(w, h, ax[0]) #one axis is enough
        
        pwd_tosave = pwd_data
        filename = 'in_out_onlyRatio_across_perNetwork_bar_iso_%s' % label_figure[i_label]
        fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
        plt.close()