# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:23:53 2019

@author: danilo
"""

import matplotlib.pyplot as plt
import numpy as np
from plotting_functions import figure_layout, set_size
  
def table_plot(str_all, filename_tosave='example', pwd_tosave=None, hem_sep=True, mouse_name=None):
##############################
    #from plotting_functions import figure_layout, set_size
    from matplotlib import cm
    x_in = 5.5 #
    y_in = 2. #
    h = 3.1
    w = x_in / (y_in/h) #2.1

    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)
    n_sbjs, n_roi = str_all.shape
    
    fig, ax = plt.subplots(figsize=(x_in, y_in))
    cm = plt.cm.get_cmap('Greens')
    plt.imshow(str_all, aspect='auto', cmap=cm)
    shift = .5
    if hem_sep:
        separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float) -1  
        xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        plt.plot([separation+shift, separation+shift], [(-shift)*np.ones(separation.shape), (n_sbjs+shift)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [(-shift)*np.ones(separation.shape), (n_sbjs+shift)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + n_roi/2]), xtick_label + xtick_label, fontsize=fontsize_main, rotation='vertical')
    else:
        #anatomical aggregation across hemi
        #separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        #functional network first attempt
        #separation = np.array([6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['DMN', 'LCN', 'MOs', 'SSs', 'GU', 'VIS', 'AUD', 'PL', 'ILA', 'ORB', 'AI', 'PTLp', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
       
        #functional network second attempt (mail Sept 24th)
        separation = np.array([8, 12, 13, 18, 21, 22, 23, 29, 30, 31, 34, 36, 37], dtype=float)*2 -1
        xtick_label = ['DMNmid', 'DMNpl', 'DMNsub', 'LCN', 'SAL', 'VIS', 'AUD', 'HPC', 'PIR', 'CTXsp', 'BF', 'THAL', 'HY']
        
        plt.plot([separation+shift, separation+shift], [(-shift)*np.ones(separation.shape), (n_sbjs+shift)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        #plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [np.min(str_all)*np.ones(separation.shape)-shift, np.max(str_all)*np.ones(len(separation))+shift], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(separation-delta_sep+shift, xtick_label, fontsize=fontsize_main)#, rotation='vertical')
        
    #cbar.set_ticks(np.linspace(-2, 2, 5, dtype=int))
    #cbar.ax.set_yticklabels(np.linspace(-2, 2, 5, dtype=int), fontsize=fontsize)    
    
    plt.xlim([-shift, n_roi-shift])
    if mouse_name:
        plt.yticks(range(n_sbjs), mouse_name, fontsize=fontsize_main)
    else:
        plt.yticks(range(n_sbjs))
    plt.ylim([-shift, n_sbjs-shift])    
    
    #ax.tick_params(labelsize=4.5)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    
    if not(pwd_tosave):
        pwd_tosave = '/home/benozzo/Desktop/' 
        filename_tosave = 'figure' 
    #plt.savefig('%s%s.pdf' % (pwd_tosave, filename_tosave), dpi=300)
    #plt.savefig('%s%s.png' % (pwd_tosave, filename_tosave), dpi=300)
    plt.savefig('%s%s.svg' % (pwd_tosave, filename_tosave), format='svg')
    plt.close()
    
def in_out_plot(str_all, filename_tosave='example', pwd_tosave=None, hem_sep=True, bar_plot=False):
##############################
    #from plotting_functions import figure_layout, set_size
    from matplotlib import cm
    x_in = 3.5 #
    y_in = 2. #
    h = 3.1
    w = x_in / (y_in/h) #2.1

    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)
    n_roi = str_all.shape[1]
    
    if not(hem_sep):
        idx_step = np.array([0, 2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=int) 
        idx_noHem = []
        for i, idx_i in enumerate(idx_step[1:]): #idx_step
            idx_noHem = np.hstack([idx_noHem, range(idx_step[i], idx_i)])
            idx_noHem = np.hstack([idx_noHem, range(idx_step[i]+n_roi/2, idx_i+n_roi/2)])
        idx_noHem = np.array(idx_noHem, dtype=int) # some regions of different hem, next each other 
        
        idx_noHem = np.array([ 7,  8, 15, 16, 17, 19, 44, 45, 52, 53, 54, 56, 2, 0, 39, 37,  1, 
                              38,  3, 40,  4, 41,  5, 42,  6, 43,  9, 46, 10, 47, 11, 48, 12, 13,
                              14, 49, 50, 51, 18, 55, 20, 57, 21, 58, 22, 59, 23, 60, 24, 25, 26,
                              27, 61, 62, 63, 64, 28, 65, 29, 30, 31, 32, 66, 67, 68, 69, 33, 70,
                              34, 35, 36, 71, 72, 73], dtype=int) #DMN and LCN first [first attmept functional network]
    
        idx_noHem = np.array([7, 8, 9, 10, 11, 15, 16, 17, 44, 45, 46, 47, 48, 52, 53, 54, 4, 18, 19, 22, 41, 55, 56, 59, 29, 66,
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
        
        #str_all = str_all[:, idx_noHem]
        #idx_noHem = np.array([0:2, 37:39, 2:4, 39:41, 4, 41, 5, 42, 6, 43, 7:9, 44:46, 9, 46, 10, 47, 11, 48, 12:15, 49:52, 15:18, 52:55, 18, 55, 19, 56, 20, 57, 21, 58, 22, 59, 23, 60, 24:28, 61:65, 28, 65, 29:33, 66:70, 33, 70, 34:37, 71:74], dtype=int)

    fig, ax = plt.subplots(figsize=(x_in, y_in))
    colors = [cm.hsv(x) for x in np.linspace(0, 1, str_all.shape[0])]
    for i_sb in range(str_all.shape[0]):
        if bar_plot:
            plt.bar(range(n_roi), str_all[i_sb], color=colors[i_sb])
        else:
            plt.plot(str_all[i_sb], '.-', color=colors[i_sb], linewidth=linewidth)
    shift = .5
    if hem_sep:
        separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float) -1  
        xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        plt.plot([separation+shift, separation+shift], [np.min(str_all)*np.ones(separation.shape), np.max(str_all)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [np.min(str_all)*np.ones(separation.shape), np.max(str_all)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + n_roi/2]), xtick_label + xtick_label, fontsize=fontsize, rotation='vertical')
    else:
        #anatomical aggregation across hemi
        #separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        #functional network first attempt
        #separation = np.array([6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['DMN', 'LCN', 'MOs', 'SSs', 'GU', 'VIS', 'AUD', 'PL', 'ILA', 'ORB', 'AI', 'PTLp', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
       
        #functional network second attempt (mail Sept 24th)
        separation = np.array([8, 12, 13, 18, 21, 22, 23, 29, 30, 31, 34, 36, 37], dtype=float)*2 -1
        xtick_label = ['DMNmid', 'DMNpl', 'DMNsub', 'LCN', 'SAL', 'VIS', 'AUD', 'HPC', 'PIR', 'CTXsp', 'BF', 'THAL', 'HY']
        
        plt.plot([separation+shift, separation+shift], [np.min(str_all)*np.ones(separation.shape), np.max(str_all)*np.ones(len(separation))], '-k', linewidth=linewidth) #vertical line
        #plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [np.min(str_all)*np.ones(separation.shape)-shift, np.max(str_all)*np.ones(len(separation))+shift], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(separation-delta_sep+shift, xtick_label, fontsize=fontsize, rotation='vertical')
        
    #cbar.set_ticks(np.linspace(-2, 2, 5, dtype=int))
    #cbar.ax.set_yticklabels(np.linspace(-2, 2, 5, dtype=int), fontsize=fontsize)
    ax.tick_params(labelsize=4.5)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    
    plt.xlim([-shift, n_roi-shift])
    #plt.ylim([n_roi-shift, -shift])    
    if not(pwd_tosave):
        pwd_tosave = '/home/benozzo/Desktop/' 
        filename_tosave = 'figure' 
    #plt.savefig('%s%s.pdf' % (pwd_tosave, filename_tosave), dpi=300)
    #plt.savefig('%s%s.png' % (pwd_tosave, filename_tosave), dpi=300)
    plt.savefig('%s%s.svg' % (pwd_tosave, filename_tosave), format='svg')
    plt.close()
    
    
def im_matrix(im_mtx, filename_tosave='example', pwd_tosave=None, vlim=None, barcolor=False, cm=None, hem_sep=True):
##############################
    #from plotting_functions import figure_layout, set_size
    
    x_in = 2.5 #
    y_in = 2.5 #
    h = 3.1
    w = x_in / (y_in/h) #2.1

    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)
    n_roi = im_mtx.shape[0]
        
    fig, ax = plt.subplots(figsize=(x_in, y_in))
    #colors = [cm.seismic(x) for x in np.linspace(0, 1, 100)]
    if not(cm):
        cm = plt.cm.get_cmap('seismic')#('YlGn')#
    if not(vlim):
        ind_max = np.unravel_index(np.argmax(np.abs(im_mtx), axis=None), im_mtx.shape)
        vlim = [np.abs(im_mtx[ind_max])*(-1), np.abs(im_mtx[ind_max])]
    
    if not(hem_sep):
        idx_step = np.array([0, 2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=int) 
        idx_noHem = []
        for i, idx_i in enumerate(idx_step[1:]): #idx_step
            idx_noHem = np.hstack([idx_noHem, range(idx_step[i], idx_i)])
            idx_noHem = np.hstack([idx_noHem, range(idx_step[i]+n_roi/2, idx_i+n_roi/2)])
        idx_noHem = np.array(idx_noHem, dtype=int) # some regions of different hem, next each other 
        
        idx_noHem = np.array([ 7,  8, 15, 16, 17, 19, 44, 45, 52, 53, 54, 56, 2, 0, 39, 37,  1, 
                              38,  3, 40,  4, 41,  5, 42,  6, 43,  9, 46, 10, 47, 11, 48, 12, 13,
                              14, 49, 50, 51, 18, 55, 20, 57, 21, 58, 22, 59, 23, 60, 24, 25, 26,
                              27, 61, 62, 63, 64, 28, 65, 29, 30, 31, 32, 66, 67, 68, 69, 33, 70,
                              34, 35, 36, 71, 72, 73], dtype=int) #DMN and LCN first [first attmept functional network]
    
        idx_noHem = np.array([7, 8, 9, 10, 11, 15, 16, 17, 44, 45, 46, 47, 48, 52, 53, 54, 4, 18, 19, 22, 41, 55, 56, 59, 29, 66,
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
        
        im_mtx = im_mtx[idx_noHem, :][:, idx_noHem]
        #idx_noHem = np.array([0:2, 37:39, 2:4, 39:41, 4, 41, 5, 42, 6, 43, 7:9, 44:46, 9, 46, 10, 47, 11, 48, 12:15, 49:52, 15:18, 52:55, 18, 55, 19, 56, 20, 57, 21, 58, 22, 59, 23, 60, 24:28, 61:65, 28, 65, 29:33, 66:70, 33, 70, 34:37, 71:74], dtype=int)
    im = plt.imshow(im_mtx, cmap=cm, vmin=vlim[0], vmax=vlim[1])
    
    if barcolor:
        #cax = fig.add_axes([0.1, .95, 0.90, 0.03])
        #cbar = fig.colorbar(im, orientation='horizontal', cax=cax)
        cbar = fig.colorbar(im, fraction=0.05, pad=0.04, shrink=.7, orientation='vertical')
        cbar.ax.tick_params(labelsize=fontsize) 
    
    shift = .5
    if hem_sep:
        separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float) -1  
        xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        plt.plot([separation+shift, separation+shift], [np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], '-k', linewidth=linewidth) #vertical line
        plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + n_roi/2]), xtick_label + xtick_label, fontsize=fontsize, rotation='vertical')
        
        plt.plot([np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], [separation+shift, separation+shift], '-k', linewidth=linewidth)
        plt.plot([np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], [separation+shift+n_roi/2, separation+shift+n_roi/2], '-k', linewidth=linewidth)
        plt.yticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + n_roi/2]), xtick_label + xtick_label, fontsize=fontsize)
    else:
        #anatomical aggregation across hemi
        #separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
        
        #functional network first attempt
        #separation = np.array([6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float)*2 -1  
        #xtick_label = ['DMN', 'LCN', 'MOs', 'SSs', 'GU', 'VIS', 'AUD', 'PL', 'ILA', 'ORB', 'AI', 'PTLp', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
       
        #functional network second attempt (mail Sept 24th)
        separation = np.array([8, 12, 13, 18, 21, 22, 23, 29, 30, 31, 34, 36, 37], dtype=float)*2 -1
        xtick_label = ['DMNmid', 'DMNpl', 'DMNsub', 'LCN', 'SAL', 'VIS', 'AUD', 'HPC', 'PIR', 'CTXsp', 'BF', 'THAL', 'HY']
        
        plt.plot([separation+shift, separation+shift], [np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], '-k', linewidth=linewidth) #vertical line
        #plt.plot([separation+shift+n_roi/2, separation+shift+n_roi/2], [np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], '-k', linewidth=linewidth) #vertical line
        delta_sep = np.diff(separation, prepend=-1) / 2
        plt.xticks(separation-delta_sep+shift, xtick_label, fontsize=fontsize, rotation='vertical')
        
        plt.plot([np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], [separation+shift, separation+shift], '-k', linewidth=linewidth)
        #plt.plot([np.zeros(separation.shape)-shift, np.ones(len(separation))*n_roi-shift], [separation+shift+n_roi/2, separation+shift+n_roi/2], '-k', linewidth=linewidth)
        plt.yticks(separation-delta_sep+shift, xtick_label, fontsize=fontsize)
    #cbar.set_ticks(np.linspace(-2, 2, 5, dtype=int))
    #cbar.ax.set_yticklabels(np.linspace(-2, 2, 5, dtype=int), fontsize=fontsize)
    ax.tick_params(labelsize=4.5)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    
    plt.xlim([-shift, n_roi-shift])
    plt.ylim([n_roi-shift, -shift])
    
    if not(pwd_tosave):
        pwd_tosave = '/home/benozzo/Desktop/' 
        filename_tosave = 'figure' 
    #plt.savefig('%s%s.pdf' % (pwd_tosave, filename_tosave), dpi=300)
    #plt.savefig('%s%s.png' % (pwd_tosave, filename_tosave), dpi=300)
    plt.savefig('%s%s.svg' % (pwd_tosave, filename_tosave), format='svg')
    plt.show()
    plt.close()
    





    