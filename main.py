# -*- coding: utf-8 -*-
"""
"""
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from statsmodels.stats import multitest as mtest 
import pickle
import numpy as np
import pdb
import os

from plotting_functions import figure_layout, set_size
from library_plot import im_matrix, in_out_plot
from function_library import func_net_ordering, mean_perNetwork_acrossSbj

from function_library_coupling import structure_functional_coupling, functional_structure_coupling, symm_level
        
def load_structural_filter(pwd_data, th_prc=None):
    #load structural data and filter at the th_prc percentile of the log transformed mtx
    
    dataset = loadmat('%sasymm_ncd_no_thr_N_74.mat' % pwd_data, squeeze_me=True)
    mtx_str = dataset['full_connectome_no_symm']
    
    idx_ = ~np.array(np.eye(mtx_str.shape[0]), dtype=bool) #idx out diag
    t = np.log(mtx_str[idx_])
    
    mtx_log = np.copy(mtx_str)
    mtx_log[idx_] = np.log(mtx_str[idx_])
    
    if th_prc:
        mtx_th_log = np.copy(mtx_log)
        mtx_th_log[mtx_log < np.percentile(t, th_prc)] = -np.Inf
        return mtx_log, mtx_th_log
    else:
        return mtx_log, mtx_str

def plot_save_results(s_in, s_out, onlySCrate_within_in, onlySCrate_within_out, onlySCrate_across_in, onlySCrate_across_out, pwd_tosave, permute_corr, n_sim):
    
    #plot matrix
    in_coupling = np.dot(s_in.T, s_in)
    im_matrix(in_coupling, filename_tosave='incoming_coupling', pwd_tosave=pwd_tosave, barcolor=True)
    in_coupling_ntw, _, _ = func_net_ordering(in_coupling)
    im_matrix(in_coupling_ntw, filename_tosave='incoming_coupling_ntw', pwd_tosave=pwd_tosave, barcolor=True, hem_sep=False)
    out_coupling = np.dot(s_out.T, s_out)
    im_matrix(out_coupling, filename_tosave='outgoing_coupling', pwd_tosave=pwd_tosave, barcolor=True)
    out_coupling_ntw, _, _ = func_net_ordering(out_coupling)
    im_matrix(out_coupling_ntw, filename_tosave='outgoing_coupling_ntw', pwd_tosave=pwd_tosave, barcolor=True, hem_sep=False)
    
    #plot curve
    in_out_plot(s_in, filename_tosave='incoming_couplig_curve', pwd_tosave=pwd_tosave)
    in_out_plot(s_out, filename_tosave='outgoing_coupling_curve', pwd_tosave=pwd_tosave)
    
    
    s_in_mean = np.zeros(n_sbj)
    s_out_mean = np.zeros(n_sbj)
    for i_sbj in range(n_sbj):
        s_in_mean[i_sbj] = np.mean(s_in[i_sbj]) 
        s_out_mean[i_sbj] = np.mean(s_out[i_sbj]) 
    [hin, x_axis_in] = np.histogram(s_in_mean, bins=6, density=True)
    [hout, x_axis_out] = np.histogram(s_out_mean, bins=6, density=True)

    w = 3.1
    h = 2.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
    fig, ax = plt.subplots(figsize=(x_in, y_in))

    plt.plot(x_axis_in[:-1], hin, '-r', label='incoming')
    plt.plot(x_axis_out[:-1], hout, '-b', label='outgoing')
    plt.plot(np.mean(s_in_mean), 0, 'sr')
    plt.plot(np.mean(s_out_mean), 0, 'sb')
    plt.legend(fontsize=fontsize)
    [_, ptest] = stats.ttest_rel(s_in_mean, s_out_mean)
    plt.title('p val: %.4f' % ptest, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    plt.savefig('%sdist_in_out.svg' % (pwd_tosave), format='svg')
    plt.close()
     
    #test single roi
    nS = s_in.shape[1]
    p_rois = np.zeros(nS)
    for i_roi in range(nS):
        _, p_rois[i_roi] = stats.ttest_rel(s_in[:, i_roi], s_out[:, i_roi])
    
    reject, p_rois_corrected, _, _ = mtest.multipletests(p_rois, alpha=.05, method='fdr_bh') #Benjamini/Hochberg 
    
    w = 6.2
    h = 2.1 #3.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
    fig, ax = plt.subplots(figsize=(x_in, y_in))
        
    plt.plot(np.mean(s_in, 0), '.-r', label='incoming', linewidth=linewidth_main)
    plt.plot(np.vstack([range(nS), range(nS)]), np.vstack([np.mean(s_in, 0)+np.std(s_in, 0), np.mean(s_in, 0)-np.std(s_in, 0)]), '-', color='r', linewidth=linewidth)
    plt.plot(np.mean(s_out, 0), '.-b', label='outgoing', linewidth=linewidth_main)
    plt.plot(np.vstack([range(nS), range(nS)]), np.vstack([np.mean(s_out, 0)+np.std(s_out, 0), np.mean(s_out, 0)-np.std(s_out, 0)]), '-', color='b', linewidth=linewidth)
    plt.legend(fontsize=fontsize)
    plt.plot(np.arange(nS)[reject], np.ones(nS)[reject]*0.02, '*k', markersize=markersize)
    plt.title('mean z-SpearCorr in/out')
    
    shift = .5
    separation = np.array([2, 4, 5, 6, 7, 9, 10, 11, 12, 15, 18, 19, 20, 21, 22, 23, 24, 28, 29, 33, 34, 37], dtype=float) -1  
    xtick_label = ['MO', 'SS', 'GU', 'VIS', 'AUD', 'ACA', 'PL', 'ILA', 'ORB', 'AI', 'RSP', 'PTLp', 'TEa', 'PERI', 'ECT', 'VISC', 'PIR', 'HPF', 'CTXsp', 'STR', 'PAL', 'THAL/HY']
    plt.plot([separation+shift, separation+shift], [0*np.ones(separation.shape), 2*np.ones(len(separation))], '-k', linewidth=.5) #vertical line
    plt.plot([separation+shift+nS/2, separation+shift+nS/2], [0*np.ones(separation.shape), 2*np.ones(len(separation))], '-k', linewidth=.5) #vertical line
    delta_sep = np.diff(separation, prepend=-1) / 2
    plt.xticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + nS/2]), xtick_label + xtick_label, fontsize=6, rotation='vertical')
    plt.xlim([-.5,nS-.5])
    plt.ylim([0,1.1])
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    
    plt.savefig('%sincominVSoutgoing.svg' % (pwd_tosave), format='svg')
    plt.close()
    
    if permute_corr:
        #if permutCorr: consider is sign VS. is not significant
        filename = 'in_out_binary_coupling_perNetwork_bar'
    else:
        filename = 'in_out_coupling_perNetwork_bar'
        
    X = [s_in.T, s_out.T]
    X_network = mean_perNetwork_acrossSbj(X, pwd_tosave=pwd_tosave, filename=filename)
    
    #select network
    from function_library import compute_stat_test
    selected_network = range(4) #range(9) #LCn, DMNmidline, DMNpostl, SAL
    selected_data = [2,3]#[0,1] #incoming and outgoing X_network
    label_data = ['incoming coupling', 'outgoing coupling']
    
    for i_data, data_i in enumerate(selected_data):
        print '----------------------'
        print label_data[i_data]
        group = []
        data_totest = []
        for ntw_i in selected_network:
            group.append(np.ones(n_sbj)*ntw_i)
            data_totest.append(X_network[data_i][ntw_i])
        compute_stat_test(data_totest, np.hstack(group))
    
    filename_to_save = '%soutput_%s.pickle' % (pwd_tosave, n_sim)
    print "Saving dataset in", filename_to_save
    pickle.dump({'s_in': s_in,
                 's_out': s_out,
                 'X_network': X_network,
                 'n_roi_min': n_roi_min
                 },
                open(filename_to_save, 'w'),
                protocol = pickle.HIGHEST_PROTOCOL)    
    
    ### only ratio within
    X = [onlySCrate_within_in.T, onlySCrate_within_out.T]
    filename = 'onlyRate_within'
    X_network = mean_perNetwork_acrossSbj(X, pwd_tosave=pwd_tosave, filename=filename)
    
    #select network
    from function_library import compute_stat_test
    selected_network = range(4) #range(9) #LCn, DMNmidline, DMNpostl, SAL
    selected_data = [2,3]#[0,1] #incoming and outgoing X_network
    label_data = ['incoming', 'outgoing']
    
    for i_data, data_i in enumerate(selected_data):
        print '----------------------'
        print label_data[i_data]
        group = []
        data_totest = []
        for ntw_i in selected_network:
            group.append(np.ones(n_sbj)*ntw_i)
            data_totest.append(X_network[data_i][ntw_i])
        compute_stat_test(data_totest, np.hstack(group))
    
    filename_to_save = '%soutput_onlyRate_within_%s.pickle' % (pwd_tosave, n_sim)
    print "Saving dataset in", filename_to_save
    pickle.dump({'s_in': s_in,
                 's_out': s_out,
                 'X_network': X_network,
                 'n_roi_min': n_roi_min
                 },
                open(filename_to_save, 'w'),
                protocol = pickle.HIGHEST_PROTOCOL) 
    
    ## only rate across
    X = [onlySCrate_across_in.T, onlySCrate_across_out.T]
    filename = 'onlyRate_across'
    X_network = mean_perNetwork_acrossSbj(X, pwd_tosave=pwd_tosave, filename=filename)
    
    #select network
    from function_library import compute_stat_test
    selected_network = range(4) #range(9) #LCn, DMNmidline, DMNpostl, SAL
    selected_data = [2,3]#[0,1] #incoming and outgoing X_network
    label_data = ['incoming', 'outgoing']
    
    for i_data, data_i in enumerate(selected_data):
        print '----------------------'
        print label_data[i_data]
        group = []
        data_totest = []
        for ntw_i in selected_network:
            group.append(np.ones(n_sbj)*ntw_i)
            data_totest.append(X_network[data_i][ntw_i])
        compute_stat_test(data_totest, np.hstack(group))
    
    filename_to_save = '%soutput_onlyRate_across_%s.pickle' % (pwd_tosave, n_sim)
    print "Saving dataset in", filename_to_save
    pickle.dump({'s_in': s_in,
                 's_out': s_out,
                 'X_network': X_network,
                 'n_roi_min': n_roi_min
                 },
                open(filename_to_save, 'w'),
                protocol = pickle.HIGHEST_PROTOCOL) 

if __name__=='__main__':
    
    #strucutre functional coupling
    permute_corr = False
    symmetric = False
    n_roi_min = 15
    n_sim = [62131, 92131, 102131, 82131]
    if permute_corr:
        n_perm = 500
    else:
        n_perm = 2

    pwd_current = '%s/' % os.getcwd()
    pwd_data = '%sdata/' % pwd_current    
    if not(os.path.isdir('%sfigure' % pwd_current)):
        os.makedirs('%sfigure' % pwd_current)
    
    for i_sim , n_sim_i in enumerate(n_sim):
        

        dataset = loadmat('%sdcmECfc_singleMouse_sim%s.mat' % (pwd_data, n_sim_i), squeeze_me=True)
        fc = dataset['fc']
        A = dataset['A'] 
        A_sparse = dataset['A_sparse']
        nS, _, n_sbj = A_sparse.shape
        p_val = .05
        
        mtx_str_log, mtx_str = load_structural_filter(pwd_data) #to be transposed to keep coherent notation with EC
        mtx_str_log = mtx_str_log.T
        mtx_str = mtx_str.T
        
        ##############
        if symmetric:
            mtx_str_log = mtx_str_log + mtx_str_log.T
            A = A + np.transpose(A, [1,0,2])
            A_sparse = A_sparse + np.transpose(A_sparse, [1,0,2])
        
        ####  only link with sc anatomical connection
        s_in = np.zeros([n_sbj, nS])
        onlySCrate_within_in = np.zeros([n_sbj, nS])
        onlySCrate_across_in = np.zeros([n_sbj, nS])
        s_out = np.zeros([n_sbj, nS])
        onlySCrate_within_out = np.zeros([n_sbj, nS])
        onlySCrate_across_out = np.zeros([n_sbj, nS])
        s_auto = np.zeros([n_sbj, nS]) #level of symmetry in the matrix
        n_nnzeros_roi = np.zeros([n_sbj, nS, 3])
        
        for i_sbj in range(n_sbj):
            print "------------------------------"
            print "sbj: %s" % i_sbj
            s_i, p_i, n_i, s_small_i, p_small_i, onlyRate_within_i, onlyRate_across_i = structure_functional_coupling(mtx_str_log, A[:,:,i_sbj], sparsified=True, n_perm=n_perm, n_roi_min=n_roi_min) #mtx_str_log
            print "incoming: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            if permute_corr: #if permutation test: binary matrix is sign vs. not sign.
                s_in[i_sbj] = p_small_i < p_val #s_in[i_sbj][p_small_i < p_val] = s_small_i[p_small_i < p_val]
            else:
                s_small_i[s_small_i == 1] = .99  #remove 1 corr to not have inf after Zfisher
                s_in[i_sbj] = np.arctanh(s_small_i) #s_small_i #
            onlySCrate_within_in[i_sbj] = onlyRate_within_i
            onlySCrate_across_in[i_sbj] = onlyRate_across_i
            n_nnzeros_roi[i_sbj, :, 0] = n_i
            
            s_i, p_i, n_i, s_small_i, p_small_i, onlyRate_within_i, onlyRate_across_i = structure_functional_coupling(mtx_str_log, A[:,:,i_sbj], incoming=False, sparsified=True, n_perm=n_perm, n_roi_min=n_roi_min) #mtx_str_log
            print "outgoing: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            if permute_corr:
                s_out[i_sbj] = p_small_i < p_val #s_out[i_sbj][p_small_i < p_val] = s_small_i[p_small_i < p_val]
            else:
                s_small_i[s_small_i == 1] = .99
                s_out[i_sbj] = np.arctanh(s_small_i) #s_small_i #
            n_nnzeros_roi[i_sbj, :, 1] = n_i
            onlySCrate_within_out[i_sbj] = onlyRate_within_i
            onlySCrate_across_out[i_sbj] = onlyRate_across_i
            
            s_i, p_i, n_i = symm_level(A_sparse[:,:,i_sbj], sparsified=True)
            print "asym level: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            s_auto[i_sbj][p_i < p_val] = s_i[p_i < p_val]
            #s_out[i_sbj] = s_i
            n_nnzeros_roi[i_sbj, :, 2] = n_i
        
        #pdb.set_trace()
        folder_name = 'sc_ec'
        if not(os.path.isdir('%sfigure/%s' % (pwd_current, n_sim_i))):
            os.makedirs('%sfigure/%s' % (pwd_current, n_sim_i))
        if not(os.path.isdir('%sfigure/%s/%s' % (pwd_current, n_sim_i, folder_name))):
            os.makedirs('%sfigure/%s/%s' % (pwd_current, n_sim_i, folder_name))
        pwd_tosave = '%sfigure/%s/%s/' % (pwd_current, n_sim_i, folder_name)
        
        plot_save_results(s_in, s_out, onlySCrate_within_in, onlySCrate_within_out, onlySCrate_across_in, onlySCrate_across_out, pwd_tosave, permute_corr, n_sim_i)
    
        ############ only link with ec connection
        s_in = np.zeros([n_sbj, nS])
        onlySCrate_within_in = np.zeros([n_sbj, nS])
        onlySCrate_across_in = np.zeros([n_sbj, nS])
        s_out = np.zeros([n_sbj, nS])
        onlySCrate_within_out = np.zeros([n_sbj, nS])
        onlySCrate_across_out = np.zeros([n_sbj, nS])
        s_auto = np.zeros([n_sbj, nS]) #level of symmetry in the matrix
        n_nnzeros_roi = np.zeros([n_sbj, nS, 3])
        
        for i_sbj in range(n_sbj):
            print "------------------------------"
            print "sbj: %s" % i_sbj
            s_i, p_i, n_i, s_small_i, p_small_i, onlyRate_within_i, onlyRate_across_i = functional_structure_coupling(mtx_str_log, A[:,:,i_sbj], sparsified=True, n_perm=n_perm, n_roi_min=n_roi_min) #mtx_str_log
            print "incoming: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            if permute_corr: #if permutation test: binary matrix is sign vs. not sign.
                s_in[i_sbj] = p_small_i < p_val #s_in[i_sbj][p_small_i < p_val] = s_small_i[p_small_i < p_val]
            else:
                s_small_i[s_small_i == 1] = .99  #remove 1 corr to not have inf after Zfisher
                s_in[i_sbj] = np.arctanh(s_small_i) #s_small_i #
            onlySCrate_within_in[i_sbj] = onlyRate_within_i
            onlySCrate_across_in[i_sbj] = onlyRate_across_i
            n_nnzeros_roi[i_sbj, :, 0] = n_i
            
            s_i, p_i, n_i, s_small_i, p_small_i, onlyRate_within_i, onlyRate_across_i = functional_structure_coupling(mtx_str_log, A[:,:,i_sbj], incoming=False, sparsified=True, n_perm=n_perm, n_roi_min=n_roi_min) #mtx_str_log
            print "outgoing: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            if permute_corr:
                s_out[i_sbj] = p_small_i < p_val #s_out[i_sbj][p_small_i < p_val] = s_small_i[p_small_i < p_val]
            else:
                s_small_i[s_small_i == 1] = .99
                s_out[i_sbj] = np.arctanh(s_small_i) #s_small_i #
            n_nnzeros_roi[i_sbj, :, 1] = n_i
            onlySCrate_within_out[i_sbj] = onlyRate_within_i
            onlySCrate_across_out[i_sbj] = onlyRate_across_i
            
            s_i, p_i, n_i = symm_level(A_sparse[:,:,i_sbj], sparsified=True)
            print "asym level: %f - min %f - max %f" % (np.mean(n_i), np.min(n_i), np.max(n_i))
            s_auto[i_sbj][p_i < p_val] = s_i[p_i < p_val]
            #s_out[i_sbj] = s_i
            n_nnzeros_roi[i_sbj, :, 2] = n_i
        
        #pdb.set_trace()
        folder_name = 'ec_sc'
        if not(os.path.isdir('%sfigure/%s' % (pwd_current, n_sim_i))):
            os.makedirs('%sfigure/%s' % (pwd_current, n_sim_i))
        if not(os.path.isdir('%sfigure/%s/%s' % (pwd_current, n_sim_i, folder_name))):
            os.makedirs('%sfigure/%s/%s' % (pwd_current, n_sim_i, folder_name))
        pwd_tosave = '%sfigure/%s/%s/' % (pwd_current, n_sim_i, folder_name)
        
        plot_save_results(s_in, s_out, onlySCrate_within_in, onlySCrate_within_out, onlySCrate_across_in, onlySCrate_across_out, pwd_tosave, permute_corr, n_sim_i)
    #end sim
    
    # plot IOI ration within
    from fig4_IOIratio_within import plot_ratio_within 
    plot_ratio_within('%sfigure/' % pwd_current)
    
    #plot IOI ration across
    from fig5_IOIratio_between import plot_ratio_between 
    plot_ratio_between('%sfigure/' % pwd_current)    