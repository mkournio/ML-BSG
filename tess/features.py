#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 21:32:52 2025

@author: michalis
"""

import lightkurve as lk
import numpy as np
from methods.tools import *
from astropy.io import fits
from scipy import stats
from scipy.stats.mstats import pearsonr, spearmanr, describe

from methods.functions import *
from methods.plot import GridTemplate, colorbar
from astropy.table import Table
from tables.io import tab_to_csv
from pandas import DataFrame
import pandas as pd
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
import matplotlib.cm as cm
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
import seaborn as sns
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.neighbors import NearestNeighbors, LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    classification_report,
    r2_score, 
    mean_absolute_error
)

import umap
from sklearn.cluster import AgglomerativeClustering, HDBSCAN
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score
)

def apply_scaler(df,
           scaler_type = 'standard',
           **kwargs):
    
    if scaler_type == 'standard':
        scaler = StandardScaler()
    elif scaler_type == 'robust':
        scaler = RobustScaler()
        
        
    feats_scaled = scaler.fit_transform(df)
    
    return feats_scaled

def apply_pca(df,
        pca_components = 5,
        **kwargs):
    
    pca = PCA(n_components=pca_components)
    pca.fit(df)
    feats_pca = pca.transform(df)
    
    cumsum_variance = np.cumsum(pca.explained_variance_ratio_)
    print(f"Total explained variance: {pca.explained_variance_ratio_.sum():.3f}")
    print(f"Cumsum variance: {cumsum_variance}")
    
    return feats_pca    

class Features:
    
    def __init__(self,
                 *args,
                 **kwargs):   
        
        self.df = pd.DataFrame(*args, **kwargs)        
        
        return
        
    def _validate(self):
        pass
        
    def get_from_sectors(self,
                         input_cat,
                         time_keys = [], 
                         freq_keys = [], 
                         rn_keys = [], 
                         bin_size = '10m',
                         calc_keys = [],
                         log_convert = [],
                         **kwargs):
       
        if len(time_keys) == 0 and len(freq_keys) == 0 and len(rn_keys) == 0:
            return  
        
        meta_keys = list(input_cat.columns)
        stars = input_cat['STAR']  
        
        column_names = np.concatenate((meta_keys,time_keys,freq_keys,rn_keys), axis=0)
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083            
        
        array = []        
        for star_index, star in enumerate(stars):

            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))                
                sectors = get_sectors_from_hdulist(ff)
                
                for s in sectors:
                    
                    t_values = np.full(len(time_keys), np.nan)
                    f_values = np.full(len(freq_keys), np.nan)
                    r_values = np.full(len(rn_keys), np.nan)
                    meta_values = np.full(len(meta_keys), np.nan, dtype='object') 
                    
                    for k_index, k in enumerate(meta_keys):
                        meta_values[k_index] = input_cat[k][star_index]    
                        
                    hdu = get_hdu_from_keys(ff, SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINSIZE = str(bin_size))[0] 
                    hdr = hdu.header
                    if len(time_keys) != 0:               
                        for k_index, k in enumerate(time_keys):                            
                            if k in hdr:
                                t_values[k_index] = hdr[k]
                                
                    if len(freq_keys) != 0:
                        hdu_f = get_hdu_from_keys(ff, SECTOR = s, HDUTYPE = 'FREQUENCIES', BINSIZE = str(bin_size))[0]
                        hdr = hdu_f.header
                        for k_index, k in enumerate(freq_keys):
                            if k in hdr:
                                f_values[k_index] = hdr[k]  
                        if len(rn_keys) != 0:
                            rn_model = hdu_f.data[-1]
                            for k_index, k in enumerate(rn_keys):                            
                                r_values[k_index] = rn_model[k]          
                                
                    sect_values = np.concatenate((meta_values,t_values,f_values,r_values), axis=0)
                    array.append(sect_values)
                    
                ff.close()
                
        t_array = list(map(list, zip(*array)))
        for n, c in zip(column_names,t_array):
            self.df[n] = c
            
        for c in calc_keys:
            if c == 'MAD_RATIO':
                self.df[c] = self.df['MAD_RAW'] / self.df['MAD']
            if c == 'JH':
                self.df[c] = self.df['Jmag'] - self.df['Hmag']
            if c == 'HK':
                self.df[c] = self.df['Hmag'] - self.df['Kmag']
            if c == 'KW4':
                self.df[c] = self.df['Kmag'] - self.df['W4mag']
            if c == 'W14':
                self.df[c] = self.df['W1mag'] - self.df['W4mag']              
            if c == 'W24':
                self.df[c] = self.df['W2mag'] - self.df['W4mag']
            if c == 'W34':
                self.df[c] = self.df['W3mag'] - self.df['W4mag']  
            if c == 'Q_JHK':
                self.df[c] = self.df['JH'] - 1.7 * self.df['HK']
                    
        for c in log_convert :
            if c in self.df.columns:
                if c in ['W0','R0']:
                    self.df[c] = np.log10(1e+7 * self.df[c] + 1)
                else:
                    self.df[c] = np.log10(self.df[c])
                    
        if 'R0' in self.df.columns:
            nan_mask = self.df['R0'] < 0.8
            try:
                self.df['TAU'][nan_mask] = np.nan
                self.df['GAMMA'][nan_mask] = np.nan
            except:
                pass
            
       # print(self[['STAR','SpC','Tmag']])                
                
        if 'save_output' in kwargs:
            self.df.to_csv(kwargs['save_output'],index=False)            
      
        return self.df    
    
    def merge_cand(self):
        
        self.df['SpC'] = [x.replace('?','') for x in self.df['SpC'] ]
        
        return        
   
    def aggregate(self, 
                  agg_cols, 
                  group_by,                  
                  agg_type='median',
                  split_cand = False,
                  **kwargs):
        
        if agg_type == 'none':
            
            return self.df       
       
        agg_cols_dict = {}
        for c in agg_cols:#:
            if agg_type == 'median':
                agg_cols_dict[c] = np.nanmedian
            elif agg_type == 'mean':
                agg_cols_dict[c] = np.nanmean
                
        self.df = self.df.groupby(group_by,as_index=False).agg(agg_cols_dict)

        return self
    
    def corner_plot(self,
                  plot_cols,
                  hue,
                  corner = True,
                  outlier_sigma = 3.,
                  **kwargs):        
        
        if isinstance(outlier_sigma, (int,float)):            
            for c in plot_cols :
                if self.df[c].dtype == np.float64 :
                    Q1 = np.nanquantile(self.df[c], 0.25)
                    Q3 = np.nanquantile(self.df[c], 0.75)
                    IQR = Q3 - Q1
                    nan_mask =  (self.df[c] < (Q1 - (outlier_sigma * IQR))) | (self.df[c] > (Q3 + (outlier_sigma * IQR)))
                    self.df[c][nan_mask] = np.nan   
        
        plt.rcParams.update(**cornerplot_kwargs)
        
        rcols = {}
        for c in plot_cols:
            rcols[c] = st(c)
        loc_df = self.df.rename(columns=rcols)        
        ax = sns.pairplot(loc_df,
                          hue = hue,
                          vars=[st(s) for s in plot_cols],
                          corner = corner,
                          palette= ['b','g','r'], 
                          markers=['s', '^', 'o'], 
                          **kwargs) 
        ax._legend.set_title("Class")
        sns.move_legend(ax, loc='center', bbox_to_anchor=(.70, .70), frameon=True)
        
        return ax  
    
    def pair_plot(self,                  
                  pair_c,
                  ax = None,
                  star_labels = False,
                  **kwargs):     
        
        variables = self.df[pair_c]
        spear_val = spearmanr(x=variables.iloc[:,0].values,y=variables.iloc[:,1].values)
        print(spear_val.statistic, spear_val.pvalue)
        
        if ax is None:
         _, ax = plt.subplots(figsize=(8, 6))
         
        plt.rcParams.update(**pairplot_kwargs)

       
        for s in set(self.df['STAR']):
            
            mask = self.df['STAR'] == s
            star = variables[mask]
            spt = self.df['SpC'][mask].iloc[0]   
            x = star.iloc[:,0].values
            y = star.iloc[:,1].values
            
            ax.plot(x, y, CLASS_M[spt], c = CLASS_C[spt])
            if star_labels:
                ax.text(x+0.015,star.iloc[:,1]+0.01, st(s))
            if '?' in spt:
                ax.plot(x, y, CLASS_M[spt], c = 'w', ms = 4)         
        
                
        ax.set_xlabel(st(pair_c[0]))
        ax.set_ylabel(st(pair_c[1]))            
        
        return ax  
    
    def _pca_transform(self, var_cols = [], **kwargs):
        
        variables = self.df[var_cols]   
        variables = variables.rename(str,axis="columns")
        scaled_variables = apply_scaler(variables,**kwargs)
        pca_variables = apply_pca(scaled_variables,**kwargs)        
               
        return pca_variables
    
    def umap_plot(self,
                  ax = None,
                  cbar_col = [],
                  umap_n = 5,
                  umap_min_dist = 0.1,
                  show_mutual = False,
                  random_state = 0,
                  **kwargs):
        
        pca_variables = self._pca_transform(**kwargs)
        umap_reducer = umap.UMAP(n_components = 2, random_state = random_state, 
                                 n_neighbors = umap_n, min_dist = umap_min_dist)
        umap_variables = umap_reducer.fit_transform(pca_variables) 
        
        if show_mutual:
            knn_tab = self.nearest_neighbors(**kwargs) 
            #print(knn_tab)

        if ax is None:
         _, ax = plt.subplots(figsize=(8, 6))
         
        plt.rcParams.update(**umap_kwargs)
         
        if len(cbar_col) == 1:            
            cbc = self.df[cbar_col]
            mn = min(cbc.values); mx = max(cbc.values)
            if cbar_col[0] == 'Q_JHK': mn = -1.3
            if cbar_col[0] == 'CROWDSAP': mn = 0.85
            cmap, cnorm = colorbar(mn,mx,cbar='rainbow')
        
        for s in set(self.df['STAR']):
            
            umask = self.df['STAR'] == s
            spt = self.df['SpC'][umask].iloc[0]            
            ustar = umap_variables[umask]          
                    
            if len(cbar_col) == 1:
                c = cmap(cnorm(cbc[umask].iloc[0]))
            else:
                c = CLASS_C[spt] 
                
            if show_mutual:
              neighbors = knn_tab['NEIGHBOR'][knn_tab['STAR'] == s]
              for n in neighbors:
                  n_neighbors = knn_tab['NEIGHBOR'][knn_tab['STAR'] == n]
                  if s in n_neighbors.to_numpy():
                      n_ustar = umap_variables[self.df['STAR'] == n]
                      ax.plot([ustar[:,0],n_ustar[:,0]],
                              [ustar[:,1],n_ustar[:,1]],'k', ls = (0, (5, 10)), lw = 0.4, zorder=1)
                      
            ax.plot(ustar[:,0],ustar[:,1], CLASS_M[spt], c = c, zorder=2)
            if '?' in spt:
                ax.plot(ustar[:,0],ustar[:,1], CLASS_M[spt], c = 'w', ms = 4, zorder=2)
                
            yl0, yl1 = ax.get_ylim()
            xl0, xl1 = ax.get_xlim()            
            ax.text(ustar[:,0] + 0.01 * (xl1-xl0),
                    ustar[:,1] - 0.03 * (yl1-yl0), st(s)[0])
            #ax.plot(ustar[:,0],ustar[:,1],'k',lw=0.1) 
            
        ax.text(0.03,0.95, 'n = %s' % umap_n, fontsize = 17, transform = ax.transAxes)

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
      #  if len(cbar_col) == 1:
      #      extend = 'neither'
      #      if cbar_col[0] == 'Q_JHK' or cbar_col[0] == 'CROWDSAP': extend = 'min'
      #      plt.colorbar(cm.ScalarMappable(norm=cnorm, cmap=cmap), ax=ax,
      #                   orientation = 'horizontal',
      #                   extend = extend)          
        
        return ax
    
    
    def knn_classify(self,
                   kn = 3,
                   class_col = 'SpC',
                   exclude_labels = [],
                   split_cand = True,
                   n_perm = 0 ,                   
                   **kwargs):
        
        x = self._pca_transform(**kwargs)     
        y = self.df[class_col].to_numpy()
        
        knn = KNeighborsClassifier(n_neighbors = kn + 1,weights="distance",metric="euclidean")
        loo = LeaveOneOut()
        y_pred = cross_val_predict(knn, x, y, cv=loo)

        acc = accuracy_score(y, y_pred)
        bal_acc = balanced_accuracy_score(y, y_pred)
        
        labels = np.unique(y)
        cmt = confusion_matrix(y, y_pred, labels=labels)
        cm_df = pd.DataFrame(cmt, 
                             index=[f"true_{label}" for label in labels],
                             columns=[f"pred_{label}" for label in labels])
 
        print("Confusion matrix:",cm_df)
        print("Report:",classification_report(y, y_pred, labels=labels))
        print("Accuracy:", acc)
        print("Balanced accuracy:", bal_acc)      
    
        # Permutation test
        if n_perm > 0:
            rng = np.random.default_rng(42)
            bal_acc_perm = np.zeros(n_perm)
            for i in range(n_perm):
                y_perm = rng.permutation(y)
                y_perm_pred = cross_val_predict(knn, x, y_perm, cv=loo)
                bal_acc_perm[i] = balanced_accuracy_score(y_perm, y_perm_pred)
            p_value = (np.sum(bal_acc_perm >= bal_acc) + 1) / (n_perm + 1)
            
            print("Permutation mean:", np.mean(bal_acc_perm))
            print("Permutation std:", np.std(bal_acc_perm))
            print("Permutation p-value:", p_value)

        return
    
    def knn_regress(self,
                    kn = 3,
                    regress_col = 'Tmag',
                    exclude_labels = [],
                    n_perm = 0 ,                   
                    **kwargs):
        
        x = self._pca_transform(**kwargs)     
        y = self.df[regress_col].to_numpy()
        
        knn = KNeighborsRegressor(n_neighbors = kn + 1, weights="distance",metric="euclidean")
        loo = LeaveOneOut()
        y_pred = cross_val_predict(knn, x, y, cv=loo)
        #for i,j in zip(y,y_pred): print(i,j)   
        
        r2 = r2_score(y, y_pred)
        mae = mean_absolute_error(y, y_pred)
        print("R2:", r2)
        print("MAE:", mae)
        
        #Permutation test
        if n_perm > 0:
            rng = np.random.default_rng(42)
            r2_perm = np.zeros(n_perm)
            mae_perm = np.zeros(n_perm)
            
            for i in range(n_perm):
                y_perm = rng.permutation(y)
                y_perm_pred = cross_val_predict(knn, x, y_perm, cv=loo)
                r2_perm[i] = r2_score(y_perm, y_perm_pred)
                mae_perm[i] = mean_absolute_error(y_perm, y_perm_pred)
                
            p_r2 = (np.sum(r2_perm >= r2) + 1) / (n_perm + 1)
            p_mae = (np.sum(mae_perm <= mae) + 1) / (n_perm + 1)
            
            print("Permutation R2 mean:", np.mean(r2_perm))
            print("Permutation R2 std:", np.std(r2_perm))
            print("Permutation p-value R2:", p_r2)
            print("Permutation MAE mean:", np.mean(mae_perm))
            print("Permutation MAE std:", np.std(mae_perm))
            print("Permutation p-value MAE:", p_mae)
        
        return
    
    def hier_clustering(self,
                      k = 'scan',
                      d_thres = 3,
                      **kwargs):
        
        x = self._pca_transform(**kwargs)
        names = [st(star)[0] for star in self.df["STAR"]]
        label_colors = {}
        for n, c in zip(names,self.df["SpC"]):
            label_colors[n] = CLASS_C[c]
            
        threshold = 9

        if isinstance(k,int) and k > 1:
            clust = AgglomerativeClustering(n_clusters = k, 
                                            distance_threshold=None,
                                            linkage= "ward", 
                                            metric= "euclidean")
            labels = clust.fit_predict(x)
            self.df['HIER_K'] = labels
            #print(self.df[['STAR','SpC','HIER_K']])
            
            Z = linkage(x, method="ward", metric="euclidean")
            
            plt.rcParams.update(**dendro_kwargs)            
            plt.figure(figsize=(12, 5))
            dendrogram(Z, labels = names,
                       color_threshold = 0,
                       above_threshold_color = 'k',
                       leaf_rotation=90, 
                       leaf_font_size=12)
            plt.ylabel(r"Ward linkage distance")
            plt.title("Hierarchical clustering")
            #plt.axhline(y = threshold, c='k', ls='--')
            ax = plt.gca()
            ax_labels = ax.get_xmajorticklabels()
            for lbl in ax_labels:
                lbl.set_color(label_colors[lbl.get_text()])
            plt.tight_layout()
            plt.show()
            
        else:
            rows = []
            for k in range(2, 12):
                clust = AgglomerativeClustering(n_clusters = k, 
                                                linkage="ward", 
                                                metric="euclidean")
                labels = clust.fit_predict(x)
                rows.append({
                    "k": k,
                    "silhouette": silhouette_score(x, labels, metric="euclidean"),
                    "calinski_harabasz": calinski_harabasz_score(x, labels),
                    "davies_bouldin": davies_bouldin_score(x, labels)
                    })
            hier_scan = pd.DataFrame(rows)
            print(hier_scan)
        
        return   
    
    def hdbscan_clustering(self,
                           min_cluster_size = 2,
                           **kwargs):
        
        x = self._pca_transform(**kwargs)

        model = HDBSCAN(min_cluster_size = min_cluster_size)
        labels = model.fit_predict(x)
        
        self.df['HDBSCAN_K'] = labels
        print(self.df[['STAR','SpC','HDBSCAN_K']])        
       
        return
    
    def nearest_neighbors(self,
                      kn = 6,
                      **kwargs):
        
        x = self._pca_transform(**kwargs)
        names = self.df["STAR"].to_numpy()
        classes = self.df["SpC"].to_numpy()
        
        nn = NearestNeighbors(n_neighbors = kn + 1, metric="euclidean")
        nn.fit(x)
        
        distances, indices = nn.kneighbors(x)
        
        rows = []
        for i in range(len(names)):
            for rank in range(1, kn + 1):
                j = indices[i, rank]
                rows.append({
                "STAR": names[i],
                "KNN_RANK": rank,
                "NEIGHBOR": names[j],
                "PCA_DIST": distances[i, rank]
                })
                
        nn_table = pd.DataFrame(rows)
        
        return nn_table
    
    def local_outlier_factor(self,
                             kn = 6,
                             **kwargs):
        
        x = self._pca_transform(**kwargs)       
        
        lof = LocalOutlierFactor(n_neighbors = 5)
        labels = lof.fit_predict(x)
        self.df['LOF'] = labels
        print(self.df[['STAR','SpC','LOF']])
        
        return
    

    
        
        
    
    '''
    
    

            
    def scatter_plot(self, x, y, mode = 'matrix', invert = [], cbar = None, alpha = None, hold = False, **kwargs):
        
        if len(x) == 0 or len(x) == 0:
            
            return      
        
        ltab = self.data.copy()
        ltab['CMARK'] = 'k'
        ltab['AMARK'] = 1.

        
        for k in x + y :
            if 'log' in k:
                ltab[k] = np.log10(ltab[k.replace('log','')])      
        
        if mode == 'zip' and len(x) != len(y):
            raise IndexError('Zip mode is activated but key vectors have not same size.')

        if mode == 'matrix':
            size_grid = {'rows_page' : len(x), 'cols_page' : len(y)}
        elif mode == 'zip':
            size_grid = {'rows_page' : 1, 'cols_page' : len(y)}
            
        g = GridTemplate(fig_xlabel='', fig_ylabel='', params = PLOT_PARAMS['cr'], mode = mode,
                         row_labels= x, col_labels= y, **dict(kwargs,**size_grid))
        
        if cbar is not None:            
            try:
                cbar_key = cbar[0]
                cbar_l, cbar_u = cbar[1:]
            except:
                raise Exception('Define a colorbar as [tab_key,low_val, high_val].')
            
            cmap, cnorm = colorbar(vmin=cbar_l, vmax=cbar_u)
            ltab['CMARK'] = np.array([cmap(cnorm(k)) for k in ltab[cbar_key]])
            
        if alpha is not None:
            try:
                a_key = alpha[0]
                a_range = alpha[1:]
            except:
                raise Exception('Define alpha as [tab_key,alpha range].')
                
            i = 0
            a_range = [0.999*np.nanmin(ltab[a_key])] + sorted(a_range) + [np.nanmax(ltab[a_key])]     
            a_val = np.linspace(0.4, 1., len(a_range) - 1)
            
            while i < len(a_range)-1:
                mask = (a_range[i] < ltab[a_key]) & (ltab[a_key] <= a_range[i+1])
                ltab['AMARK'][mask] = a_val[i]
                #print(a_range[i],a_val[i],a_range[i+1])

                i += 1
        
                
        if mode == 'zip':
            
            for c1, c2 in zip(x,y) :
                
                ax = g.GridAx()
                self._panel(ax,ltab,c1,c2)
                
                if c1 in invert: ax.invert_yaxis()
                if c2 in invert: ax.invert_xaxis()
                
        elif mode == 'matrix':
            
            for c1 in x:

                for c2 in y:

                    ax = g.GridAx()
                    self._panel(ax,ltab,c1,c2)
                    
                    if c1 in invert: 
                        ax.invert_yaxis()
                        invert.remove(c1)
                    
                    if c2 in invert: 
                        ax.invert_xaxis()
                        invert.remove(c2)
                        
        if 'cmap' in locals(): 
            g.add_colorbar(cnorm,cmap,label=cbar_key,extend='min')    
        
        if not hold:
            g.close_plot()
            
        return               

    def _panel(self, ax, tab, k1, k2, **kwargs):
                
        if ax == None:
            
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots()
            
            
        LBVs = [('LBV' in x) & ('?' not in x) for x in tab['SpC']]
        BREs = [('B[e]SG' in x) & ('?' not in x) for x in tab['SpC']]
      #  SMB = [(x == 'SMB') and ('MW' in y) for x,y in zip(tab['REF'],tab['GAL'])]
        MW  = [('MW' in x) for x in tab['GAL']]
        MC = [('LMC' in x) or ('SMC' in x) for x in tab['GAL']]
        
        e_kwargs = {'elinewidth' : 0.5, 'capsize' : 0, 'ls' : 'none'}        
        
        ax.scatter(tab[k2][MC],tab[k1][MC],c=tab['CMARK'][MC],alpha=tab['AMARK'][MC],**plot_mcs)        
        ax.scatter(tab[k2][MW],tab[k1][MW],c=tab['CMARK'][MW],alpha=tab['AMARK'][MW],**plot_mw)  
        
        ax.plot(tab[k2][LBVs],tab[k1][LBVs],**plot_LBV)
        ax.plot(tab[k2][BREs],tab[k1][BREs],**plot_BREs)
        
        return
    
        
    def get_from_primary_headers(self, hdr_keys, update_table = True):         
      
        if len(hdr_keys) == 0:
            return

        elif any([x in self.data.columns for x in hdr_keys]):
            raise Exception('One of header keys already exists as input column. Aborting.')

        print('Creating columns from header keys: {}'.format(hdr_keys))

        stars = self.input_cat['STAR']
        array = np.full((len(stars),2*len(hdr_keys)), np.nan)
        
        for star_index, star in enumerate(stars):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0])) 
                hdr = ff[0].header                
                
                s_keys = np.full(len(hdr_keys), np.nan)
                s_keys_err = np.full(len(hdr_keys), np.nan) 
                
                for k_index, k in enumerate(hdr_keys):
                    if k in hdr:
                        s_keys[k_index] = hdr[k]                        
                        try:
                            s_keys_err[k_index] = hdr.comments[k]
                        except:
                            pass

                array[star_index] = np.concatenate((s_keys,s_keys_err), axis=0) 

                ff.close()
                
        t_array = list(map(list, zip(*array))) 
        
        hdr_keys_err = ['e_'+x for x in hdr_keys]
        columns_names = np.concatenate((hdr_keys,hdr_keys_err), axis=0)
        
        if  update_table :
            
            self.data.add_columns(t_array,names=columns_names)
            
            return 
        
        else:           
            
            return Table(t_array,names=columns_names)''' 