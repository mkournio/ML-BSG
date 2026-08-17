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
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    classification_report,
    r2_score, 
    mean_absolute_error
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

def transform_variables(input_variables, **kwargs):
    
    variables= input_variables.rename(str,axis="columns")
    scaled_variables = apply_scaler(variables,**kwargs)
    pca_variables = apply_pca(scaled_variables,**kwargs)
    
    return pca_variables       

class Features(DataFrame):
    
    def __init__(self,
                 *args,
                 **kwargs):   
        
        super().__init__(*args,**kwargs)
        
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
            self[n] = c
            
        for c in calc_keys:
            if c == 'MAD_RATIO':
                self[c] = self['MAD_RAW'] / self['MAD']
            if c == 'JH':
                self[c] = self['Jmag'] - self['Hmag']
            if c == 'HK':
                self[c] = self['Hmag'] - self['Kmag']
            if c == 'KW4':
                self[c] = self['Kmag'] - self['W4mag']
            if c == 'W14':
                self[c] = self['W1mag'] - self['W4mag']              
            if c == 'W24':
                self[c] = self['W2mag'] - self['W4mag']
            if c == 'W34':
                self[c] = self['W3mag'] - self['W4mag']  
            if c == 'Q_JHK':
                self[c] = self['JH'] - 1.7 * self['HK']
                    
        for c in log_convert :
            if c in self.columns:
                if c in ['W0','R0']:
                    self[c] = np.log10(1e+7 * self[c] + 1.)
                else:
                    self[c] = np.log10(self[c])
                    
        if 'R0' in self.columns:
            nan_mask = self['R0'] < 0.8
            try:
                self['TAU'][nan_mask] = np.nan
                self['GAMMA'][nan_mask] = np.nan
            except:
                pass
            
        #print(self[['STAR','SpC','Tmag']])                
                
        if 'save_output' in kwargs:
            self.to_csv(kwargs['save_output'],index=False)    
      
        return self    
   
    def _aggregate(self, 
                  cols, 
                  group_by,
                  mode='median',
                  **kwargs):
        
        if mode == 'none':
            
            return self.copy()
        
        df_copy = self.copy()
        
        agg_cols ={}
        for c in cols:#:
            if mode == 'median':
                agg_cols[c] = np.nanmedian
            elif mode == 'mean':
                agg_cols[c] = np.nanmean
                
        agg_self = df_copy.groupby(group_by,as_index=False).agg(agg_cols)
        
        if 'save_output' in kwargs:
            agg_self.to_csv(kwargs['save_output'],index=False)

        return agg_self 
    
    def pair_plot(self,
                  plot_cols,
                  hue,
                  corner = True,
                  split_cand = True,
                  aggregate_type = 'none',
                  outlier_sigma = 5.,
                  **kwargs):
        
        self_l = self._aggregate(plot_cols,
                                 group_by = ['STAR','SpC'],
                                 mode=aggregate_type)            

        
        if isinstance(outlier_sigma, (int,float)):            
            for c in plot_cols :
                if self_l[c].dtype == np.float64 :
                    Q1 = np.nanquantile(self_l[c], 0.25)
                    Q3 = np.nanquantile(self_l[c], 0.75)
                    IQR = Q3 - Q1
                    nan_mask =  (self_l[c] < (Q1 - (outlier_sigma * IQR))) | (self_l[c] > (Q3 + (outlier_sigma * IQR)))
                    self_l[c][nan_mask] = np.nan  
        
        import seaborn as sns
        
        ax = sns.pairplot(self_l,
                          hue = hue,
                          vars=plot_cols,
                          corner = corner,
                          **kwargs)     
        
        return ax  
    
    def pair_plot_single(self,
                         pair_cols,
                         aggregate_type = 'median',
                         **kwargs):
        
        self_l = self._aggregate(pair_cols,
                                 group_by = ['STAR','SpC'],
                                 mode=aggregate_type)
        
        variables = self_l[pair_cols]
       
        fig, ax = plt.subplots(figsize=(8, 6))
        for s in set(self_l['STAR']):
            
            mask = self_l['STAR'] == s
            star = variables[mask]
            spt = self_l['SpC'][mask].iloc[0]
            
            if 'YHG' in spt:
                c = 'g'
                m = 'o'
            elif 'B[e]SG' in spt :
                c = 'b' 
                m = 's'
            elif 'LBV' in spt:
                c = 'orange'
                m = '^'

                
            ax.plot(star.iloc[:,0],star.iloc[:,1], m, c = c, ms = 9)
            if not aggregate_type == 'none':
                ax.text(star.iloc[:,0]+0.03,star.iloc[:,1]+0.03, s, size=8)
            if '?' in spt:
                ax.plot(star.iloc[:,0],star.iloc[:,1], m, c = 'w', ms = 4)
                
        ax.set_xlabel(pair_cols[0], fontsize=10)
        ax.set_ylabel(pair_cols[1], fontsize=10)            
        
        return ax  
    
    def knn_classify(self,
                   var_cols,
                   class_col = 'SpC',
                   split_cand = True,
                   aggregate_type = 'none',
                   n_perm = 0 ,                   
                   **kwargs):
        
        self_l = self._aggregate(var_cols,
                                 group_by = ['STAR','SpC'],
                                 mode=aggregate_type)         
        if split_cand:
            self_l['SpC'] = [x.replace('?','') for x in self_l['SpC'] ]           
                     
        # for c in var_cols + cbar_col:    
            #     mask = np.isnan(self_l[c])
            #     self_l = self_l[~mask]       
        x = transform_variables(self_l[var_cols],**kwargs)      
        y = self_l[class_col].to_numpy()
        
        knn = KNeighborsClassifier(n_neighbors=3,weights="distance",metric="euclidean")
        loo = LeaveOneOut()
        y_pred = cross_val_predict(knn, x, y, cv=loo)
        #for i,j in zip(y,y_pred): print(i,j)        
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
                    var_cols,
                    regress_col = 'Tmag',
                    split_cand = True,
                    aggregate_type = 'none',
                    n_perm = 0 ,                   
                    **kwargs):
        
        self_l = self._aggregate(var_cols + [regress_col],
                                 group_by = ['STAR','SpC'],
                                 mode=aggregate_type)         
        if split_cand:
            self_l['SpC'] = [x.replace('?','') for x in self_l['SpC'] ] 
            
        x = transform_variables(self_l[var_cols],**kwargs)
        y = self_l[regress_col].to_numpy()
        
        knn = KNeighborsRegressor(n_neighbors=3,weights="distance",metric="euclidean")
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

        
        
    
    
        
        

    
    def umap_plot(self,
                  var_cols,
                  ax = None,
                  cbar_col = [],    
                  aggregate_type = 'none',
                  **kwargs):
        
        import umap
        
        self_l = self._aggregate(var_cols + cbar_col,
                                 group_by = ['STAR','SpC'],
                                 mode=aggregate_type)
        
       # for c in var_cols + cbar_col:    
       #     mask = np.isnan(self_l[c])
       #     self_l = self_l[~mask]
            
        variables = self_l[var_cols]
        
        variables= variables.rename(str,axis="columns") 
        scaled_variables = apply_scaler(variables,**kwargs)
        pca_variables = apply_pca(scaled_variables,**kwargs)
        
        umap_reducer = umap.UMAP(n_components = 2, random_state=0, 
                                 n_neighbors = kwargs.get('n_neighbors',12),
                                 min_dist = kwargs.get('min_dist',0.1))
        umap_variables = umap_reducer.fit_transform(pca_variables)
        
        if ax is None:
         _, ax = plt.subplots(figsize=(8, 6))
         
        if len(cbar_col) == 1:            
            cbc = self_l[cbar_col]
            mn = min(cbc.values); mx = max(cbc.values)
            if cbar_col[0] == 'Q_JHK': mn = -1.3
            if cbar_col[0] == 'CROWDSAP': mn = 0.85
            cmap, cnorm = colorbar(mn,mx,cbar='rainbow')
        
        for s in set(self_l['STAR']):
            
            umask = self_l['STAR'] == s
            ustar = umap_variables[umask]
            spt = self_l['SpC'][umask].iloc[0]
            
            if 'YHG' in spt:
                c = 'r'
                m = 'o'
            elif 'B[e]SG' in spt :
                c = 'b' 
                m = 's'
            elif 'LBV' in spt:
                c = 'green'
                m = '^'
                
            ax.plot(ustar[:,0],ustar[:,1],'k',lw=0.1) 
            if len(cbar_col) == 1:
                c = cmap(cnorm(cbc[umask].iloc[0]))
                
            ax.plot(ustar[:,0],ustar[:,1], m, c = c, ms = 10)
            ax.text(ustar[0,0]-0.05,ustar[0,1]-0.18, s, size=6)
            if '?' in spt:
                ax.plot(ustar[:,0],ustar[:,1], m, c = 'w', ms = 4)
                
        ax.set_xlabel('UMAP 1', fontsize=10)
        ax.set_ylabel('UMAP 2', fontsize=10)
        if len(cbar_col) == 1:
            extend = 'neither'
            if cbar_col[0] == 'Q_JHK' or cbar_col[0] == 'CROWDSAP': extend = 'min'
            plt.colorbar(cm.ScalarMappable(norm=cnorm, cmap=cmap), ax=ax,
                         orientation = 'horizontal',
                         extend = extend)          
        
        return ax
    
        
        
    
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