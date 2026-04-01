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
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA


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
                         log_convert = [],
                         meta_keys = [],
                         **kwargs):
       
        if len(time_keys) == 0 and len(freq_keys) == 0 and len(rn_keys) == 0:
            return        

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
                    
                    if len(meta_keys) != 0:
                        for k_index, k in enumerate(meta_keys):
                            if k in input_cat.columns:
                                meta_values[k_index] = input_cat[k][star_index]    
                        
                    hdu = get_hdu_from_keys(ff, SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINSIZE = str(bin_size))[0] 
                    hdr = hdu.header
                    if len(time_keys) != 0:               
                        for k_index, k in enumerate(time_keys):                            
                            if k in hdr:
                                t_values[k_index] = hdr[k]
                                
                    hdu_f = get_hdu_from_keys(ff, SECTOR = s, HDUTYPE = 'FREQUENCIES', BINSIZE = str(bin_size))[0]
                    hdr = hdu_f.header                                 
                    if len(freq_keys) != 0:
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
            
        for c in log_convert :
            if c in self.columns:
                self[c] = np.log10(self[c])
        
        if 'save_output' in kwargs:
            self.to_csv(kwargs['save_output'],index=False)    
      
        return self    
   
    def _aggregate(self, 
                  cols, 
                  group_by,
                  mode='median',
                  **kwargs):
        
        agg_cols ={}
        for c in cols:#:
            if mode == 'median':
                agg_cols[c] = np.nanmedian
            elif mode == 'mean':
                agg_cols[c] = np.nanmean
                
        agg_self = self.groupby(group_by,as_index=False).agg(agg_cols)
        
        if 'save_output' in kwargs:
            agg_self.to_csv(kwargs['save_output'],index=False)

        return agg_self 
    
    def pair_plot(self,
                  plot_cols,
                  hue,
                  corner = True,
                  split_cand = True,
                  outlier_sigma = 5,
                  **kwargs):
        
        self_l = self.copy()
        
        if split_cand:
            self_l['SpC'] = [x.replace('?','') for x in self_l['SpC'] ]
            mask = ['?' not in x for x in self_l['SpC']]
            self_l = self_l[mask]
        
        if isinstance(outlier_sigma, (int,float)):            
            for c in plot_cols :
                if self_l[c].dtype == np.float64 :
                    Q1 = np.nanquantile(self_l[c], 0.25)
                    Q3 = np.nanquantile(self_l[c], 0.75)
                    IQR = Q3 - Q1
                    nan_mask =  (self_l[c] < (Q1 - (outlier_sigma * IQR))) | (self_l[c] > (Q3 + (outlier_sigma * IQR)))
                    self_l[c][nan_mask] = None  
        
        import seaborn as sns
        
        ax = sns.pairplot(self_l,
                          hue = hue,
                          vars=plot_cols,
                          corner = corner,
                          **kwargs)     
        
        return ax  
    
    def umap_plot(self,
                  var_cols,
                  meta_cols = [],    
                  aggregate = False,
                  **kwargs):
        
        import umap
        
        if aggregate:
            self_l = self._aggregate(cols = var_cols,
                             group_by = ['STAR','SpC']) 
        else:
            self_l = self.copy()
        
        for c in var_cols:    
            mask = np.isinf(self_l[c]) | np.isnan(self_l[c])
            self_l = self_l[~mask]
            
        variables = self_l[var_cols]
        scaled_variables = apply_scaler(variables,**kwargs)
        pca_variables = apply_pca(scaled_variables,**kwargs)
        
        umap_reducer = umap.UMAP(n_components = 2, random_state=0, 
                                 n_neighbors = kwargs.get('n_neighbors',12),
                                 min_dist = kwargs.get('min_dist',0.1))
        umap_variables = umap_reducer.fit_transform(pca_variables)
        
        fig, ax = plt.subplots(figsize=(8, 6))
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
                c = 'lime'
                m = '^'
                
            ax.plot(ustar[:,0],ustar[:,1],c='0.8')    
            ax.plot(ustar[:,0],ustar[:,1], m, c = c, ms = 9)
            ax.text(ustar[0,0]+0.03,ustar[0,1]+0.03, s, size=8)
            if '?' in spt:
                ax.plot(ustar[:,0],ustar[:,1], m, c = 'w', ms = 4)
                
        ax.set_xlabel('UMAP 1', fontsize=10)
        ax.set_ylabel('UMAP 2', fontsize=10)            
        
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