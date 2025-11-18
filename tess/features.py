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

class Features(object):
    
    def __init__(self,
                 data, 
                 **kwargs):
        
        self._validate()
        self.data = data
        
    def _validate(self):
        pass
    
    def cols_from_headers(self, hdr_keys, update_table = True):         
      
        if len(hdr_keys) == 0:
            return

        elif any([x in self.data.columns for x in hdr_keys]):
            raise Exception('One of header keys already exists as input column. Aborting.')

        print('Creating columns from header keys: {}'.format(hdr_keys))

        stars = self.data['STAR']
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
                
        k_array = list(map(list, zip(*array))) 
        
        hdr_keys_err = ['e_'+x for x in hdr_keys]
        columns_names = np.concatenate((hdr_keys,hdr_keys_err), axis=0)
        
        if  update_table :
            
            self.data.add_columns(k_array,names=columns_names)
            
            return 
        
        else:           
            
            return Table(k_array,names=columns_names)         
            
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