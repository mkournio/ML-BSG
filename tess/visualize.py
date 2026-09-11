#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 23:49:43 2025

@author: michalis
"""
import lightkurve as lk
import numpy as np
import os
from constants import *
from methods.functions import *
from methods.plot import *
from methods.tools import *
from astropy.io import fits

class Visualize(GridTemplate):
    
    # Class for visualizing objects
    
    def __init__(self, 
                 data, 
                 plot_key = 'flux',
                 rows_page = PLOT_XLC_NROW,
                 cols_page = PLOT_XLC_NCOL,
                 **kwargs):
        
     self._validate()
     self.data = data
     self.plot_key = plot_key
     
     super().__init__(rows_page = rows_page, 
                      cols_page = cols_page,
                      fig_xlabel = PLOT_XLABEL[plot_key], 
                      fig_ylabel = PLOT_YLABEL[plot_key], 
                      **kwargs) 
     
     return
     
    def _validate(self):
        pass
    
    def lightcurves(self, 
                    stitched = False, 
                    models = False,
                    bin_size = '10m', 
                    lc_props = None,
                    preview = None,
                    dict_log = {},
                    **kwargs):
        
        if bin_size not in ['10m','30m']:
            raise Exception('Set bin_size among 10m and 30m. Aborting..')
            
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083
            
        #if 'output_path' in kwargs:
         #   path_to_output_fits = kwargs['output_path']
            
        ltab = self.data.copy()
       # log_file = open("log_vis_ems", "w")
        
        for l in ltab:
            
            star = l['STAR']
            tic = l['TIC']
            try:
                spc = l['SpC']
            except:
                spc = ''
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0:
                
                print('LC plotting {} TIC {}'.format(star,tic))
                
                hdulist = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                sectors = get_sectors_from_hdulist(hdulist)
                if 'custom_sect' in kwargs:
                    if star in kwargs['custom_sect']:
                        sectors = kwargs['custom_sect'][star]
             
                hdu_raw = [get_hdu_from_keys(hdulist, SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINNING = 'F')[0] for s in sectors]
                hdu_bin = [get_hdu_from_keys(hdulist, SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINNING = 'T', BINSIZE = str(bin_size))[0] for s in sectors]
                if models:
                    hdu_mods = [get_hdu_from_keys(hdulist, SECTOR = s, HDUTYPE = 'FREQUENCIES', BINNING = 'T', BINSIZE = str(bin_size))[0] for s in sectors]
                    grouped_hdu_mods = group_consecutive_hdus(hdu_mods,sectors)
                    
                if 'B[e]SG' in spc:
                    lc_type = 'B[e]SG'
                elif 'LBV' in spc:
                    lc_type = 'LBV'
                else:
                    lc_type = 'raw'
                    
                if stitched:
                    
                    minmax = get_minmax_flux(hdu_bin, flux_key = self.plot_key)
                    grouped_hdu_raw = group_consecutive_hdus(hdu_raw,sectors)
                    
                    ax_scaling = get_ax_scaling(grouped_hdu_raw)
                    axes = self.GridAx(divide=True, ax_scaling = ax_scaling)
                    plot_lc_multi(axes, grouped_hdu_raw, m='.',  flux_key = self.plot_key, lc_type = 'raw')
                    
                    grouped_hdu_bin = group_consecutive_hdus(hdu_bin,sectors)
                    plot_lc_multi(axes, grouped_hdu_bin, flux_key = self.plot_key, lc_type = 'binned')
                    
                    if models and self.plot_key == 'dmag':
                        plot_mod_multi(axes, grouped_hdu_mods, ref_hdus = grouped_hdu_bin, m='-', lw=0.6)
                     
                    vlines = []   
                    for obj in dict_log:
                        if obj == star:
                            vlines = dict_log[obj]
                            
                    add_plot_features(axes, mode = self.plot_key,
                                      upper_left='{} (TIC {})'.format(star,tic), 
                                    #  upper_right='CROWD {:.2f}'.format(crowdsap),
                                      lower_left=spc, y_min_max = minmax, vlines = vlines)
                    
                else:
                     
                     for i in range(len(hdu_raw)):
                         
                         r = hdu_raw[i]
                         b = hdu_bin[i]
                         sect = sectors[i]                           
                                                                                               
                         if preview != None:
                             pr_ind = 0
                             with open(preview) as pv:
                                 for l in pv:
                                     l = l.split()
                                     if l[0] == self.plot_key and int(l[1]) == tic and int(l[2]) == sect:                                         
                                         pr_ind = 1
                                         break
                             if pr_ind == 0:
                                 continue
                         
                         ax = self.GridAx()
                         
                         prop_args = []
                         if lc_props != None:
                             with open(lc_props) as pf:
                                 for l in pf:
                                     l = l.split()
                                     if int(l[0]) == tic and int(l[1]) == sect:
                                             prop_args = [float(x) for x in l[2:]]

                         
                         plot_lc_single(ax, r, m='.', flux_key = self.plot_key, lc_type = lc_type,**kwargs)
                         plot_lc_single(ax, b, flux_key = self.plot_key, lc_type = 'binned', prop_args = prop_args,**kwargs)
                         if models and self.plot_key == 'dmag':
                             mod_hdu = hdu_mods[i]
                             plot_mod_single(ax, mod_hdu, ref_hdu = b, ls='--', lw=1.5)
                  
                         add_plot_features(ax, mode = self.plot_key,
                                           upper_left=st(star)[0], lower_left=spc,
                                           lower_right='{} ({})'.format(tic,sect))
                         
             #   log_file.write('{:30s} {:+.8f} {:+.8f} {:10s} {} {}\n'.format(star,l['RA'],l['DEC'],spc,tic,get_filename(self.filename,self.output_format)))
 
        self.close_plot()
     #   log_file.close()
        
        return
    
    def periodograms(self,
                     bin_size = '10m',
                     preview = None,
                     **kwargs
                     ):
        
        if bin_size not in ['10m','30m']:
            raise Exception('Set bin_size among 10m and 30m. Aborting..')
            
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083       
    
        ltab = self.data.copy()
        
        for l in ltab:
            
            star = l['STAR']
            tic = l['TIC']
            spc = l['SpC']
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0:
                
                print('LS plotting {} TIC {}'.format(star,tic))
                
                hdulist = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                sectors = get_sectors_from_hdulist(hdulist)
                hdu_pgs = [get_hdu_from_keys(hdulist, SECTOR = s, HDUTYPE = 'PERIODOGRAMS', BINNING = 'T', BINSIZE = str(bin_size))[0] for s in sectors]
                hdu_rns = [get_hdu_from_keys(hdulist, SECTOR = s, HDUTYPE = 'FREQUENCIES', BINNING = 'T', BINSIZE = str(bin_size))[0] for s in sectors]
                
                if 'B[e]SG' in spc:
                    c_class = 'B[e]SG'
                elif 'LBV' in spc:
                    c_class = 'LBV'
                else:
                    c_class = 'any'
                
                for i in range(len(hdu_pgs)):
                      
                      pg = hdu_pgs[i]
                      rn = hdu_rns[i]
                      sect = sectors[i]
                      
                      if preview != None:
                          pr_ind = 0
                          with open(preview) as pv:
                              for l in pv:
                                  l = l.split()
                                  if l[0] == self.plot_key and int(l[1]) == tic and int(l[2]) == sect:
                                      pr_ind = 1
                                      break
                              if pr_ind == 0:
                                  continue                      
                      
                      ax = self.GridAx()
                      
                      plot_ls_single(ax, pg, model = rn, c_class=c_class, **kwargs)
                      
                      add_plot_features(ax, mode = self.plot_key,
                                        upper_left=st(star)[0], lower_left=spc,
                                        upper_right='{} ({})'.format(tic,sect))
                
        self.close_plot()
        
        return
                

   
    
    
    
    