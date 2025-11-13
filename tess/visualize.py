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
                      fig_xlabel= PLOT_XLABEL['lc'], 
                      fig_ylabel=PLOT_YLABEL[plot_key], 
                      **kwargs)	 
     
     return
     
    def _validate(self):
        pass
    
    def lightcurves(self, 
                    stitched = False, 
                    bin_size = None, 
                    **kwargs):
        
        if bin_size not in ['10m','30m']:
            raise Exception('Set bin_size among 10m and 30m. Aborting..')
            
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083
            
        ltab = self.data.copy()
        log_file = open("log_vis", "w")
        
        for l in ltab:
            
            star = l['STAR']
            tic = l['TIC']
            spc = l['SpC']
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0:
                
                print('plotting {} TIC {}'.format(star,tic))
                
                hdulist = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                sectors = get_sectors_from_hdulist(hdulist)
                hdu_raw = [get_hdu_from_keys(hdulist, SECTOR = s, BINNING = 'F') for s in sectors]
                hdu_bin = [get_hdu_from_keys(hdulist, SECTOR = s, BINNING = 'T', BINSIZE = str(bin_size)) for s in sectors]
                    
                if stitched:
                    
                    minmax = get_minmax_flux(hdu_bin, flux_key = self.plot_key)
                    grouped_hdu_raw = group_consecutive_hdus(hdu_raw,sectors)
                    
                    ax_scaling = get_ax_scaling(grouped_hdu_raw)
                    axes = self.GridAx(divide=True, ax_scaling = ax_scaling)
                    plot_lc_multi(axes, grouped_hdu_raw, m='.',  flux_key = self.plot_key, lc_type = 'raw')
                    
                    grouped_hdu_bin = group_consecutive_hdus(hdu_bin,sectors)
                    plot_lc_multi(axes, grouped_hdu_bin, flux_key = self.plot_key, lc_type = 'binned')
                        
                    add_plot_features(axes, mode = self.plot_key,
                                      upper_left='{} (TIC {})'.format(star,tic), 
                                    #  upper_right='CROWD {:.2f}'.format(crowdsap),
                                      lower_left=spc, y_min_max = minmax)
                    
                else:
                    
                     for r,b,sect in zip(hdu_raw,hdu_bin,sectors):
                         
                         ax = self.GridAx()
                         
                         plot_lc_single(ax, r, m='.', flux_key = self.plot_key, lc_type = r.header['PIPELINE'])
                         plot_lc_single(ax, b, flux_key = self.plot_key, lc_type = 'binned')
                         
                         add_plot_features(ax, mode = self.plot_key,
                                           upper_left=star, lower_left=spc,
                                           lower_right='{} ({})'.format(tic,sect))
                         
                log_file.write('{:30s} {:+.8f} {:+.8f} {:10s} {} {}\n'.format(star,l['RA'],l['DEC'],spc,tic,get_filename(self.filename,self.output_format)))
                
        self.close_plot()
        log_file.close()
        
        return