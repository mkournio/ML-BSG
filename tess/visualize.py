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
    
    # Class for 

    def __init__(self, data, rows_page = 5, **kwargs):

     self.data = data
   #  self._validate()

     super().__init__(rows_page = PLOT_XLC_NROW, cols_page = PLOT_XLC_NCOL,
                      fig_xlabel= PLOT_XLABEL['lc'], fig_ylabel='', **kwargs)	 
     #self._extract_LCs(**kwargs)
     
  
    def lightcurves(self, **kwargs):
       
       lc_files = os.listdir(path_to_fits)       
       stars = self.data['STAR']
       tics = self.data['TIC']
       spcs = self.data['SpC']
       ras = self.data['RA']
       decs = self.data['DEC']

       for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
           
        filename = [f for f in lc_files if star in f]
        if len(filename) > 0:
            hdulist = fits.open(os.path.join(path_to_fits,filename[0]))
            sectors = get_sectors_from_hdulist(hdulist)
            for s in sectors:
                hdu = get_hdu_from_keys(hdulist, SECTOR = s, TTYPE2 = 'DMAG')
                ax = self.GridAx()
                plot_lc_single(hdu, ax=ax, flux_key ="dmag", lc_type = 'spoc_binned')
        else:
            print('Star {} not found in database'.format(star))
            
       self._page_close()
       
       return

        
       '''
        
        if len(files) > 0 : 
            
            s_sects,s_files = zip(*sorted(zip(sects,files)))
            
            filename = os.path.join(path_to_fits,'{}_{}.fits'.format(star,tic))
            if kwargs.get('save_fits'):
             ff = FitsObject(filename)

            for f, sect in zip(s_files,s_sects):
           
             print('{}: extracting Sector {} of TIC {} (SPOC)'.format(star,sect,tic))
  
             fits = os.listdir(path_to_spoc_files + f)[0]
             lc = lk.TessLightCurveFile(os.path.join(path_to_spoc_files + f,fits)).remove_nans().remove_outliers()
             if kwargs.get('save_fits'): 
                 ff.add_lc(lc, header_source = lc, flux_key = 'flux', binning = 'RAW')
            
             ax_lc = self.GridAx()
             
             plot_lc_single(lc, ax=ax_lc, lc_type = 'spoc', m = '.')
             if time_bin_size is not None:
               lc = lc.bin(time_bin_size = time_bin_size)
               plot_lc_single(lc, ax=ax_lc, lc_type = 'spoc_binned')
               '''
