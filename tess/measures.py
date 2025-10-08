#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 11:13:10 2025

@author: michalis
"""

import lightkurve as lk
import numpy as np
from methods.tools import *
from astropy.io import fits
from scipy import stats
from methods.functions import *


class Measures(object):
    
    def __init__(self,
                 data, 
                 **kwargs):
        
        self._validate()
        self.data = data
        
    def _validate(self):
        pass
    
    def calculate(self, *measures, binned = True, stitched = False, **kwargs):
        
        if len(measures) == 0:
            return
        
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                sectors = get_sectors_from_hdulist(ff)
                
                if binned:
                    hdus = [get_hdu_from_keys(ff, SECTOR = s, BINNING = 'T') for s in sectors]
                else:
                    hdus = [get_hdu_from_keys(ff, SECTOR = s, BINNING = 'F') for s in sectors]
       
                for hdu, sect in zip(hdus,sectors):
                    
                    values = self._lc_measures(hdu,*measures)
                    for m, v in zip(measures,values):
                        hdu.header[m] = v
                        
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()
                
        return
    
    def combine_from_headers(self, *measures, mode='avg'):
        
        return
    
    def _lc_measures(self, lc, *measures, flux_key = 'dmag'):

        mval = np.full(len(measures), np.nan)

        if isinstance(lc, lk.LightCurve):
            time = lc.time.value
            flux = lc[flux_key].value
            try:
                flux_err = lc[flux_key+"_err"].value
            except:
                flux_err = np.zeros(len(flux))                
        elif isinstance(lc, fits.BinTableHDU):
            
            time = lc.data['time']
            flux = lc.data[flux_key]
            try:
                flux_err = lc.data[flux_key+"_err"]
            except:
                flux_err = np.zeros(len(flux))                 
        elif isinstance(lc, np.ndarray):
            time = lc[0]
            flux = lc[1]
            if lc.shape[0] > 2:
                flux_err = lc[2]
            else:
                flux_err = np.zeros(len(flux))                
        else:
            raise TypeError('Object light curve does not have supportive format!')
            
        nan_mask = np.isnan(flux)
        flux = flux[~nan_mask]

        for i, m in enumerate(measures):
            
            if m == 'std':
                mval[i] = get_std(flux)
            elif m == 'mad':
                mval[i] = get_mad(flux)
            elif m == 'skw':
                mval[i] = get_skew(flux)
            elif m == 'eta':
                mval[i] = get_eta(flux)
            elif m == 'psi':
                mval[i] = get_psi_sq(flux)         
            
                
        return mval 