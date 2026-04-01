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
from astropy.table import Table

class TimeDomain(object):
    
    def __init__(self,
                 data,
                 measures,
                 **kwargs):
        
        self._validate()
        self.data = data.copy()           
        self.measures = measures
        
        return

    def _validate(self):
        pass
    
    def calculate(self, 
                  bin_size, 
                  stitched = False, 
                  **kwargs):
 
        if len(self.measures) == 0:
            return
        
        if bin_size not in ['raw','10m','30m']:
            raise Exception('Set bin_size among raw, 10m, 30m. Aborting..')
      
        if bin_size == '10m':
            bin_size_d = 0.00694
        elif bin_size == '30m':
            bin_size_d = 0.02083        
       
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                hdus = get_hdu_from_keys(ff[1:], HDUTYPE = 'LIGHTCURVE', BINSIZE = str(bin_size_d))
                print(f'Extracting TD metrics for {star} - {tic}')

                for hdu in hdus:           
                    
                    values = self.td_lc(hdu,self.measures)
                    
                    for m, v in zip(self.measures,values):
                        
                        if m == 'MSE':                            
                            mse = self.lc_mse(hdu)
                            for im, vm in enumerate(mse):
                                hdu.header[f'MSE{im}'] = vm
                        else:
                            if np.isnan(v):
                                v = None
                            hdu.header[m] = v
                        
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                
                ff.close()
                
        return
    
    def header_combine(self, mode ='average', bin_size = 'raw'):
        
        if len(self.measures) == 0:
            return
        
        if not mode in ['average','median']:
            raise Exception('Set mode among average and median. Aborting..')
            
        if bin_size not in ['raw','10m','30m']:
            raise Exception('Set bin_size among raw, 10m, 30m. Aborting..')
            
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083
        
        print('Combining to primary headers the metrics: {}'.format(self.measures))
        
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                sectors = get_sectors_from_hdulist(ff)
                
                if bin_size == 'raw':                    
                    hdus = [get_hdu_from_keys(ff, SECTOR = s,  TTYPE2 = 'flux', BINNING = 'F')[0] for s in sectors]
                else:
                    hdus = [get_hdu_from_keys(ff, SECTOR = s,  TTYPE2 = 'flux', BINSIZE = str(bin_size))[0] for s in sectors]

                s_values = []
                crowdsap_v = []
                
                for hdu, sect in zip(hdus,sectors):
                    
                    sect_values = np.full(len(self.measures), np.nan)                    
                   
                    for i, m in enumerate(self.measures):                        
                        if m in hdu.header and hdu.header[m] != 'None':
                            sect_values[i] = hdu.header[m]  
                        
                    s_values.append(sect_values)
                    crowdsap_v.append(hdu.header['CROWDSAP'])
                        
                if mode == 'average':
                    c_values = np.nanmean(s_values,axis=0)
                    c_values_err = np.nanstd(s_values,axis=0)
                elif mode == 'median':
                    c_values = np.nanmedian(s_values,axis=0)
                    c_values_err = np.nanstd(s_values,axis=0)     
                    
                if isinstance(c_values,np.ndarray) :
                    for m, v, e_v in zip(self.measures, c_values, c_values_err):
                        if np.isnan(v): v = None; e_v=None
                        ff[0].header[m] = (v, e_v)
                    ff[0].header['SCOMBINE'] = mode
                    ff[0].header['SBINSIZE'] = bin_size                   
                    ff[0].header['MINCROWD'] = np.nanmin(crowdsap_v)
                    ff[0].header['AVECROWD'] = np.nanmean(crowdsap_v)
                    
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()                                              
                
        return     
           
    @staticmethod
    def td_lc(f, measures, flux_key = 'dmag'):
        
        lc = get_lc_from_filename(f, flux_key = flux_key)
        lc = lc.remove_outliers(sigma=3.)
        
        mval = np.full(len(measures), np.nan)
        
        for i, m in enumerate(measures):
            
            if m == 'STD':
                mval[i] = get_std(lc.flux)
            elif m == 'MAD':
                mval[i] = get_mad(lc.flux)
            elif m == 'IQR':
                mval[i] = get_iqr(lc.flux)
            elif m == 'SKW':
                mval[i] = get_skew(lc.flux)
            elif m == 'KRT':
                mval[i] = get_kurt(lc.flux)
            elif m == 'ZCR':
                mval[i] = k_cross(lc.flux, kappa = 1)[0]
            elif m == 'ETA':
                mval[i] = get_eta(lc.flux)
            elif m == 'PSI':
                mval[i] = get_psi_sq(lc.flux)
                
        return mval
    
    @staticmethod    
    def lc_mse(f, flux_key = 'dmag'):
        
        lc = get_lc_from_filename(f, flux_key = flux_key)
        lc = lc.remove_outliers(sigma=3.)
        mse = get_mse(lc.flux)
        
        return mse
        
    
class FrequencyDomain(object):
    
    def __init__(self,
                 data,
                 measures,
                 **kwargs):
        
        self.data = data.copy()
        self.measures = measures
        self._validate()
        
        return
    
    def _validate(self):
        
        pass
    
    def calculate(self, 
                  bin_size, 
                  min_freq = 2/27.,
                  stitched = False, 
                  **kwargs):
 
        if len(self.measures) == 0:
            return
        
        if bin_size not in ['raw','10m','30m']:
            raise Exception('Set bin_size among raw, 10m, 30m. Aborting..')
      
        if bin_size == '10m':
            bin_size_d = 0.00694
        elif bin_size == '30m':
            bin_size_d = 0.02083        
       
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                hdus = get_hdu_from_keys(ff[1:], HDUTYPE = 'FREQUENCIES', BINSIZE = str(bin_size_d))
                print(f'Extracting FD metrics for {star} - {tic}')

                for hdu in hdus:
                    
                    values = self.fd_lc(hdu,self.measures, min_freq = min_freq)
                    
                    for m, v in zip(self.measures,values):
                        
                        if m == 'SEN':
                            
                            hdu_per = get_hdu_from_keys(ff[1:], 
                                                        SECTOR = hdu.header['SECTOR'],
                                                        HDUTYPE = 'PERIODOGRAMS', 
                                                        BINSIZE = str(bin_size_d))                            
                            if len(hdu_per) > 0:
                                pg = hdu_per[0].data
                                pg = pg[ pg['frequency'] >= min_freq ]

                                act_m = hdu.data[-1]
                                rn_prop = [act_m['W0'],act_m['R0'],act_m['TAU'],act_m['GAMMA']]
                                snr_spectrum = pg['ampl_ini'] / lorentz(pg['frequency'], *rn_prop)
                                v = get_sen(snr_spectrum**2)
                            
                        if np.isnan(v): 
                            v = None
                            
                        hdu.header[m] = v
                       
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
               
                ff.close()
                
        return
    
    @staticmethod
    def fd_lc(f, measures, min_freq = 2/27., flux_key = 'dmag'):
        
        if isinstance(f, fits.BinTableHDU):
            data = f.data.copy()
        else:
            return
        
        data = data[data['frequency'] >= min_freq]

        mval = np.full(len(measures), np.nan)

        for i, m in enumerate(measures):
            
            if m == 'TOP':
                mval[i] = get_top(data, minf = min_freq)
            elif m == 'HPR':
                mval[i] = get_hpr(data)
            elif m == 'WFM':
                mval[i] = get_wfm(data, minf = min_freq)
            elif m == 'WFD':
                mval[i] = get_wfd(data, minf = min_freq)
                
        return mval 
