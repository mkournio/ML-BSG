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
        
        if 'MSE' in measures:
            measures.remove('MSE')
            measures += ['MSP','MSD','MSC','MSS']  
            
        self.measures = measures
        
        return

    def _validate(self):
        pass
    
    def calculate(self, stitched = False, **kwargs):
 
        if len(self.measures) == 0:
            return
        
        print('Calculating the time-domain metrics: {}'.format(self.measures))
        
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                for hdu in ff[1:]:
                    
                    values = self.single_lc(hdu,self.measures)
                    
                    for m, v in zip(self.measures,values):
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
    
    def reset_headers(self, primary = True, sectors = False):
        
        if (len(self.measures) == 0) or (not any([primary,sectors])):
            
            return 
        
        print('Reseting time-domain metrics..')
        
        extra_keys = ['SCOMBINE','SBINSIZE']
        measures = self.measures + extra_keys 
        
        stars = self.data['STAR']
        tics = self.data['TIC']        
        for star, tic in zip(stars,tics):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                hdus = ff
                if primary and not sectors:
                    hdus = [ff[0]]
                elif sectors and not primary:
                    hdus = ff[1:]
                
                for hdu in hdus:
                        hdr = hdu.header
                        for m in measures:
                            if m in hdr: del hdr[m]
                    
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()
                   
        return        
           
    @staticmethod
    def single_lc(f, measures, flux_key = 'dmag'):
        
        lc = f.copy()
        
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
        flux = flux[20:-20]
        
        mval = np.full(len(measures), np.nan)
        
        if 'MSP' in measures:
            msp, msd, msc, mss = get_mse(flux)

        for i, m in enumerate(measures):
            
            if m == 'STD':
                mval[i] = get_std(flux)
            elif m == 'MAD':
                mval[i] = get_mad(flux)
            elif m == 'IQR':
                mval[i] = get_iqr(flux)
            elif m == 'SKW':
                mval[i] = get_skew(flux)
            elif m == 'KRT':
                mval[i] = get_kurt(flux)
            elif m == 'ZCR':
                mval[i] = k_cross(flux, kappa = 1)[0]
            elif m == 'ETA':
                mval[i] = get_eta(flux)
            elif m == 'PSI':
                mval[i] = get_psi_sq(flux)
            elif m == 'MSP':
                mval[i] = msp
            elif m == 'MSD':
                mval[i] = msd
            elif m == 'MSC':
                mval[i] = msc
            elif m == 'MSS':
                mval[i] = mss
                
        return mval 
    
class FrequencyDomain(object):
    
    def __init__(self,
                 data,
                 **kwargs):
        
        self.data = data.copy()
        self.validate()
        
        return
    
    def _validate(self):
        
        pass