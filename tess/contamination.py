#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 12:05:11 2026

@author: michalis
"""

import TESS_Localize as tl
import astropy.units as u 
from methods.functions import *
from astropy.coordinates import SkyCoord
from constants import paths
import astropy.units as u
import matplotlib.pyplot as plt
from methods.tools import *
import copy
import os

class Contamination(object):
    
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
                  bin_size = '10m', 
                  mode = 'crowdsap',
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
                sectors = get_sectors_from_hdulist(ff)
                
                print(f'Contamination: Running {mode} for {star} - {tic}')
                
                if mode == 'crowdsap':
                    
                    hdus = [get_hdu_from_keys(ff[1:], SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINSIZE = str(bin_size_d))[0] for s in sectors]
                    for hdu in hdus:
                        
                        sec = hdu.header['SECTOR']
                        tic = hdu.header['TICID']
                        CR = Crowdsap(tic, sec, query_service = 'mast', **kwargs)#.run()
                    
                
                
                elif mode == 'localize':
                    
                    hdus_raw = [get_hdu_from_keys(ff[1:], SECTOR = s, HDUTYPE = 'LIGHTCURVE', BINNING = 'F')[0] for s in sectors]
                    hdus = [get_hdu_from_keys(ff[1:], SECTOR = s, HDUTYPE = 'FREQUENCIES', BINSIZE = str(bin_size_d))[0] for s in sectors]
                    
                    for hdu, hdu_r in zip(hdus,hdus_raw):
                        
                        sec = hdu.header['SECTOR']
                        tic = hdu.header['TICID']
                        print()
                        trend = lk.LightCurve(time = hdu_r.data['time'], flux = hdu_r.data['trend'])
                        
                        TL = Localize(tic,sec, trend=trend)
                        TL.run(hdu.data)
                        tld = TL.run_d; print(tld)
                        
                        for m in self.measures:
                            
                            if m == 'TLPW':
                                hdu.header[m] = tld['pval_0.01'][0]
                            if m == 'TLLW':
                                hdu.header[m] = tld['lhood_0.01'][0]
                            if m == 'TLPS':
                                hdu.header[m] = tld['pval_0.1'][0]
                            if m == 'TLLS':
                                hdu.header[m] = tld['lhood_0.1'][0]   
                            
                #ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                
                ff.close()
                
        return
    
class Crowdsap(object):
    
    def __init__(self,
                 tic,
                 sector,
                 **kwargs):
        
        self.tic = tic
        self.sector = sector
        #self.tpf = query_tpf(tic,sector,**kwargs)
        result = get_gaia(tic).to_pandas()
        #mask = result['RPmag'] < 18
        #result = result[mask]
        
        g_target = result.iloc[0]
        
        for r in [2,3,4,5]:
            mask = result['_r'] < r * TESS_pix_size
            mresult = result[mask]
            
            RPdiff = mresult['RPmag'] - min(mresult['RPmag'])
            RPdiff = 10 ** (-0.4 * RPdiff)
            
            print(sector,1 / np.nansum(RPdiff))
        #print(np.nansum(RPdiff))
        
        
        
        return
    
    def run(self,
            r = 10,
            gaia_cat = 'I/355/gaiadr3',
            **kwargs):
        
        
        
        '''
        mask = np.zeros(self.tpf[0].shape[1:], dtype='bool')
        
        for i in range(mask.shape[0]):
            for j in range(nmask.shape[1]):
                if nmask[i,j]:
                    tpf_row = ref_p[0] + i
                    tpf_col = ref_p[1] + j
                    cond_box = (self.RA_pix>tpf_col) & (self.RA_pix<tpf_col+1) & (self.DE_pix>tpf_row) & (self.DE_pix<tpf_row+1)
                    
                    loc_diff = self.RPdiff[cond_box & (self.RPdiff != 0)]
                    loc_diff = loc_diff[~np.isnan(loc_diff)]
                    if len(loc_diff) > 0 and min(loc_diff) < min_diff :
                        min_diff = min(loc_diff)
                        if update and min_diff < min_thres :
                               for di in [max(0,i-1),i,min(i+1,nmask.shape[0]-1)]:
                                   for dj in [max(0,j-1),j,min(j+1,nmask.shape[1]-1)]:
                                       nmask[di,dj] = False

        #print(mask)
        '''
        
        return
    
    
    
    
    
class Localize(object):
    
    def __init__(self,
                 tic,
                 sector,
                 trend = None,
                 **kwargs):
        
        self.tic = tic 
        self.sector = sector

        self.tpf = self._query(**kwargs)        
        self.tpf = self._corr_tpf(trend)
        
        #sap_corr = self.tpf.to_lightcurve(aperture_mask=self.tpf.pipeline_mask)
        #sap_corr.plot()
        #plt.show()
        #dmag_from_tpf = -2.5 * np.log10(sap_corr.flux.value / np.nanmedian(sap_corr.flux.value))
        #plt.plot(sap_corr.time.value,dmag_from_tpf); plt.gca().invert_yaxis() 

        self.gaia_source, self.gmag = self._get_gaia(**kwargs)
        print(f'TIC {tic} s{sector} - Gaia {self.gaia_source}, G = {self.gmag}')
        
        self.run_d = {}
        self.run_d['tic'] = self.tic
        self.run_d['sector'] = self.sector  
        self.run_d['gsource'] = self.gaia_source
        self.run_d['gmag'] = self.gmag
        
        return
    
    def _corr_tpf(self, 
                  trend):
        
        if isinstance(trend, lk.LightCurve):
            
            tpf_corr = copy.deepcopy(self.tpf)            
            qmask = np.asarray(tpf_corr.quality_mask, dtype=bool)
            
            t_tpf = np.asarray(tpf_corr.time.value, dtype=float)
            t_trend = np.asarray(trend.time.value, dtype=float)
            trend_flux = np.asarray(trend.flux.value, dtype=float)

            mask = np.isfinite(t_trend) & np.isfinite(trend_flux)
            
            t_trend = t_trend[mask]
            trend_flux = trend_flux[mask]            
            trend_flux = trend_flux / np.nanmedian(trend_flux)
            
            trend_interp = np.interp(t_tpf, t_trend, trend_flux, left = np.nan, right = np.nan)
            
            flux = np.asarray(tpf_corr.hdu[1].data["FLUX"], dtype=float)
            flux_err = np.asarray(tpf_corr.hdu[1].data["FLUX_ERR"], dtype=float)
            
            flux[qmask] = flux[qmask] / trend_interp[:, None, None]
            flux_err[qmask] = flux_err[qmask] / trend_interp[:, None, None]
            
            tpf_corr.hdu[1].data["FLUX"][:] = flux
            tpf_corr.hdu[1].data["FLUX_ERR"][:] = flux_err 
            
            sap_corr = tpf_corr.to_lightcurve(aperture_mask=tpf_corr.pipeline_mask)
            sap_clean = sap_corr.remove_outliers()
            
            keep = np.isin(
                np.asarray(tpf_corr.time.value, dtype=float),
                np.asarray(sap_clean.time.value, dtype=float))
            
            return tpf_corr[keep] 
        
        else:
            
            return self.tpf  
        
        
    
    def _get_gaia(self,
                  **kwargs):
                 
        result = Vizier(columns=['+_r']).query_region(
            SkyCoord.from_name('TIC %s' % str(self.tic)),
            catalog=['I/355/gaiadr3'], radius= 20 * u.arcsec)[0]
        
        return result['Source'][0], result['Gmag'][0]
    
    def _query(self,
               **kwargs):
        
        search = lk.search_targetpixelfile('TIC %s' % str(self.tic))
        sect_mask = search.mission == 'TESS Sector %02d' % self.sector
        search = search[sect_mask]        
        
        if any(search.author == 'SPOC'):
            search = search[search.author == 'SPOC']
        else:
            search = search[search.author == 'TESS-SPOC']
            
        tpf = search.download(quality_bitmask='default') 
        
        return tpf
    
    def run(self, 
            freq_param,
            pow_lims = [1e-2, 0.1],
            minf = 0.06,
            pca = 1,
            **kwargs):
        
        if isinstance(pow_lims,float):
            pow_low, pow_up = [pow_lims], [np.inf]
        else:
            pow_low = pow_lims
            pow_up = np.append(pow_lims[1:],np.inf)            
       
        self.run_d[f'pca'] = pca
        self.run_d[f'minf'] = minf   
        
        tl_kwargs ={
            'targetpixelfile': self.tpf,
            'gaia' : True,
            'magnitude_limit': self.gmag + 3.,
            'frequnit': 1 / u.day,
            'principal_components': pca            
            }
        
        for pl, pu in zip(pow_low,pow_up):
            
            freq_list = freq_indep_sn(freq_param, pow_lim = [pl,pu], minf = minf)
            print(pl, pu, freq_list)
            
            if len(freq_list) > 0:
                tl_kwargs['frequencies'] = freq_list                 
            else:
                self.run_d[f'pval_{pl}'] = [None]
                self.run_d[f'lhood_{pl}'] = [None]
                
                continue
                
            try:
                tess_loc = tl.Localize(**tl_kwargs)
            except:
                tl_kwargs['principal_components'] = 0                
                tess_loc = tl.Localize(**tl_kwargs)
                print('Moved to PCA = 0..')
            
            plot_path = path_to_output_localize + f'{self.tic}_s{self.sector}_pl{pl}_pc{pca}'
           
            tess_loc.plot_lc(save = f'{plot_path}_LC.png')
            tess_loc.plot(method='snr', save = f'{plot_path}_SNR.png')
            plt.close('all')
            
            gaia_list = tess_loc.starfit
            gmask = gaia_list['source'] == str(self.gaia_source)
            gaia_list = gaia_list[gmask]
            
            self.run_d[f'pval_{pl}'] = gaia_list['pvalue'].to_numpy()
            self.run_d[f'lhood_{pl}'] = gaia_list['relative likelihood'].to_numpy()
        
        return