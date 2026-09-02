#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 12:05:11 2026

@author: michalis
"""

import TESS_Localize as tl
import lightkurve as lk
import astropy.units as u 
from methods.functions import *
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

class Localize(object):
    
    def __init__(self,
                 tic,
                 sector,
                 **kwargs):
        
        self.tic = tic 
        self.sector = sector
        
        self.tpf = self._query(**kwargs)
        self.gaia_source, self.gmag = self._get_gaia(**kwargs)
        print('Gaia source %s, G = %s' % (self.gaia_source, self.gmag))
        
        return
    
    def _get_gaia(self,
                  **kwargs):
                 
        result = Vizier(columns=['+_r']).query_region(
            SkyCoord.from_name('TIC %s' % str(self.tic)),
            catalog=['I/355/gaiadr3'], radius= 20 * u.arcsec)[0]
        
        return result['Source'][0], result['Gmag'][0]
    
    def _query(self,
               **kwargs):
        
        search = lk.search_targetpixelfile('TIC %s' % str(self.tic))
        sect_mask = search.mission == 'TESS Sector %s' % self.sector
        search = search[sect_mask]
        
        if any(search.author == 'SPOC'):
            search = search[search.author == 'SPOC']
        else:
            search = search[search.author == 'TESS-SPOC']
            
        tpf = search.download(quality_bitmask='default')  
        
        
        return tpf
    
    def run(self, 
            freq_param,
            pow_lims = [0.01, 0.1, 0.3],
            minf = 0.09,
            pca = 1,
            **kwargs):
        
        run_d = {}
        run_d[f'gmag'] = self.gmag
        run_d[f'pca'] = pca
        run_d[f'minf'] = minf
        
        for pl in pow_lims:
            
            freq_list = freq_indep_sn(freq_param, pow_lim = pl, minf = minf)
            
            tess_loc = tl.Localize(targetpixelfile = self.tpf,
                                   gaia = True,
                                   magnitude_limit = self.gmag + 3.,
                                   frequencies = freq_list,
                                   frequnit = 1 / u.day,
                                   principal_components = pca)
            
            tess_loc.plot_lc(save = f'{self.tic}_s{self.sector}_pl{pl}_pc{pca}.png')
          #  plt.savefig(f'{self.tic}_s{self.sector}_pl{pl}_pc{pca}', format='png')
            
            gaia_list = tess_loc.starfit
            gmask = gaia_list['source'] == str(self.gaia_source)
            gaia_list = gaia_list[gmask]
            
            run_d[f'pval_{pl}'] = gaia_list['pvalue'].to_numpy()
            run_d[f'lhood_{pl}'] = gaia_list['relative likelihood'].to_numpy()
        
        return run_d