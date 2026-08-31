#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 12:05:11 2026

@author: michalis
"""

import TESS_Localize as tl
import lightkurve as lk
import astropy.units as u 
    
class Localize(object):
    
    def __init__(self,
                 tic,
                 **kwargs):
        
        self.tic = tic 
        
        return
    
    def _query(self,
               sector,
               ):
        
        search = lk.search_targetpixelfile('TIC %s' % str(self.tic))
        sect_mask = search.mission == 'TESS Sector %s' % sector
        search = search[sect_mask]
        
        if any(search.author == 'SPOC'):
            search = search[search.author == 'SPOC']
        else:
            search = search[search.author == 'TESS-SPOC']
            
        tpf = search.download(quality_bitmask='default')
        
        return tpf
    
    def run(self, 
            freq_list,
            mag_limit = 10,
            pca = 3,
            **kwargs):
        
        tpf = self._query(**kwargs)       
        gaia_list = tl.Localize(targetpixelfile = tpf, 
                              gaia = True, 
                              magnitude_limit = mag_limit,
                              frequencies = freq_list, 
                              frequnit = 1 / u.day, 
                              principal_components = pca)
        
        return gaia_list