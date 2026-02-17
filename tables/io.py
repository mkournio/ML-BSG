#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:48:29 2025

@author: michalis
"""

from astropy.io import ascii
import numpy as np

def tab_to_csv(table,filename,**kwargs):
    
    ctab = table.copy().to_pandas()
    ctab.to_csv(filename)
    
    return

def ls_params(output='params',params={},pg_tab=[],meta={},**kwargs):
    
    nterms = params['nterms'].value
    nmod = params['nmod'].value
     
    with open(output, 'w') as file:
        
        if len(meta) > 0: 
            for m in meta:
                file.write(f'#{m} = {meta[m]}\n')
                
        header = ['##STEP']
        header.extend(['FREQ0','P0','SNR','OFFSET'])
        for h in range(1,params['nterms'].value+1):
            header.append(f'AMP{h}')
            header.append(f'PHS{h}')
        header.extend(['W0','R0','TAU','GAMMA',
                       'e_W0','e_R0','e_TAU','e_GAMMA'])
            
        header = ' '.join([f'{x:12s}' for x in header])
        file.write(header+'\n')       
               
        for m in range(nmod+1):
            
            m_par = [m]
            
            if f'f_{m}' in params:
                m_par.extend([params[f'f_{m}'].value,
                     params[f'maxp_{m}'].value,
                     params[f'snr_{m}'].value,
                     params[f'offset_{m}'].value])
                for h in range(1,nterms+1):
                    m_par.append(params[f'ampl_{m}_{h}'].value)
                    m_par.append(params[f'phase_{m}_{h}'].value)
            else:
                m_par.extend(np.full(4 + 2*nterms, np.nan))
                
            if f'W0_{m}' in params:
                m_par.extend([params[f'W0_{m}'].value,
                              params[f'R0_{m}'].value,
                              params[f'TAU_{m}'].value,
                              params[f'GAMMA_{m}'].value,
                              params[f'W0_{m}'].value-params[f'W0_{m}'].min,
                              params[f'R0_{m}'].value-params[f'R0_{m}'].min,
                              params[f'TAU_{m}'].value-params[f'TAU_{m}'].min,
                              params[f'GAMMA_{m}'].value-params[f'GAMMA_{m}'].min
                              ])
            else:
                m_par.extend(np.full(8, np.nan))

            model = ' '.join([f'{x:12.4e}' for x in m_par])
            file.write(model+'\n')
        file.write('###\n')
            
        if len(pg_tab) == 3:
            pg_header = ['##FREQ[c/d]','AMP_INI[mag]','AMP_END[mag]']
            pg_header = ' '.join([f'{x:12s}' for x in pg_header])
            file.write(pg_header+'\n')
            
            for v, p0, p1 in zip(pg_tab[0],pg_tab[1],pg_tab[2]):
                pg_i = [v, p0, p1]
                pg_i = ' '.join([f'{x:12.4e}' for x in pg_i])
                file.write(pg_i+'\n')
        file.write('###\n')            
            
    file.close()
    
    return
