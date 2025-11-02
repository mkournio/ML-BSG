import numpy as np
from methods import *

# INITIAL CATALOGS TO COMBILE FROM SORTED BY PREFERENCE IN DOUBLE ENTRIES
INI_CATS = ['EMS','SIMBAD-LBVs','deBURGOS+24','WESMAYER+22','SEARLE+08','McERLEAN+99','CROWTHER+06','FIRNSTEIN+12', 'HAUCKE+19','SIMBAD-BSGs']
        
XM_CATS = {
          'vizier:IV/38/tic' : ['TIC','Tmag'],
          'vizier:I/352/gedr3dis' : ['rpgeo','b_rpgeo','B_rpgeo'],
          'vizier:I/311/hip2' : ['Plx','e_Plx'],
          'vizier:I/355/gaiadr3' : ['BPmag','Gmag','RPmag','e_BPmag','e_Gmag','e_RPmag','RUWE'],
          'vizier:II/246/out' : ['Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag'],
          'vizier:II/328/allwise' : ['W1mag','W2mag','W3mag','W4mag','e_W1mag','e_W2mag','e_W3mag','e_W4mag']
          }
          
#RENAME COLUMN NAMES INTO GLOBAL ONES          
RENAME_COLS = {
          'TEFF' : ['Teff'],
          'e_TEFF' : ['e_Teff'],
          'E_TEFF' : ['E_Teff'], 
          'LOGG' : ['logg'],
          'e_LOGG' : ['e_logg'],  
          'E_LOGG' : ['E_logg'],           
          'RA'   : ['_RAJ2000','RA_d'],
          'DEC'  : ['_DEJ2000','DEC_d'],
          'STAR' : ['Name'],
          'VSINI': ['vsini'],
          'SpC'  : ['SPTYPE','SP_TYPE'],
          'GDIST': ['rpgeo'], 'e_GDIST': ['b_rpgeo'], 'E_GDIST': ['B_rpgeo']
          }
          
APPEND_COLS= {
          'EMS': {**sbcoord_d('RA','DEC', keys=['STAR']), **pow10('TEFF','TEFF2',keys=['LTEFF','LTEFF2'])},
          'FIRNSTEIN+12': {'REF': 'F12', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'SEARLE+08': {'REF': 'S08', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'McERLEAN+99': {'REF': 'M99', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'WESMAYER+22': {'REF': 'W22', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'CROWTHER+06': {'REF': 'C06', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
        #  'FRASER+10' : {'REF': 'F10', 'GAL': 'MW', **coord_h2d('RA','DEC',keys=['RAh','DEC']), 'e_TEFF': 1000, 'e_LOGG': 0.1},
          'HAUCKE+19' : {'REF': 'H19', 'GAL': 'MW', **coord_h2d('RA','DEC',keys=['RAh','DEC'])},
          'deBURGOS+24' : {'REF': 'D24', 'GAL': 'MW', **meancol('e_TEFF', 'e_LOGG', keys=[['e_TEFF','E_TEFF'],['e_LOGG','E_LOGG']])},
          'vizier:I/311/hip2' : {'HDIST': lambda x: 1./(1e-3*x['Plx'])},
          'SIMBAD-BSGs': {'REF': 'SMB'},
          'SIMBAD-LBVs': {'REF': 'SMB', **galloc('GAL',keys=['RA','DEC'])},
          'combined': {
                      **galcoord('GLON','GLAT',keys=['RA','DEC']),
                      **dist('DIST',keys=['GAL','GDIST']),
                      **absmag('MJ','MH','MK','MG',keys=['Jmag','Hmag','Kmag','Gmag','DIST']),
                      **diffcol('BR','JK',keys=[['BPmag','RPmag'],['Jmag','Kmag']]),
                      **slogl('SLOGL',keys=['TEFF','LOGG']),
                      **spc2t('SpCt',keys=['SpC'])
                      }
          }	
          
   

#KEY COLUMN FOR CHECKING DOUBLE ENTRIES
DOUBLE_KEYS = ['STAR']

MERGE_KEYS = ['REF']


