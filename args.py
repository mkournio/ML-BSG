import numpy as np
from functions import *

# INITIAL CATALOGS TO COMBILE FROM SORTED BY PREFERENCE IN DOUBLE ENTRIES
INI_CATS = ['EMS']#,'vizier:J/A+A/687/A228/tabled1','WESMAYER+22','SEARLE+08','McERLEAN+99','CROWTHER+06','FIRNSTEIN+12', 'FRASER+10', 'HAUCKE+19','SIMBAD-BSGs']
        
XM_CATS = {
         # 'vizier:IV/38/tic' : ['TIC','Tmag'],
          #'vizier:I/352/gedr3dis' : ['rpgeo','b_rpgeo','B_rpgeo'],
         # 'vizier:I/311/hip2' : ['Plx','e_Plx'],
        #  'vizier:I/355/gaiadr3' : ['BPmag','Gmag','RPmag','e_BPmag','e_Gmag','e_RPmag','RUWE'],
         # 'vizier:II/246/out' : ['Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag'],
          #'vizier:II/328/allwise' : ['W1mag','W2mag','W3mag','W4mag','e_W1mag','e_W2mag','e_W3mag','e_W4mag']
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
          'EMS': {**sbcoord_d('RA','DEC',keys=['STAR']), **pow10('TEFF','TEFF2',keys=['LTEFF','LTEFF2'])},
          'FIRNSTEIN+12': {'REF': 'F12', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'SEARLE+08': {'REF': 'S08', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'McERLEAN+99': {'REF': 'M99', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'WESMAYER+22': {'REF': 'W22', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'CROWTHER+06': {'REF': 'C06', 'GAL': 'MW', **sbcoord_d('RA','DEC',keys=['STAR'])},
          'FRASER+10' : {'REF': 'F10', 'GAL': 'MW', 'e_TEFF': 1000, 'E_TEFF': 1000, 'e_LOGG': 0.1, 'E_LOGG': 0.1},
          'HAUCKE+19' : {'REF': 'H19', 'GAL': 'MW', 'E_TEFF': lambda x: x['e_TEFF'], 'E_LOGG': lambda x: x['e_LOGG']},
          'vizier:J/A+A/687/A228/tabled1' : {'REF': 'D24', 'GAL': 'MW'},
          'vizier:I/311/hip2' : {'HDIST': lambda x: 1./(1e-3*x['Plx'])},
          'SIMBAD-BSGs': {'REF': 'SMB'},
          'combined': {
                      **galcoord('GLON','GLAT',keys=['RA','DEC']),
                   #   'MG' : lambda x: x['Gmag'] - (5 * np.log10(x['GDIST'])) + 5,
                    #  'MK' : lambda x: x['Kmag'] - (5 * np.log10(x['GDIST'])) + 5,
                     # 'BR' : lambda x: x['BPmag']-x['RPmag'],
                      #'JK' : lambda x: x['Jmag']-x['Kmag'],
                      **slogl('SLOGL',keys=['TEFF','LOGG'])
                      }
          }	
          
   

#KEY COLUMN FOR CHECKING DOUBLE ENTRIES
DOUBLE_KEYS = ['STAR']

MERGE_KEYS = ['REF']


