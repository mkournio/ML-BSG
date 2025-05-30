from astroquery.simbad import Simbad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii

GALAXIES = {
           'LMC': [80, -69, 11],
           'SMC': [21, -73, 13],
           'M31': [10.4, 41, 2],
           'M33': [23.4, 30.6, 1],
           'NGC3109': [150.75, -26.15, 1],
           'NGC6822': [296.2, -14.8, 1],
           'NGC2403': [114.2, 65.6, 1],
           'NGC300': [13.7, -37.70, 1],
           'NGC55': [3.8, -39.2, 1],
           'IC1613': [16.23, 2.13, 1],
           'WLM': [0.5, -15.5, 1],
           'SEXTA': [152.7, -4.7, 1],
           'LEOA': [149.9, 30.7, 1],
           'M101': [210.8, 54.3, 1]
           }
       
def find_id(x):
 x = x.split('|')
 for c in ['HD','HR','BD','MCW','HIP','Hen','LS','ALS','2MASS','TIC']:
   for i in x :
     i=''.join(i.split())
     if i.startswith(c): 
      return i 
 return i 

def find_gal(ra,dec):
 for g, p in GALAXIES.items(): 
  if ((ra-p[0])**2) + ((dec-p[1])**2) < (p[2]**2):  
    return g    
 return 'MW'

'''
Simbad.ROW_LIMIT = -1
Simbad.add_votable_fields('ids','sptype','rv_value','v*','flux(V)','ra(d)', 'dec(d)','bibcodelist(1850-2023)')
script = '((sptype > O9.5 & sptype < A0) | (sptypes in ("B1/2","B2/3","B0/1","B5/7","B0.5/1"))) & sptypes in ("Ia","Ib","I0","0","Iab","II","Ib/II","Iab-Ib","Ia/ab","Ib/II","Iab/b","Ia/ab","Ib-II","I/II","Iap","Ibp","II:")'
tab = Simbad.query_criteria(script)

tab['STAR'] = [find_id(x) for x in tab['IDS']]
tab['GAL'] = [find_gal(x,y) for x,y in zip(tab['RA_d'],tab['DEC_d'])]

tab = tab['GAL','STAR','RA_d', 'DEC_d','SP_TYPE','RV_VALUE','FLUX_V','V__vartyp','BIBLIST_1850_2023','MAIN_ID']
tab.sort(['GAL','RA_d','DEC_d'])
tab_pd = pd.DataFrame(np.array(tab))
#tab_pd.to_csv('SIMBAD-BSGs2', sep=',', index=False)
'''
tab = ascii.read('SIMBAD-BSGs')

plt.plot(tab['RA_d'],tab['DEC_d'],'.'); plt.xlabel('RA (deg)'); plt.ylabel('DEC (deg)')
for g, p in GALAXIES.items(): 
 c = plt.Circle((p[0], p[1]), p[2], color='r', fill=False)
 plt.gca().add_patch(c)
 plt.text(p[0], p[1], g, c='r', fontsize=8)
plt.show()
