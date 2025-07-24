from tables import *
from methods import *
import args 
import matplotlib.pyplot as plt
import numpy as np
from astroquery.simbad import Simbad
#from evolution import *
import os
from astropy.coordinates import SkyCoord
from constants import *
import pandas as pd
from tess import *


xc = PreProcess(args).compile()
GDW =  [(not any(y in x for y in ('III','V','O9','A0','A1','A2','A3','ON','OC'))) or ('LBV' in x) for x in xc['SpC']]
xc = xc[GDW]

cm = PostProcess(xc).xmatch().append('combined'); cm.sort(['STAR'])
EMS = [('LBV' in x) or ('B[e]SG' in x) for x in cm['SpC']]
cm = cm[EMS]
#cm = cm[(cm['GAL'] == 'MW')]; cm.sort(['STAR'])
#
cm = cm[cm['STAR']=='V439 Cyg']

#cm.sort(['Tmag']); print(cm['STAR','RA','DEC','TIC','Tmag'].pprint(max_lines=-1,max_width=-1))
#LCs = Extract(data=cm, plot_key='dmag', plot_name='ems_tess', output_format='eps')
#LCs.lightcurves(time_bin_size = 0.02, extract_field = True, save_fits = True)
Visualize(data=cm, plot_key='nflux').lightcurves()

'''

#print(cm['STAR','BPmag','RPmag','Jmag','Kmag','BR','JK'].pprint(max_lines=-1))
#mast_query(cm, download_dir='data/', products = ["Lightcurve","ffi"])


LBV = [('LBV' in x) & ('?' not in x) for x in cm['SpC']]
cLBV = [('LBV?' in x) for x in cm['SpC']]
BRC = [('B[e]SG' in x) & ('?' not in x) for x in cm['SpC']]
cBRC = [('B[e]SG?' in x) for x in cm['SpC']]

plt.plot(cm['RA'],cm['DEC'],**plot_all)
plt.plot(cm['RA'][LBV],cm['DEC'][LBV],**plot_LBV)
plt.plot(cm['RA'][cLBV],cm['DEC'][cLBV],**plot_cLBV)
plt.plot(cm['RA'][BRC],cm['DEC'][BRC],**plot_BRC)
plt.plot(cm['RA'][cBRC],cm['DEC'][cBRC],**plot_cBRC)

for g, p in GALAXIES.items(): 
 c = plt.Circle((p['RA'], p['DEC']), p['RAD'], color='r', fill=False)
 plt.gca().add_patch(c)
 plt.text(p['RA'], p['DEC'], g, c='r', fontsize=8)



#print(cm['STAR','SpC','RA','DEC','Gmag','GDIST','MG','SLOGL','REF'].pprint(max_lines=-1))
#plt.plot(cm['TEFF'],cm['SLOGL'],'rx'); plt.xlabel(r'T$_{eff}$'); plt.ylabel(r'log(T$_{eff}^4$/g [$L_{\odot}$])'); plt.xlabel(r'T$_{eff}$ [K]'); plt.gca().invert_xaxis() ; plt.show()

cm = cm[(cm['MK']<-2.5)]; cm.sort(['RA'])
print(cm['STAR','SpC','RA','DEC','REF','DIST','GDIST','MK','TIC'].pprint(max_lines=-1))

LBV = [('LBV' in x) & ('?' not in x) for x in cm['SpC']]
cLBV = [('LBV?' in x) for x in cm['SpC']]
BRC = [('B[e]SG' in x) & ('?' not in x) for x in cm['SpC']]
cBRC = [('B[e]SG?' in x) for x in cm['SpC']]
DR24 = ['D24' in x for x in cm['REF']]

plt.plot(cm['BR'],cm['MG'],**plot_all); plt.xlabel(r'G$_{B}$-G$_{R}$'); plt.ylabel(r'M$_{G}$')  
plt.plot(cm['BR'][DR24],cm['MG'][DR24],'rx'); 
plt.plot(cm['BR'][LBV],cm['MG'][LBV],**plot_LBV) 
plt.plot(cm['BR'][cLBV],cm['MG'][cLBV],**plot_cLBV) 
plt.plot(cm['BR'][BRC],cm['MG'][BRC],**plot_BRC) 
plt.plot(cm['BR'][cBRC],cm['MG'][cBRC],**plot_cBRC)
plt.gca().invert_yaxis() 
plt.show()

'''
'''
plt.plot(cm['JK'],cm['MK'],**plot_all); plt.xlabel(r'J-K'); plt.ylabel(r'M$_{K}$')  
plt.plot(cm['JK'][DR24],cm['MK'][DR24],'rx')
plt.plot(cm['JK'][LBV],cm['MK'][LBV],**plot_LBV)
plt.plot(cm['JK'][cLBV],cm['MK'][cLBV],**plot_cLBV) 
plt.plot(cm['JK'][BRC],cm['MK'][BRC],**plot_BRC)
plt.plot(cm['JK'][cBRC],cm['MK'][cBRC],**plot_cBRC)
plt.gca().invert_yaxis() 

plt.show()
#gal = SkyCoord(ra=np.array(cm['RA'])*u.degree, dec=np.array(cm['DEC'])*u.degree, frame='icrs').galactic
#plt.plot(cm['RUWE'],gal.b,'o')
#plt.plot(cm['RUWE'],cm['GLAT'],'rx'); plt.xlabel(r'RUWE'); plt.ylabel(r'GLAT'); plt.axvline(1.4)   
#plt.show()

#plt.plot(cm['TEFF'],cm['SLOGL'],'rx'); plt.xlabel(r'T$_{eff}$'); plt.ylabel(r'log(T$_{eff}^4$/g [$L_{\odot}$])'); plt.xlabel(r'T$_{eff}$ [K]')  
#plt.plot(np.log10(cm['TEFF']),cm['LOGG'],'rx'); plt.xlabel(r'T$_{eff}$'); plt.ylabel(r'log($g$)'); plt.xlabel(r'log(T$_{eff}$ [K])')



#for tr in _tracks:
 #m_ini = tr.split('p')[0][1:]
 #time, logl, logt, logg, GR_V, GB_V, G_V, V_K, J_K, MV = np.loadtxt(EVOL_TRACKS_PATH + tr, skiprows=2, usecols = (1,3,4,-16, -3,-4,-5,-6, -8, -13), unpack = True)
 
# if 86 > int(m_ini) > 6 : 
#  sll = slogl(10**logt,logg)[1:180]
 # t = 10**logt[1:180]
  #plt.plot(t,sll,'k'); plt.text(t[0],sll[0]-0.04,str(m_ini))
  #plt.plot(logt[0:180],logg[0:180],'k'); plt.text(logt[0],logg[0]-0.02,str(m_ini))
 
 #GB_GR = GB_V - GR_V
 #MK = MV - V_K
 #MG = G_V + MV
 #if 86 > int(m_ini) > 6 : 
 # plt.plot(J_K[1:230],MK[1:230],'k'); plt.text(J_K[0],MK[0]-0.1,str(m_ini)) 
 # plt.plot(GB_GR[1:250],MG[1:250],'k'); plt.text(GB_GR[0],MG[0]-0.1,str(m_ini))
 

#plt.gca().invert_yaxis() 
#plt.show()


'''


