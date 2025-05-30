from tables import *
from functions import *
import args 
import matplotlib.pyplot as plt
import numpy as np
from astroquery.simbad import Simbad
#from evolution import *
import os
from astropy.coordinates import SkyCoord

#EM = EvolutionaryModels()
#data['MEVOL'] = EM.interp2d('logTe', 'logg', 'Mass', data['TEFF'], data['LOGG'], post_rsg = False, method='linear')

###### SPECTROSCOPIC HR DIAGRAM
#EM.plot_spectroHR(data, output_format = 'pdf', output_name = 's_hrdiag', hold = False, inter = False)
#x_range = np.arange(3.6,4.9,0.1)
#f = lambda x, gamma : 4 * x - np.ma.log10(gamma / (KAPPA * SBOLTZ / Ccm))
#g1, = EM.ax.plot(x_range,f(x_range, 1),'r-', lw = 2)
#g06, = EM.ax.plot(x_range,f(x_range, 0.6),'r:', lw = 2)
#EM.ax.legend([g1, g06], [ r'$\Gamma = 1$',r'$\Gamma = 0.6$'])444
#EM.panel.PanelClose()

EVOL_TRACKS_PATH = 'input/tracks/r04/'
dict_tr = {}
_tracks = os.listdir(EVOL_TRACKS_PATH)

xc = PreProcess(args).compile()
GDW =  [not any(y in x for y in ('III','V','O9','A0','A1','A2','A3','ON','OC')) for x in xc['SpC']]
xc = xc[GDW]; 

cm = PostProcess(xc).xmatch()#.append('combined')
print(cm)
#print(cm['STAR','RA','DEC','SLOGL','GLON','GLAT'])
#LBV = ['lbv' in x for x in cm['SpC']]
#BRC = ['B[e]SG' in x for x in cm['SpC']]

#plt.plot(cm['RA'],cm['DEC'],'.')
#plt.plot(cm['RA'][LBV],cm['DEC'][LBV],'ko',mew=2,ms=14,mfc='None');
#plt.plot(cm['RA'][BRC],cm['DEC'][BRC],'ks',mew=2,ms=14,mfc='None'); 
#plt.show()
#print(cm['STAR','SpC','RA','DEC','Gmag','GDIST','MG','SLOGL','REF'].pprint(max_lines=-1))
#plt.plot(cm['TEFF'],cm['SLOGL'],'rx'); plt.xlabel(r'T$_{eff}$'); plt.ylabel(r'log(T$_{eff}^4$/g [$L_{\odot}$])'); plt.xlabel(r'T$_{eff}$ [K]'); plt.gca().invert_xaxis() ; plt.show()
'''
cm = cm[(cm['GAL']=='MW') & (cm['MK']<-2.5)]; cm.sort(['RA'])
print(cm['STAR','SpC','RA','DEC','Gmag','GDIST','MG','Tmag','t1','SLOGL'])
DR24 = ['D24' in x for x in cm['REF']]
LBV = ['lbv' in x for x in cm['SpC']]
BRC = ['B[e]SG' in x for x in cm['SpC']]

#plt.plot(cm['BR'],cm['MG'],'o'); plt.xlabel(r'G$_{B}$-G$_{R}$'); plt.ylabel(r'M$_{G}$')  
#plt.plot(cm['BR'][DR24],cm['MG'][DR24],'rx'); 
#plt.plot(cm['BR'][LBV],cm['MG'][LBV],'ko',mew=2,ms=14,mfc='None'); 
#plt.plot(cm['BR'][BRC],cm['MG'][BRC],'ks',mew=2,ms=14,mfc='None'); 

plt.plot(cm['JK'],cm['MK'],'o'); plt.xlabel(r'J-K'); plt.ylabel(r'M$_{K}$')  
plt.plot(cm['JK'][DR24],cm['MK'][DR24],'rx'); 
plt.plot(cm['JK'][LBV],cm['MK'][LBV],'ko',ms=14,mfc='None'); 
plt.plot(cm['JK'][BRC],cm['MK'][BRC],'ks',ms=14,mfc='None'); 

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
 

plt.gca().invert_yaxis() 
plt.show()

'''




