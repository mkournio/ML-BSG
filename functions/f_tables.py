from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.table import Table,  hstack, vstack, Column, MaskedColumn
from astroquery.simbad import Simbad
import numpy as np
from astropy.io import ascii

from constants import *

def dict_func(func):
         def wrapper(*args,**kwargs):
           
           dict_f = {x:y for x,y in zip(args,func(*args,**kwargs))}
           
           return dict_f
                   
         return wrapper
         
         
#Functions for tables and arrays 
@dict_func      
def slogl(*nc,keys):

      tk, gk = keys
      sL_sun = (5778**4.) / (10**4.44)     
      
      return [lambda x: np.log10( ((x[tk]**4.)/ (10**x[gk])) / sL_sun )]

@dict_func      
def sbcoord_d(*nc,keys):

      id_k = keys[0]
      def simb_q(ids):

        s_obj = Simbad.query_objects(list(ids))
        
        return to_deg(s_obj['RA'],s_obj['DEC']) 
             
      func_ra = lambda x: simb_q(x[id_k])[0]
      func_dec = lambda x: simb_q(x[id_k])[1]

      return [func_ra, func_dec]
 
@dict_func           
def dist(*nc,keys):

      g_k, d2_k = keys

      def gal_to_d(gal,dist2):
             
         d_arr = np.ma.array(np.zeros(len(gal)),mask=True)
         for i, j, k in zip(range(len(d_arr)),gal,dist2):
             try:
              d_arr[i] = int(GALAXIES[j]['D'])
             except:
              d_arr[i] = k        

         return d_arr
      
      return [lambda x: gal_to_d(x[g_k],x[d2_k])]

@dict_func            
def absmag(*nc,keys):

      d = keys[-1]
      keys = keys[:-1]
      
      func_v = []
      for k in keys:
        f = lambda x, k=k: x[k] - (5 * np.log10(x[d])) + 5
        func_v.append(f)
               
      return func_v    
 
@dict_func     
def galloc(*nc,keys):

      ra_k, dec_k = keys
 
      def coord_to_gal(ra,dec):
      
       gal=[]
       for r, d in zip(ra,dec):
        
        g_star='MW'
        for g, p in GALAXIES.items(): 
         if ((r-p['RA'])**2) + ((d-p['DEC'])**2) < (p['RAD']**2):  
          g_star = g
          break
        gal.append(g_star) 

       return gal
    
      return [lambda x: coord_to_gal(x[ra_k],x[dec_k])]
      
@dict_func            
def galcoord(*nc,keys):

      ra_k, dec_k = keys
      
      def gconv(ra,dec):
      
       gal = SkyCoord(ra=np.array(ra)*u.degree, dec=np.array(dec)*u.degree, frame='icrs').galactic
       
       return gal.l, gal.b
      
      func_l = lambda x: gconv(x[ra_k],x[dec_k])[0]
      func_b = lambda x: gconv(x[ra_k],x[dec_k])[1]      
      
      return [func_l, func_b]

@dict_func     
def pow10(*nc,keys):
      
      func_v = []
      for k in keys:
        f = lambda x, k=k : 10**x[k]
        func_v.append(f)
               
      return func_v
              
def merged_col(col1,col2):

      mcol = Column(col2, dtype='object') 
      for i in range(len(col1)):
       ind = np.where(col1 == col1[i])
       mcol[i] = '|'.join(np.array(col2[ind]))

      return mcol  
            
def radmass(logg,rad):
		
	logg_arr = np.array(logg.filled(np.nan))
	rad_arr = np.array(rad.filled(np.nan))

	return (10**logg_arr) * ((rad_arr * RSUN_TO_CM)**2) / (G_ACC * MSUN_TO_GR)

def sedscal(rad,dist):

	rad_arr = np.array(rad.filled(np.nan))		
	dist_arr = np.array(dist.filled(np.nan))

	return np.log10(rad_arr / dist_arr)

def float_list(a):
	return [float(x) for x in a]

def to_deg(ra,dec):
	
	c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

	return c.ra.degree, c.dec.degree

def mask_outliers(m_arr, m = 3.):

    x = np.ma.copy(m_arr)	

    d = np.abs(x - np.ma.median(x))
    mdev = np.ma.median(d)
    s = d / (float(mdev) if mdev else 1.)

    x.mask[s>m] = True

    return x
	

def fill_array(a, fillvalue = 0):

	import itertools

	return np.array(list(itertools.izip_longest(*a, fillvalue=fillvalue))).T

