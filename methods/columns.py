from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.table import Table,  hstack, vstack, Column, MaskedColumn
import numpy as np
from astropy.io import ascii
from .functions import *
from constants import *

def dict_func(func):
         def wrapper(*args,**kwargs):
           
           if len(args) != len(func(*args,**kwargs)):
              raise ValueError("Defined columns do not match the output of the called function")
              
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
             
      func_ra = lambda x: simb_q(x[id_k])[0]
      func_dec = lambda x: simb_q(x[id_k])[1]

      return [func_ra, func_dec]
      
@dict_func      
def diffcol(*nc,keys):

      func_v = []
      for kp in keys:
        m1_k, m2_k = kp
        f = lambda x, m1_k = m1_k, m2_k = m2_k : x[m1_k] - x[m2_k] 
        func_v.append(f)
    
      return func_v

@dict_func      
def meancol(*nc,keys):

     func_v = []
     for kp in keys:      
        f = lambda x, kp = kp: np.array(x[kp].to_pandas().mean(axis=1))
        func_v.append(f) 
            
     return func_v

@dict_func      
def coord_h2d(*nc, keys):

     ra_k, dec_k = keys
     
     func_ra = lambda x: to_deg(x[ra_k],x[dec_k])[0]
     func_dec = lambda x: to_deg(x[ra_k],x[dec_k])[1]
     
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
      np.seterr(divide = 'ignore') 
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
 
@dict_func     
def spc2t(*nc,keys):
    
    spc_key = keys[0]

    return [lambda x: sptype_to_temp(x[spc_key])]

def merged_col(col1,col2):

      mcol = Column(col2, dtype='object') 
      for i in range(len(col1)):
       ind = np.where(col1 == col1[i])
       mcol[i] = '|'.join(np.array(col2[ind]))

      return mcol  
            



