from constants import *
import numpy as np

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
   
def simb_q(ids):

       from astroquery.simbad import Simbad
       
       s_obj = Simbad.query_objects(list(ids))
        
       return to_deg(s_obj['RA'],s_obj['DEC'])   
         
       
def sedscal(rad,dist):

	rad_arr = np.array(rad.filled(np.nan))		
	dist_arr = np.array(dist.filled(np.nan))

	return np.log10(rad_arr / dist_arr)

def float_list(a):
	return [float(x) for x in a]

def to_deg(ra,dec):

	from astropy.coordinates import SkyCoord
	from astropy import units as u
	
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

def radmass(logg,rad):
		
	logg_arr = np.array(logg.filled(np.nan))
	rad_arr = np.array(rad.filled(np.nan))

	return (10**logg_arr) * ((rad_arr * RSUN_TO_CM)**2) / (G_ACC * MSUN_TO_GR)

############# TESS LIGHTCURVES

def fit_pol(x,y,yerr,deg,unit = 'mag'):

	idx = np.isfinite(x) & np.isfinite(y)
	x = x[idx]; y = y[idx]; yerr = yerr[idx]

	p = np.poly1d(np.polyfit(x, y, deg))
	yfit = p(x)

	yfit_norm = y/yfit
	e_yfit_norm = yerr/yfit

	dm = -2.5 * np.log10(yfit_norm)
	e_dm = 2.5 * yerr / (LN10*y)

	if unit == 'mag':
		return x, yfit, dm, e_dm
	elif unit == 'norm':
		return x, yfit, ynorm, e_yfit_norm

def save_two_col(x,y,filename):


	with open(filename,'w') as f:
		for i,j in zip(x,y):
    			f.write("%s %s\n" % (i,j))
	f.close()

	return

def save_three_col(x,y,z,filename):

	with open(filename,'w') as f:
		for i,j,k in zip(x,y,z):
    			f.write("%s %s %s\n" % (i,j,k))
	f.close()

	return


