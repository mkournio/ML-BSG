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

def round_to(x,ref):
	import math

	rounddown = int(math.floor(x / ref)) * ref
	roundup = int(math.ceil(x / ref)) * ref	

	return rounddown, roundup

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

def remove_slow_lcs(files):
    
    for f in files:
        if "a_fast" in f:
            try:
                files.remove(f.replace("a_fast","s"))
            except:
                pass
            
    return files

def normalize(lc, deg = 2, flux_key ="pdcsap_flux", coeff = None):
    
    import lightkurve as lk
    
    if isinstance(lc, lk.LightCurve):
        time = lc.time.value
        flux = lc[flux_key].value
        try:
            flux_err = lc[flux_key+"_err"].value
        except:
            flux_err = np.zeros(len(flux))
            
    elif isinstance(lc, np.ndarray):
        time = lc[0]
        flux = lc[1]
        if lc.shape[0] > 2:
            flux_err = lc[2]
        else:
            flux_err = np.zeros(len(flux))
    
    mask_nan = np.isfinite(flux)
    flux = np.ma.array(flux, mask=~mask_nan)
    
    if coeff is None:
        coeff = np.ma.polyfit(time, flux, deg)
     
    p = np.poly1d(coeff)
    flux_fit = p(time)
    
    nflux = flux/flux_fit
    e_nflux = flux_err/flux_fit
    
    nflux = np.where(nflux.mask,np.nan,nflux)
    dm = -2.5 * np.log10(nflux)
    e_dm = 2.5 * flux_err / (LN10*flux)
    
    new_lc = lk.LightCurve(time=time, flux=flux, flux_err = flux_err)
    new_lc.add_columns([nflux,e_nflux,dm,e_dm,flux_fit],
                       names=['nflux','nflux_err','dmag','dmag_err','fitmodel'])

    return new_lc, coeff

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


