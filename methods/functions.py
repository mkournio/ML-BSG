from constants import *
import numpy as np
import lightkurve as lk
from scipy import stats
import re
import EntropyHub as EH
import matplotlib.pyplot as plt
from lightkurve.correctors import DesignMatrix,RegressionCorrector
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

TESS_pix_size = 21
GAIA_UPMARK = 64

class Gaia(object):

    def __init__(self,tpf,cat='Gaia3'):
        
        self.tpf = tpf
        self.ra = tpf.ra
        self.dec = tpf.dec
        if cat=='Gaia3' :
            self.cat = 'I/355/gaiadr3'
        elif cat=='Gaia2' :
            self.cat = 'I/345/gaia2'
        else :
            raise ValueError("Choose between Gaia3 and Gaia2.")
            
        self.get_prop()
        
        return

    def get_prop(self):

        DS3 = self._query()

        RA_pix,DE_pix = self.tpf.wcs.all_world2pix(DS3.RA_ICRS,DS3.DE_ICRS,0.01)
        RA_pix += self.tpf.column ; DE_pix += self.tpf.row
        RA_pix += 0.5 ; DE_pix += 0.5
        
        cond_box = (RA_pix>self.tpf.column) & (RA_pix<self.tpf.column+self.tpf.shape[2]) & \
                   (DE_pix>self.tpf.row) & (DE_pix<self.tpf.row+self.tpf.shape[1])

        DS3 = DS3[cond_box]
        self.RA_pix = RA_pix[cond_box]
        self.DE_pix = DE_pix[cond_box]
        self.BPRP  = DS3['BP-RP'].values
        RP    = DS3['RPmag'].values
        GaIDs  = DS3['Source'].values
        
        self.gaia_ind = None
        try:
            GA_query = Vizier(catalog='I/355/gaiadr3').query_region(SkyCoord(ra=self.ra,dec=self.dec,unit=(u.deg,u.deg),
                                                                             frame='icrs'), radius = 10 * u.arcsec)[0]
            GA_star =  GA_query['Source'][0]
            self.gaia_ind = np.where(GaIDs == GA_star)
            self.star_row = int(DE_pix[self.gaia_ind] - self.tpf.row)
            self.star_col = int(RA_pix[self.gaia_ind] - self.tpf.column)
            self.RPdiff =  RP - np.ma.min(GA_query['RPmag'])
        except:
            print('No Gaia counterparts retrieved')
            self.star_row = None
            self.star_col = None
            self.RPdiff =  RP - min(RP)
            
        return

    def _query(self):

        DS3 = Vizier(catalog=self.cat,columns=['*','+_r']); DS3.ROW_LIMIT = -1
        DS3_query = DS3.query_region(SkyCoord(ra=self.tpf.ra,dec=self.tpf.dec,unit=(u.deg,u.deg), frame='icrs'),
                                     radius = max(self.tpf.shape[1:]) * TESS_pix_size * u.arcsec)
        DS3 = DS3_query[0].to_pandas()

        return DS3

    def update_mask(self, check_mask, ref_p, min_thres = 3, update = False):

       min_diff = 100.
       
       if check_mask is not None :
           nmask = check_mask.copy()
           for i in range(nmask.shape[0]):
               for j in range(nmask.shape[1]):
                   if nmask[i,j]:
                       tpf_row = ref_p[0] + i
                       tpf_col = ref_p[1] + j
                       cond_box = (self.RA_pix>tpf_col) & (self.RA_pix<tpf_col+1) & (self.DE_pix>tpf_row) & (self.DE_pix<tpf_row+1)
                       
                       loc_diff = self.RPdiff[cond_box & (self.RPdiff != 0)]
                       loc_diff = loc_diff[~np.isnan(loc_diff)]
                       if len(loc_diff) > 0 and min(loc_diff) < min_diff :
                           min_diff = min(loc_diff)
                           if update and min_diff < min_thres :
                               for di in [max(0,i-1),i,min(i+1,nmask.shape[0]-1)]:
                                   for dj in [max(0,j-1),j,min(j+1,nmask.shape[1]-1)]:
                                       nmask[di,dj] = False
       else:
           nmask = None
           
       return nmask, min_diff
   
    def plot(self,ax):
        
        gaia_sizes = GAIA_UPMARK / 2**self.RPdiff
        ax.scatter(self.RA_pix,self.DE_pix,s=gaia_sizes, marker='.', c='c')
        
        if self.gaia_ind != None:
            ax.scatter(self.RA_pix[self.gaia_ind], self.DE_pix[self.gaia_ind], 
                       s=GAIA_UPMARK, marker='x', color='c', linewidths=2)
            
        return ax

   
def getmask(tpf, star_row = None, star_col = None, thres = 0.1):
    
    if star_row == None:
        star_row = int(0.5 * tpf.shape[1])
    if star_col == None:
        star_col = int(0.5 * tpf.shape[2])
        
    mask = np.zeros(tpf[0].shape[1:], dtype='bool')
    mask[star_row][star_col] = True
    
    flux_matr = np.nanmedian(tpf.flux, axis=0)
    cen_flux = flux_matr[star_row][star_col]
    
    rad = 1
    mask_size = len(mask)
    
    for c in np.arange(star_col,mask_size,1):
        
        col_break = -1
        
        for r in np.arange(star_row,-1,-1):
            
            max_flx = -1            
            if flux_matr[r][c] > thres * cen_flux:                
                mask[r][c] = True
                max_flx = 1
                col_break = 1
                
            if max_flx < 0.:
                break
            
        for r in np.arange(star_row,mask_size,1):
            
            max_flx = -1            
            if flux_matr[r][c] > thres * cen_flux:
                mask[r][c] = True
                max_flx = 1
                col_break = 1
                
            if max_flx < 0.: 
                break
            
        if col_break < 0.: 
            break
        
    for c in np.arange(star_col,-1,-1) :
        
        col_break = -1
        
        for r in np.arange(star_row,-1,-1):
            
            max_flx = -1
            if flux_matr[r][c] > thres * cen_flux:
                mask[r][c] = True
                max_flx = 1
                col_break = 1
                
            if max_flx < 0.: 
                break
            
        for r in np.arange(star_row,mask_size,1):
            
            max_flx = -1
            if flux_matr[r][c] > thres * cen_flux:
                mask[r][c] = True
                max_flx = 1
                col_break = 1
                
            if max_flx < 0.: 
                break
            
        if col_break < 0.:
            break
        
    return mask	

def dmatr(matr, pca_num):
    
    return DesignMatrix(matr,name='regressors').pca(pca_num).append_constant()

def lccor(tpf, mask, bkg_mask, pca_num, **kwargs):    
   
    lc = tpf.to_lightcurve(aperture_mask=mask)
    flux_mask = (lc.flux_err > 0) & (~np.isin(lc.quality,[1,4,16,32,1024,2048,16384]))
    lc = lc[flux_mask]
    
    rgr = tpf.flux[flux_mask][:, bkg_mask]
    dm = dmatr(rgr,pca_num)
    
    lcc = RegressionCorrector(lc).correct(dm)
    
    if 'gaps' in kwargs:
        gaps = kwargs['gaps']
        lcc.flux = set_nans(lcc.time.value, lcc.flux.value, gaps)

    return lcc
    
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
	
def sptype_to_temp(sp_v):
    
    t_v = []
    
    for sp in sp_v:
        
        try:
            t = float(re.findall(r"(?:\d*\.*\d+)", sp)[0])
        except:
            t = np.nan
            
        t_v.append(t)
        
    return t_v

def mask_outliers(m_arr, m = 3.):

    x = np.ma.copy(m_arr)	

    d = np.abs(x - np.ma.median(x))
    mdev = np.ma.median(d)
    s = d / (float(mdev) if mdev else 1.)

    x.mask[s>m] = True

    return x
	
def set_nans(time,flux,gaps):    
  
    if len(gaps) != 4:
        
        return flux
    
    gb,gl,gu,ge = gaps
    time_dif = time[1:] - time[:-1]
    ind_nan = np.argmax(time_dif)
    flux[ind_nan-gl:ind_nan+gu] = np.nan
    
    if gb > 0: flux[:gb] = np.nan
    if ge > 0: flux[-ge:] = np.nan
    
    return flux


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

# TODO - make the following as static methods 

def get_std(flux):
    
    return np.std(flux)

def get_mad(flux):
    
    return np.median(np.absolute(flux - np.median(flux)))

def get_skew(flux):
    
    return stats.skew(flux)

def get_eta(flux):

	succd = np.sum(np.diff(flux)**2) / (len(flux) - 1) 

	return succd / np.var(flux)

def get_iqr(flux):
    
    q75, q25 = np.percentile(flux, [75 ,25])
    
    return q75 - q25

def get_kurt(flux):
    
    return stats.kurtosis(flux)

def get_mse(flux, m = 2, tau = 15, tol = 0.2):
    
    if len(flux) > 40000:
        
        return  np.full(4, np.nan)
        
    Mobj = EH.MSobject('SampEn', m = m, r = tol * np.nanstd(flux), Logx = np.exp(1), Vcp = False)
    mse = EH.MSEn(flux, Mobj, Scales = tau, Methodx = 'coarse', RadNew = 0, Plotx = False)
    
    msx = np.arange(1,tau+1,1)
    msy = mse[0]
    
    idx = np.isfinite(msx) & np.isfinite(msy)
    msx = msx[idx]
    msy = msy[idx]
    
    std = get_std(msy)
    z2 = np.polyfit(msx,msy,deg=2)
    z1 = np.polyfit(msx,msy,deg=1)
    
    #f2 = np.poly1d(z2)
    #f1 = np.poly1d(z1)
    #plt.plot(msx,f2(msx))
    #plt.plot(msx,f1(msx))
    
    return  np.mean(msy**2), std, 2*z2[0], z1[0]

def k_cross(flux, kappa = 5):
    
    i = 0
    d_k = []
    nflux = np.copy(flux)
    while i < kappa:
        
        nflux = nflux - np.mean(nflux)
        
        clipped = np.copy(nflux)
        clipped[clipped >= 0] = 1
        clipped[clipped < 0] = 0
        
        d_k.append( np.sum(np.diff(clipped)**2)/len(clipped) )
        
        nflux = np.diff(nflux)
        
        i += 1
        
    return np.array(d_k)

def get_psi_sq(flux, kappa = 5):
    
    std = np.std(flux)
    d_star = k_cross(flux, kappa)
    
    delta_star = np.append(d_star[0], np.diff(d_star))
    
    wn = np.random.normal(0, std, size=len(flux))
    d_gauss = k_cross(wn,kappa)
    
    delta_gauss = np.append(d_gauss[0], np.diff(d_gauss))
    
    return np.sum( ( (delta_star - delta_gauss)**2 ) / delta_gauss)	

def save_two_col(x,y,filename):


	with open(filename,'w') as f:
		for i,j in zip(x,y):
    			f.write("%s %s\n" % (i,j))
	f.close()

	return

def save_three_col(lc, filename, meta, units = 'mag'):
    
    x = lc.time.value
    if units == 'mag':
        y = lc.dmag.value
        z = lc.dmag_err.value
        header = f'##BTJD-REF\tDMAG\tDMAG_ERR\n'
        
    with open(filename,'w') as f:
        
        for m in meta:
            f.write(f'#{m} = {meta[m]}\n')
            
        f.write(header)   
        for i,j,k in zip(x,y,z):
            f.write("%s %s %s\n" % (i,j,k))
            
    f.close()
    
    return


