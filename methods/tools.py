from constants import *
from astropy.io import fits
import datetime
import numpy as np
import lightkurve as lk
import warnings
import os

def check_header_key(hdu,key,val):
    if key in hdu.header:
        return hdu.header[key] == val
    else:
        print('Header does not contain {} keyword'.format(str(key)))
        return False
    
def get_hdu_from_keys(hdulist,**kwargs):
    
    if kwargs == {}:
        return None
    
    hdu=[]
    for i in hdulist:
        
        f = 0
        hdr = i.header
        for key, value in kwargs.items():
            if key in hdr and hdr[key] == value:
                f += 1
        if f == len(kwargs):
            hdu.append(i)
            
    if len(hdu) == 0:
        
        return None
    else:
        if len(hdu) > 1:
            warnings.warn('More than 1 headers found; returned the first.')
            
        return hdu[0]  
      
def get_sectors_from_hdulist(hdulist,**kwargs):
    
    sectors = []
    for h in hdulist:
        try:
            sectors.append(int(h.header['SECTOR']))
        except:
            pass
        
    return sorted(set(sectors))

def group_consecutive_hdus(hdulist,sectors):
    
    import more_itertools as mit    

    if len(sectors) > 1:
        return [list(group) for group in mit.consecutive_groups(hdulist, ordering = lambda x: sectors[hdulist.index(x)])]
    else:
        return [hdulist]
    
def get_ax_scaling(hdulist):
    
    h0 = len(hdulist[0])
    ax_scaling = []
    for h in hdulist[1:]:
        ax_scaling.append(100 * len(h) / h0)
        
    return ax_scaling        
    
def get_minmax_flux(hdulist, flux_key):
    
    mins=[]; maxs=[]       
    for hdu in hdulist:
        std = np.nanstd(hdu.data[flux_key])
        med = np.nanmedian(hdu.data[flux_key])
        mins.append(med-3*std)
        maxs.append(med+3*std)

    return min(mins), max(maxs) 

def get_fits_name(star,tic):
    
    return os.path.join(path_to_output_fits,'{}_{}.fits'.format(star,tic))


class FitsObject(object):

      from astropy.io import fits

      def __init__(self,fits_file = None,**kwargs):
           
       self.hdu_list = []
       self.fits_file = fits_file
       
       hdu = fits.PrimaryHDU()
       p_head = {
                "DATE": datetime.datetime.now().strftime("%Y-%m-%d"),
                "CREATOR": "KOURNIOTIS"       
                }       
       p_head.update(kwargs)
       for key in p_head:
         hdu.header["{}".format(key).upper()] = p_head[key]                    
       self.hdu_list.append(hdu)
       
       return
        
      def add_lc(self,
                 lc,
                 header_source = None,
                 **kwargs):
          
       if isinstance(lc, lk.LightCurve):
           time = lc.time.value
           flux = lc.flux.value
           try:
               flux_err = lc.flux_err.value
           except:
               flux_err = None
       elif isinstance(lc, np.ndarray):
           time = lc[0]
           flux = lc[1]
           if lc.shape[0] > 2:
               flux_err = lc[2]
           else:
               flux_err = None
       
       if isinstance(flux, np.ma.MaskedArray):
           flux = np.where(flux.mask,np.nan,flux)

       cols = []
       cols.append(fits.Column(name='time',format="D",unit=flux_units['time'],array=time))
       cols.append(fits.Column(name='flux',format="E",unit=flux_units['flux'],array=flux))
       if flux_err is not None :
        cols.append(fits.Column(name='flux_err',format="E",unit=flux_units['flux_err'],array=flux_err))
        
       #try:
       extra_cols = [c for c in lc.columns if c not in ['time','flux','flux_err']]
       for c in extra_cols:
               cols.append(fits.Column(name=c,format="E",unit=flux_units[c],array=lc[c]))
       #except:
        #   pass           
           
       coldefs = fits.ColDefs(cols)
       hdu = fits.BinTableHDU.from_columns(coldefs) 

       if header_source is not None:
           hdu.header = self._create_header_from_original(hdu.header,header_source)
       
       for key in kwargs:
         hdu.header["{}".format(key).upper()] = kwargs[key]
         
       self.hdu_list.append(hdu)      

       return
   
      def add_ls(self,
                 ls,
                 header_source = None,
                 **kwargs):
          
          
       return
   
      def add_field_from_tpf(self,tpf_file,timestamp = 0): 
       
       tpf = fits.open(tpf_file)          
       hdu = fits.ImageHDU(data=tpf[1].data['FLUX'][timestamp],header=tpf[1].header)
       self.hdu_list.append(hdu)      

       return
   
      def add_aperture_from_spoc(self,tess_lc):
       
       aperture = np.array(tess_lc.hdu[2].data)
       sector = tess_lc.hdu[0].header['SECTOR']
       # Convert flags to boolean
       #aperture = [[np.binary_repr(x, width = 8)[-2] == '1' for x in row] for row in aperture]
       hdu = fits.ImageHDU(data=aperture, header = tess_lc.hdu[2].header)
       hdu.header['SECTOR'] = str(sector)

       self.hdu_list.append(hdu)
       
       return          
   
      def _create_header_from_original(self, hdr,lc_fits):
       
       header_p = lc_fits[0].header
       header_b = lc_fits[1].header
           
       for key in tess_header_keys_p :
        try:
         hdr[key] = header_p[key]
        except:
         pass              
           
       for key in tess_header_keys_b :
        try:
         hdr[key] = header_b[key]
        except:
         pass    
            
       return hdr
      
      def close(self,overwrite=True):
      
       hdu_objs = fits.HDUList(self.hdu_list)  
       if self.fits_file is not None:
           hdu_objs.writeto(self.fits_file, overwrite=overwrite)
       hdu_objs.close()
       
       return 
     



'''

def colorbar(fig,vmin,vmax,label,cbar='jet',**kwargs):

		import matplotlib.cm as cm
		from matplotlib.colors import Normalize

		cmap = cm.get_cmap(cbar).reversed()
		cmap.set_under('white')
		norm = Normalize(vmin,vmax)

		cbar_ax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
		cmmapable = cm.ScalarMappable(norm, cmap)
		cmmapable.set_array([])
		cbar=fig.colorbar(cmmapable, cax=cbar_ax,**kwargs) 
		cbar.ax.set_title(label, size=CBAR_TITLE_SIZE)
		cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)

		return cmap, norm
'''
