from constants import *
from astropy.io import fits
import datetime
import numpy as np
import lightkurve as lk
import warnings
import os
import astropy.units as u
from astropy.io import fits


def check_header_key(hdu,key,val):
    if key in hdu.header:
        return hdu.header[key] == val
    else:
        print('Header does not contain {} keyword'.format(str(key)))
        return False
    
def get_hdu_from_keys(hdulist,**kwargs):
    
    if kwargs == {}:
        return []
    
    hdu=[]
    for i in hdulist:
        
        f = 0
        hdr = i.header
        for key, value in kwargs.items():
            if (key in hdr) and (hdr[key] == value):
                    f +=1
        if f == len(kwargs):
            hdu.append(i)          
            
    return hdu
      
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
   
      @staticmethod   
      def append_lombscargle(hdu_list,
                             pg_tab,
                             rn_tab = [],
                             header_source = None,
                             **kwargs):
          
          if len(pg_tab) < 2:
              return

          cols = []
          cols.append(fits.Column(name='frequency',format="E",unit='1/d',array=pg_tab[0]))
          for i in range(1,len(pg_tab)):
              cols.append(fits.Column(name='p_{}'.format(i-1),format="E",unit='mag',array=pg_tab[i]))
              
          coldefs = fits.ColDefs(cols)
          hdu = fits.BinTableHDU.from_columns(coldefs)
          
          hdu.header['HDUTYPE'] = 'LOMBSCARGLE'
          if header_source is not None:
              hdu.header = FitsObject.create_header_from_selection(hdu.header, header_source, keys = pg_header_keys)
              
          for i in range(len(rn_tab)):              
              for e, rnk in enumerate(['WHITE','RED','TAU','GAMMA']):
                  v = rn_tab[i][e]; e_v = rn_tab[i][e+4]
                  if np.isnan(v): v = None; e_v=None
                  hdu.header["{}_{}".format(rnk,i)] = (v, e_v)
                  
          #print(hdu.header)
          hdu_list.append(hdu)
          
          return
      
      @staticmethod 
      def append_frequencies(hdu_list,
                             frequencies,
                             t0,
                             header_source = None,
                             **kwargs):
          
          nterms = int(len(frequencies[0][3:])/2)
          tr_freqs = list(map(list, zip(*frequencies)))

          cols = []
          cols.append(fits.Column(name='frequency',format="E", unit='1/d',array=tr_freqs[0]))
          cols.append(fits.Column(name='snr',format="E",array=tr_freqs[1]))
          cols.append(fits.Column(name='offset',format="E",unit='mag',array=tr_freqs[2]))
  
          t = 0
          while t < nterms:
              cols.append(fits.Column(name='amplitude_{}'.format(t),format="E",unit='mag',array=tr_freqs[t+3]))
              cols.append(fits.Column(name='phase_{}'.format(t),format="E",unit='rad',array=tr_freqs[t+3+nterms]))
              t += 1
          
          coldefs = fits.ColDefs(cols)
          hdu = fits.BinTableHDU.from_columns(coldefs)
          
          hdu.header['HDUTYPE'] = 'FREQUENCIES'
          hdu.header['TSTART'] = t0
          if header_source is not None:
              hdu.header = FitsObject.create_header_from_selection(hdu.header, header_source, keys = pg_header_keys)
 
          hdu_list.append(hdu)
          #print(hdu.data)          
          
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
   
      @staticmethod   
      def create_header_from_selection(hdr,source_hdr,keys):
          
       for key in keys:
           
           if key in source_hdr:
               hdr[key] = source_hdr[key]        
          
       return hdr
   
    
      
      def close(self,overwrite=True):
      
       hdu_objs = fits.HDUList(self.hdu_list)  
       if self.fits_file is not None:
           hdu_objs.writeto(self.fits_file, overwrite=overwrite)
       hdu_objs.close()
       
       return 
     

class FitsList(object):
    
    def __init__(self, tab):
        
        self.tab = tab.copy()
        
        
        return
    
    def add_header_keys(self, key_dict, primary= False, replace = False):
        
        if len(key_dict) < 1:
            
            return
        
        print('Adding keys from headers..')

        ltab = self.tab
        for star, tic in zip(ltab['STAR'],ltab['TIC']):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]            
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                if primary :
                    hdus = [ff[0]]
                else:
                    hdus = ff[1:]                    
                    
                for hdu in hdus:
                    hdr = hdu.header
                    for k in key_dict:                        
                        if k in hdr and not replace:
                            pass
                        else:
                            hdr[k] = key_dict[k]
                
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()
                
        return
    
    def remove_header_keys(self, keys = [], primary = False):
        
        if len(keys) < 1:
            
            return 
        
        print('Removing keys from headers..')
        
        ltab = self.tab       
        for star, tic in zip(ltab['STAR'],ltab['TIC']):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]            
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                
                if primary :
                    hdus = [ff[0]]
                else:
                    hdus = ff[1:]  
                
                for hdu in hdus:
                    hdr = hdu.header
                    for k in keys:
                        if k in hdr: 
                            del hdr[k]
                            
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()
                   
        return    
    
    def remove_hdu(self, hdutypes = []):
        
        if len(hdutypes) < 1:
            
            return
        
        ltab = self.tab
        for star, tic in zip(ltab['STAR'],ltab['TIC']):
            
            filename = [f for f in os.listdir(path_to_output_fits) if star in f]            
            if len(filename) > 0 :
                
                ff = fits.open(os.path.join(path_to_output_fits,filename[0]))
                new_list = [x for x in ff if (('HDUTYPE' not in x.header) or (x.header['HDUTYPE'] not in hdutypes))]
                ff = fits.HDUList(new_list)
                
                ff.writeto(os.path.join(path_to_output_fits,filename[0]), overwrite=True)
                ff.close()
                
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
