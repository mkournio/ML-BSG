from constants import *
from astropy.io import fits
import datetime
import numpy as np

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
                 time,
                 flux,
                 flux_err=None,
                 flux_column_name='flux',
                 tess_lc_file = None,
                 **kwargs):
       
       if isinstance(flux, np.ma.MaskedArray):
           flux = np.where(flux.mask,np.nan,flux)
      
       cols = []
       flux_column_name = flux_column_name.upper()
       cols.append(fits.Column(name='TIME',format="D",unit='d',array=time))
       cols.append(fits.Column(name=flux_column_name,format="E",unit=flux_unit[flux_column_name],array=flux))
       if flux_err is not None :
        cols.append(fits.Column(name=flux_column_name+'_ERR',format="E",unit=flux_unit[flux_column_name],array=flux_err))
        
       coldefs = fits.ColDefs(cols)
       hdu = fits.BinTableHDU.from_columns(coldefs) 

       if tess_lc_file is not None:
           hdu.header = self._create_header_from_original(hdu.header,tess_lc_file)
       
       for key in kwargs:
         hdu.header["{}".format(key).upper()] = kwargs[key]
         
       self.hdu_list.append(hdu)      

       return
   
      def add_field_from_tpf(self,tpf): 
          
       hdu = fits.ImageHDU(data=tpf[1].data['FLUX'][0],header=tpf[1].header)
       hdu.header['TYPE'] = "FIELD"
       self.hdu_list.append(hdu)      

       return
   
      def add_aperture_from_spoc(self,tess_lc_file):
       
       aperture = np.array(tess_lc_file.hdu[2].data)
       
       # Convert flags to boolean
       aperture = [[np.binary_repr(x, width = 8)[-2] == '1' for x in row] for row in aperture]
       hdu = fits.ImageHDU(data=aperture,header=tess_lc_file.hdu[2].header)

       self.hdu_list.append(hdu)
       
       return          
   
      def _create_header_from_original(self, hdr,tess_lc_file):
          
       header_p = tess_lc_file.hdu[0].header
       header_b = tess_lc_file.hdu[1].header
           
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
