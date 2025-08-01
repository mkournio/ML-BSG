import lightkurve as lk
import numpy as np
import os
from constants import *
from methods.functions import *
from methods.plot import *
from methods.tools import FitsObject

class Extract(GridTemplate):
    
    # Class for 

    def __init__(self, data, plot_key = 'flux', **kwargs):

     self.data = data
     self.plot_key = plot_key
     self._validate()

     self.fSPOC = os.listdir(path_to_spoc_files)
     self.fFFI = os.listdir(path_to_tesscut_files)
     super().__init__(rows_page = PLOT_XLC_NROW, cols_page = PLOT_XLC_NCOL,
                      fig_xlabel= PLOT_XLABEL['lc'], fig_ylabel=PLOT_YLABEL[plot_key], **kwargs)	 
     #self._extract_LCs(**kwargs)

    def _validate(self):
     pass

    def lightcurves(self, mode = 'SPOC', **kwargs):        
        
     if kwargs.get('save_lcs') :         
         if not os.path.exists(path_to_lcs): 
            os.makedirs(path_to_lcs)
         
     if mode == 'SPOC':
         self._extract_spoc(**kwargs)     
     else:
         pass
     
     self.close_plot()
     
     return
     
    def _extract_spoc(self, time_bin_size = None,**kwargs):
        
        stars = self.data['STAR']
        tics = self.data['TIC']
        spcs = self.data['SpC']
        ras = self.data['RA']
        decs = self.data['DEC']

        for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
            
         files = [f for f in self.fSPOC if str(tic) in f]
         sects = [int(f.split('-')[1][1:]) for f in files]
         
         if len(files) > 0 : 
             
             s_sects,s_files = zip(*sorted(zip(sects,files)))
             
             filename = os.path.join(path_to_fits,'{}_{}.fits'.format(star,tic))
             if kwargs.get('save_fits'):
              ff = FitsObject(filename)

             for f, sect in zip(s_files,s_sects):
            
              print('{}: extracting Sector {} of TIC {} (SPOC)'.format(star,sect,tic))
   
              fits = os.listdir(path_to_spoc_files + f)[0]
              lc = lk.TessLightCurveFile(os.path.join(path_to_spoc_files + f,fits)).remove_nans().remove_outliers()
             
              ax_lc = self.GridAx()
              
              # Normalize raw by the polynomial fitting the binned light curve              
              bcoeff = None
              if time_bin_size is not None:
                lcb = lc.bin(time_bin_size = time_bin_size) 
                n_lcb, bcoeff = normalize(lcb)
                
              n_lc, _ = normalize(lc, coeff = bcoeff)
                 
              if kwargs.get('save_fits'): 
                    ff.add_lc(n_lc, header_source = lc, binning = 'F')
                    if time_bin_size is not None:
                        ff.add_lc(n_lcb, header_source = lc, binning = 'T', binsize = str(time_bin_size))
     
              # Plot the lightcurves
              plot_lc_single(ax_lc, n_lc, flux_key = self.plot_key, lc_type = 'spoc', m = '.')              
              if time_bin_size is not None:
                  plot_lc_single(ax_lc, n_lcb, flux_key = self.plot_key, lc_type = 'spoc_binned')
                  
              if self.plot_key == 'flux':
                  plot_lc_single(n_lc, ax=ax_lc, flux_key = 'fitmodel', m = '--', lc_type = 'fit')
              
              add_plot_features(ax_lc, mode = self.plot_key, upper_left=star,
                                lower_left=spc,lower_right='{} ({})'.format(tic,sect))        
     
 
             if kwargs.get('save_fits'):
                 
               if kwargs.get('extract_field'):
                   
                # Adds the field from the latest sector LC file   
                ff.add_aperture_from_spoc(lc)
                
                ffi_files = [j for j in self.fFFI if (str(ra)[:5] in j) and (str(dec)[:5] in j)]
                tpf_v = [f for f in ffi_files if int(f.split('-')[1][1:]) == sect]
                if len(tpf_v) > 0:
                    tpf_file = tpf_v[0]
                else:
                    tpf_file = ffi_files[0]
           
                ff.add_field_from_tpf(path_to_tesscut_files + tpf_file)
                
               ff.close(overwrite=True)         

        return 
                  
        
'''        
    def test(self):
		stars_id = self.data['STAR']
		stars_ra = self.data['RA']
		stars_dec = self.data['DEC']
		stars_tic = self.data['TIC']

		for s in range(len(stars_id)) :

			star = stars_id[s]
			ra = stars_ra[s]
			dec = stars_dec[s]
			tic = stars_tic[s]
			#print('Extraction for %s - TIC %s' % (star,tic))

			try:
				thmask = self.data['MASK'][s]
			except:
				thmask = FFI_MASK			

			spoc_files = [f for f in self.fSPOC if str(tic) in f]
			spoc_sects = [int(f.split('-')[1][1:]) for f in spoc_files]

			ffi_files = [j for j in self.fFFI if (str(ra)[:5] in j) and (str(dec)[:5] in j)]
			ffi_sects = [int(f.split('-')[1][1:]) for f in ffi_files]

			if not self.all_type:
			 try:
                          sspoc = self.data['SSPOC'][s]
                          sffi = self.data['SFFI'][s]
                          sspoc = sect_unmasked(sspoc); sffi = sect_unmasked(sffi);
                          spoc_files, spoc_sects = [x for x, y in zip(spoc_files,spoc_sects) if y in sspoc], sspoc
                          ffi_files, ffi_sects = [x for x, y in zip(ffi_files,ffi_sects) if y in sffi], sffi
			 except:
                          pass

			star_sects =  sorted(set(spoc_sects + ffi_sects))
			if len(star_sects) > 0: self._extract_single(star,ra,dec,star_sects,spoc_files,ffi_files,thmask)

		return

 	def _extract_single(self, star, ra, dec, star_sects, spoc_files, ffi_files, thmask):

		sect_old = -2 
		for sect in star_sects:

			ax_tpf = self.GridAx()
			ax_lc = self.GridAx()
			ax_lc_n = self.GridAx(); ax_lc_n.invert_yaxis()

			spoc_sect = [f for f in spoc_files if int(f.split('-')[1][1:]) == sect]
			ffi_sect = [f for f in ffi_files if int(f.split('-')[1][1:]) == sect]
			pca=1; ndeg=2

			spoc_mask = None
			spoc_ccdv = None

			if bool(spoc_sect):
				print('%s, sector %s: SPOC lightcurve processing..' % (star,sect))

				spoc_lc_file = lk.TessLightCurveFile(TESS_SPOC + spoc_sect[0])
				spoc_mask = spoc_lc_file.pipeline_mask
				spoc_ccdv = [spoc_lc_file.hdu[2].header['CRVAL2P'], spoc_lc_file.hdu[2].header['CRVAL1P']]

				ax_tpf.set_ylabel('%s (%s)' % (star,sect), fontsize=16); ax_tpf.set_xlabel('')

				spoc_lc = spoc_lc_file.get_lightcurve('PDCSAP_FLUX')

				x_act = spoc_lc.time
				y_act = spoc_lc.flux
				yerr_act = spoc_lc.flux_err

				ltype = 'SPOC' 		
				x_n, yfit,dm,e_dm= fit_pol(x_act,y_act,yerr_act,ndeg,unit='mag')
				ax_lc.plot(x_n,yfit,LC_COLOR[ltype],ls='--')
				ax_lc.plot(x_act,y_act,LC_COLOR[ltype])
				ax_lc_n.plot(x_n,dm,LC_COLOR[ltype])

				if self.save_files:
				 save_three_col(x_n,dm,e_dm,TESS_LC_PATH+'/%s_%s_DMAG_%s' % (star,sect,ltype))

				spoc_lc_file.hdu.close()

			if bool(ffi_sect):
			 try:			  
				print('%s, sector %s: FFI, extracting and processing..' % (star,sect))
				tpf = lk.TessTargetPixelFile(TESS_TPF_PATH + ffi_sect[0])
				bk_mask = ~tpf.create_threshold_mask(THR_BCKG, reference_pixel=None)
				bk_lc = tpf.to_lightcurve(aperture_mask=bk_mask).remove_nans()

				G = Gaia(tpf,ra,dec,cat='Gaia3')
				mask = getmask(tpf,G.star_row,G.star_col,thres=thmask); 
				if sect == 11 and star == 'HD125545'  : 
					print 'Manual mask for HD125545'
					mask = getmask(tpf,G.star_row,G.star_col,thres=0.27);

				#save_mask(mask,self.EXTLCpath+'/%s_%s_MASK' % (star,tpf.sector))
				#_, min_diff_FFI = G.update_mask(min_thres = 5, check_mask = mask, ref_p = [tpf.row, tpf.column], update = False)
				#_, min_diff_SPOC = G.update_mask(min_thres = 5, check_mask = spoc_mask, ref_p = spoc_ccdv, update = False)
				#print 'Min RPdiff for %s from sector %s: (FFI, SPOC) (%.2f, %.2f)' % (star,sect,min_diff_FFI, min_diff_SPOC)				
				tpf.plot(ax=ax_tpf,aperture_mask=mask,mask_color='#FD110D',show_colorbar=False,\
					bkg_mask=bk_mask,bkg_mask_color='w',title='',spoc_mask=spoc_mask,spoc_ccdv=spoc_ccdv)

				gaia_sizes = GAIA_UPMARK / 2**G.RPdiff
				ax_tpf.scatter(G.RA_pix,G.DE_pix,s=gaia_sizes, marker='.', c='c') 
				if G.gaia_ind != None: 
					ax_tpf.scatter(G.RA_pix[G.gaia_ind], G.DE_pix[G.gaia_ind], s=GAIA_UPMARK, marker='x', color='c', linewidths=2)

				ax_tpf.text(0.05,0.90,'%s' % thmask,color='c',size=12,transform=ax_tpf.transAxes)
				ax_tpf.set_ylabel('%s (%s)' % (star,sect), fontsize=16); ax_tpf.set_xlabel('')
				ax_tpf.set_xlim(tpf.column,tpf.column+tpf.shape[2])
				ax_tpf.set_ylim(tpf.row,tpf.row+tpf.shape[1])

				raw_lc = tpf.to_lightcurve(aperture_mask=mask).remove_nans()
				corrected_lc1 = lccor(tpf, mask, bk_mask, pca)
				x_act = corrected_lc1.time
				y_act = corrected_lc1.flux
				yerr_act = corrected_lc1.flux_err

				ltype = 'FFI'
				x_n, yfit,dm,e_dm= fit_pol(x_act,y_act,yerr_act,ndeg,unit='mag')
				ax_lc.plot(x_act,y_act,LC_COLOR[ltype])
				ax_lc.plot(x_n,yfit,LC_COLOR[ltype],ls='--')
				ax_lc_n.plot(x_n,dm,LC_COLOR[ltype])

				if self.save_files:
				 save_three_col(x_n,dm,e_dm,TESS_LC_PATH+'/%s_%s_DMAG_%s' % (star,sect,ltype))

				tpf.hdu.close()
		         except:
				print('FFI READING WAS SKIPPED - %s !' % ffi_sect[0])
				pass 
		return

'''



