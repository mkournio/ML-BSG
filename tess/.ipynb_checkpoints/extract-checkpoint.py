#from plot_methods import GridTemplate
#from functions import *
import lightkurve as lk
import numpy as np
import os 


class GridTemplate(object):
        pass

class TessLCs(GridTemplate):

	def __init__(self, data, plots_per_grid = 5, all_type = True, save_files = False, **kwargs):

		self.data = data
		self.all_type = all_type
		self.save_files = save_files				
		self._validate()

		#fl_id = 0 
		#while os.path.exists(TESS_LC_PATH+'EXTLC_%s' % fl_id): fl_id += 1
		#self.EXTLCpath = TESS_LC_PATH+'EXTLC_%s' % fl_id
		#os.mkdir(self.EXTLCpath)

		#self.fSPOC = os.listdir(TESS_SPOC)
		#self.fFFI = os.listdir(TESS_TPF_PATH)
		#super(ExtractLC,self).__init__(plots_per_grid, cols_page = PLOT_XLC_NCOL, fig_xlabel=PLOT_XLABEL['lc'], fig_ylabel='', **kwargs)	
		#self._extract_LCs(**kwargs)
		#self.GridClose()

	def _validate(self):
		pass

	def _extract_LCs(self, **kwargs) :

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





