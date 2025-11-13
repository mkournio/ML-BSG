import lightkurve as lk
import numpy as np
import os
from constants import *
from methods.functions import *
from methods.plot import *
from methods.tools import FitsObject, get_fits_name
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def max_freq_sn(pg, sn_window, freq_mask):
    
    v_freq = pg.frequency.value
    v_pow = pg.power.value
    maxf = pg.frequency_at_max_power.value
    
    mask_low = (maxf - sn_window/2. < v_freq) & (v_freq < maxf - freq_mask)
    mask_up = (v_freq > maxf + freq_mask) & (maxf + sn_window/2. > v_freq) 
    sn_mask = mask_low | mask_up
    
    return max(pg.power) / np.std(v_pow[sn_mask])

def sinusoidal(freq,off,ampl,phase):
    
    if len(ampl) != len(phase):
        raise IndexError('Inconsistent lengths for amplitudes and phases.')
        
    freq = freq * np.array(range(1,len(ampl)+1))
   
    func = lambda x, x0: off + sum(ampl * np.sin(2 * np.pi * freq * (x-x0) + phase))  
    vfunc = np.vectorize(func)
    
    return lambda t: vfunc(t,t[0])


def lorentz(X,w,zero,tau,gamma):
    
    return w + ( zero / ( 1 + (2 * np.pi * tau * X)**gamma))

class Extract(GridTemplate):
    
    def __init__(self, 
                 data, 
                 plot_key = 'flux', 
                 **kwargs):

     self.data = data.copy()
     self.plot_key = plot_key
     self._validate()

     self.fSPOC = os.listdir(path_to_spoc_files)
     self.fTSPOC = os.listdir(path_to_tess_spoc_files)
     self.fFFI = os.listdir(path_to_tesscut_files)
     
     # probably the following line should move  to lightcurves functions
     super().__init__(rows_page = PLOT_XLC_NROW, cols_page = PLOT_XLC_NCOL,
                      fig_xlabel= PLOT_XLABEL['lc'], fig_ylabel=PLOT_YLABEL[plot_key], **kwargs)
     
     self.close_plot()

     return
 
    def _validate(self):
     pass

    def lightcurves(self, 
                    time_bin = [], 
                    mode = 'ALL', 
                    **kwargs):        
        
     if kwargs.get('save_lcs') :         
         if not os.path.exists(path_to_lcs): 
            os.makedirs(path_to_lcs)
            
     stars = self.data['STAR']
     tics = self.data['TIC']
     spcs = self.data['SpC']
     ras = self.data['RA']
     decs = self.data['DEC']
     
     for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
         
         spoc_files = [f for f in self.fSPOC if str(tic) in f]
         spoc_files = remove_slow_lcs(spoc_files)
         sec_spoc = [int(f.split('-')[1][1:]) for f in spoc_files]
         
         tspoc_files = [f for f in self.fTSPOC if str(tic) in f]
         sec_tspoc = [int(f.split('-')[2][1:5]) for f in tspoc_files]
         
         sect_all = sorted(set(np.concatenate((sec_spoc,sec_tspoc),axis=0)))
       
      #   print(star,tic,spoc_files,sec_spoc)
       #  print(star,tic,tspoc_files,sec_tspoc)
         
         if ((len(spoc_files) > 0) or (len(tspoc_files) > 0)) and kwargs.get('save_fits'):
                 self.ff = FitsObject(get_fits_name(star,tic))                 
       
         for sect in sect_all:             
             
             sect = 's%04d' % sect
             spoc_f = [f for f in spoc_files if sect in f]
             tspoc_f = [f for f in tspoc_files if sect in f]
             
             if len(spoc_f) > 0:                 
                 ax_lc = self._extract_lc(spoc_f[0], time_bin, **kwargs)
                 print('{}: extracted Sector {} of TIC {} (SPOC)'.format(star,sect,tic))
             elif (len(tspoc_f) > 0) and (mode == 'ALL'):
                 ax_lc = self._extract_lc(tspoc_f[0], time_bin, **kwargs)
                 print('{}: extracted Sector {} of TIC {} (TESS-SPOC)'.format(star,sect,tic))
                 
             add_plot_features(ax_lc, mode = self.plot_key, upper_left=star,
                                      lower_left=spc,lower_right='{} ({})'.format(tic,sect))
                 
         if hasattr(self,'ff'):
             
             # TO DO
             if kwargs.get('extract_field'):
                 
                 # Adds the field from the latest sector LC file
                 self.ff.add_aperture_from_spoc(lc)
                 
                 ffi_files = [j for j in self.fFFI if (str(ra)[:5] in j) and (str(dec)[:5] in j)]
                 tpf_v = [f for f in ffi_files if int(f.split('-')[1][1:]) == sect]
                 
                 if len(tpf_v) > 0:
                     tpf_file = tpf_v[0]
                 else:
                     tpf_file = ffi_files[0]
                     
                 self.ff.add_field_from_tpf(path_to_tesscut_files + tpf_file)
                 
             self.ff.close(overwrite=True) 
             del self.ff

     self.close_plot()
     
     return
     
    def _extract_lc(self, 
                    f, 
                    time_bin, 
                    **kwargs):
        
        if f.startswith('tess'):
            pipeline = 'spoc'
            lc_file = os.listdir(path_to_spoc_files + f)[0]
            path_to_input_file = os.path.join(path_to_spoc_files + f,lc_file)
        elif f.startswith('hlsp'):
            pipeline = 'tess-spoc'
            lc_file = os.listdir(path_to_tess_spoc_files + f)[0]
            path_to_input_file = os.path.join(path_to_tess_spoc_files + f,lc_file)
            
        lc = lk.TessLightCurveFile(path_to_input_file).remove_nans().remove_outliers()
        
        # Opening fits using astropy because LightKurve method for
        # retrieving header is depraced
        lc_fits = fits.open(path_to_input_file)
        
        ax_lc = self.GridAx()
        
        # Normalize raw by the polynomial fitting the binned light curve
        bcoeff = None
        binned_lcs = []
        
        if len(time_bin) > 0:
            time_bin = sorted(time_bin, reverse = True)            
            lcb = lc.bin(time_bin_size = time_bin[0])
            n_lcb, bcoeff = normalize(lcb)
            binned_lcs.append(n_lcb)
            
            for t in time_bin[1:]:
                lcb = lc.bin(time_bin_size = t)
                n_lcb, _ = normalize(lcb, coeff = bcoeff)
                binned_lcs.append(n_lcb)           
            
        n_lc, _ = normalize(lc, coeff = bcoeff)
        if kwargs.get('save_fits'):
            self.ff.add_lc(n_lc, header_source = lc_fits, binning = 'F', pipeline = pipeline)
            for t, n_lcb in zip(time_bin,binned_lcs) :
                self.ff.add_lc(n_lcb, header_source = lc_fits, binning = 'T', binsize = str(t), pipeline = pipeline)
                
        # Plot the lightcurves
        plot_lc_single(ax_lc, n_lc, flux_key = self.plot_key, lc_type = pipeline, m = '.')
        if len(time_bin) > 0:
            plot_lc_single(ax_lc, binned_lcs[0], flux_key = self.plot_key, lc_type = 'binned')
            
        if self.plot_key == 'flux':
            plot_lc_single(n_lc, ax=ax_lc, flux_key = 'fitmodel', m = '--', lc_type = 'fit')            
         
        lc_fits.close()
        
        return ax_lc
    
    def header_key(self, 
                   key = 'TIMEDEL', 
                   mode = 'ALL',  
                   **kwargs):
               
        stars = self.data['STAR']
        tics = self.data['TIC']
        spcs = self.data['SpC']
        ras = self.data['RA']
        decs = self.data['DEC']
        key_vect = []
        
        for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
            
            spoc_files = [f for f in self.fSPOC if str(tic) in f]
            spoc_files = remove_slow_lcs(spoc_files)
            sec_spoc = [int(f.split('-')[1][1:]) for f in spoc_files]
            
            tspoc_files = [f for f in self.fTSPOC if str(tic) in f]
            sec_tspoc = [int(f.split('-')[2][1:5]) for f in tspoc_files]
            
            sect_all = sorted(set(np.concatenate((sec_spoc,sec_tspoc),axis=0)))
          
         
            for sect in sect_all:             
                
                sect = 's%04d' % sect
                spoc_f = [f for f in spoc_files if sect in f]
                tspoc_f = [f for f in tspoc_files if sect in f]
                
                if len(spoc_f) > 0:
                    lc_file = os.listdir(path_to_spoc_files + spoc_f[0])[0]
                    path_to_input_file = os.path.join(path_to_spoc_files + spoc_f[0],lc_file)
                elif (len(tspoc_f) > 0) and (mode == 'ALL'):
                    lc_file = os.listdir(path_to_tess_spoc_files + tspoc_f[0])[0]
                    path_to_input_file = os.path.join(path_to_tess_spoc_files + tspoc_f[0],lc_file)
                    
                lc_fits = fits.open(path_to_input_file)
                key_vect.append(lc_fits[1].header[key])                
                lc_fits.close()
                
        
        
        print('Min {:.4f} - Max {:.4f} - Median {:.4f}'.format(min(key_vect),max(key_vect),np.nanmedian(key_vect)))
     #   print(sum(np.array(key_vect) >= 0.))
        plt.hist(key_vect); plt.show()
       
        return
    
    def periodograms(self,
                     stitched = False,
                     prewhitening = False,
                     **kwargs):
        
        stars = self.data['STAR']
        tics = self.data['TIC']
        
        for star, tics in zip(stars,tics):
            
            fits_files = [j for j in os.listdir(path_to_output_fits) if (star in j)]
            
            if len(fits_files) > 0:
                
                ff = fits.open(get_fits_name(star,tic))
               
                if stitched:
                    
                    # TO DO - calculate LS for stitched, 
                    # need to decide how it can be saved on the fits files                   
                    pass
                
                else:
                    for hdu in ff[1:]:
                        
                        self.lomb_scargle(hdu,**kwargs)                        
                        #add new hdu LS in fits file
                        
                    ff.writeto(get_fits_name(star,tic), overwrite=True)
                    ff.close()
                    
        return

    @staticmethod
    def lombscargle(f,
                  prewhitening = False,
                  flux_key = 'dmag', 
                  nterms = 1,
                  red_noise = True,
                  lc_plot = False,
                  ls_plot = False,
                  **kwargs):
        
        lc = f.copy()
        
        if isinstance(lc, lk.LightCurve):
            lc = lc.remove_nans()
        else: 
            if isinstance(lc, fits.BinTableHDU):
                time = lc.data['time']
                flux = lc.data[flux_key]
            elif (isinstance(lc, np.ndarray)) or (isinstance(lc, list)):
                time = lc[0]
                flux = lc[1]
            else:
               raise TypeError('Object light curve does not have supportive format!') 
                    
            flux = flux[20:-20]; time = time[20:-20]            
            lc = lk.LightCurve(time=time,flux=flux).remove_nans()
            
        ls_method = 'fast'
        if nterms > 1: ls_method = 'fastchi2'
        fit_model = []
            
        lc_ini = lc.copy()
        pg = lc_ini.to_periodogram(ls_method=ls_method,nterms=nterms)
        pg_ini = pg.copy()       

        model = pg.model(time=lc_ini.time,frequency=pg.frequency_at_max_power)        
        theta  = pg.model_params(frequency=pg.frequency_at_max_power)
        sn_window = kwargs.pop('sn_window',1.)
        freq_mask = kwargs.pop('freq_mask',1/27.)
        sn_val = max_freq_sn(pg, sn_window, freq_mask)
        
        fit_model.append([pg.frequency_at_max_power.value,sn_val.value,theta])
        #sinf=sinusoidal(m_freqs[-1], m_offsets[-1], m_amplitudes[-1], m_phases[-1])
        #siny = sinf(lc.time.value)
        
        if prewhitening:
            
            prw_sn = kwargs.pop('prw_sn',4.)
            prw_maxn = kwargs.pop('prw_maxn',30.)
            
            ind = 0
            while (ind < prw_maxn) and (sn_val >= prw_sn):
          
                prw_res = lc.copy()
                prw_res.flux = prw_res.flux - model.flux
                
                pg = prw_res.to_periodogram(ls_method=ls_method,nterms=nterms)                
                sn_val = max_freq_sn(pg, sn_window, freq_mask)
                if sn_val < prw_sn:
                    break
                
                model.flux = model.flux + pg.model(time=prw_res.time, frequency=pg.frequency_at_max_power).flux                
                theta  = pg.model_params(frequency=pg.frequency_at_max_power)  
                fit_model.append([pg.frequency_at_max_power.value,sn_val.value,theta])                
                #sinf=sinusoidal(m_freqs[-1], m_offsets[-1], m_amplitudes[-1], m_phases[-1])
                #siny += sinf(lc.time.value)              

                ind += 1
                
        if red_noise:
            popt, perr = Extract.rednoise(pg, fit_scale='log')
            print(popt,perr)
            if prewhitening: 
                popt_ini, perr_ini = Extract.rednoise(pg_ini, fit_scale='log') 
                print(popt_ini,perr)

            
        ax_lc = None
        if lc_plot:
            _, ax_lc = plt.subplots()
            
            ax_lc.plot(lc.time.value,lc.flux.value)
            ax_lc.plot(model.time.value,model.flux.value,'r')
            #ax_lc.plot(lc.time.value,siny,'g--')
                
        ax_ls = None
        if ls_plot:
            _, ax_ls = plt.subplots()
            x = pg.frequency.value
            y = np.log10(pg.power.value)
            ax_ls.plot(x,y,'k')
            
            if prewhitening:
                x = pg_ini.frequency.value
                y = np.log10(pg_ini.power.value)
                ax_ls.plot(x,y,'b')
                if red_noise: ax_ls.plot(x,np.log10(lorentz(x,*popt_ini)),'r--')
                
            if red_noise: ax_ls.plot(x,np.log10(lorentz(x, *popt)),'r-')
   
            ax_ls.set_xscale('log')
            ax_ls.set_xlim([0.05,x[-1]])                
            ax_ls.set_ylim([-6.5,None])                
            
        return fit_model, lc.time[0].value
    
    @staticmethod
    def rednoise(p,
                 fit_scale = 'log',
                 low_lim = 2/27., 
                 up_lim = np.inf,
                 **kwargs):
        
        pg = p.copy()
        
        nan_ind = np.isnan(pg.power)
        x = pg.frequency[~nan_ind].value 
        y = pg.power[~nan_ind].value
        fit_func = lorentz        
        if fit_scale == 'log':
            y = np.log10(y)
            fit_func = lambda x,w,z,t,g: np.log10(lorentz(x,w,z,t,g))
       
        mask = (np.array(x) > low_lim) & (np.array(x) < up_lim)   
        popt, pconv = curve_fit(fit_func,x[mask],y[mask], p0=[1e-5, 1e-3, 0.1, 3.])#,bounds=bounds)
        perr = np.sqrt(np.diag(pconv))
        
        return popt, perr
    
    
        
                  
        
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



