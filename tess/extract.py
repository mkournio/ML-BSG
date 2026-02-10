import lightkurve as lk
import numpy as np
import os
from constants.io import *
from methods.functions import *
from methods.plot import *
from methods.tools import FitsObject, get_fits_name
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings; warnings.filterwarnings("ignore")
from methods.tools import *
from tables.io import ls_params
from lmfit import minimize, Parameters, report_fit
from matplotlib import patches

def max_freq_sn(pg, sn_window, freq_mask, mode='std'):
    
    v_freq = pg.frequency.value
    v_pow = pg.power.value
    maxf = pg.frequency_at_max_power.value
    
    mask_low = (maxf - sn_window/2. < v_freq) & (v_freq < maxf - freq_mask)
    mask_up = (v_freq > maxf + freq_mask) & (maxf + sn_window/2. > v_freq) 
    sn_mask = mask_low | mask_up
    
    if mode == 'std':
        
        return max(pg.power) / np.std(v_pow[sn_mask])
    
    elif mode == 'mean':
        
        return max(pg.power) / np.mean(v_pow[sn_mask])

def noise_stats(pg, frequency, sn_window, freq_mask):
    
    v_freq = pg.frequency.value
    v_pow = pg.power.value
    
    mask_low = (frequency - sn_window/2. < v_freq) & (v_freq < frequency - freq_mask)
    mask_up  = (frequency + freq_mask < v_freq) & (v_freq < frequency + sn_window/2.) 
    mask = mask_low | mask_up
    
    return np.mean(v_pow[mask]), np.std(v_pow[mask])

def percentile_mask(tpf,q):
    
    medf = np.nanmedian(tpf.flux,axis=0).value
    medf_val = medf.flatten()
    
    p = np.nanpercentile(medf_val, q)
    
    return medf < p

def plot_bkg_aperture(ax,tpf,bkg_mask):
    
    b_mask = tpf._parse_aperture_mask(bkg_mask)
    in_aperture = np.where(b_mask)
    ap_row = in_aperture[0] + tpf.row - 0.5
    ap_col = in_aperture[1] + tpf.column - 0.5 
    for ii in range(len(ap_row)):
        
        rect=patches.Rectangle((ap_col[ii],ap_row[ii]),1,1, fill=False, hatch="\\", color='w')
        ax.add_patch(rect)
        
    return ax
    

def fourier_series(t,params):
    
    p = params.valuesdict()
    
    if p['nmod'] == 0:
        
        return np.zeros(len(t))
    
    func = np.zeros(len(t))
    for m in range(p['nmod']):
        
        func += p[f'offset_{m}']
        for h in range(1,p['nterms']+1):

            func += p[f'ampl_{m}_{h}'] * np.sin(2 * np.pi * h * p[f'f_{m}'] * (t-t[0]) + p[f'phase_{m}_{h}'])
            
    return func

def fit_residuals(params, x, y, y_err = []): 
    
    if len(y_err) == 0:
        y_err = np.ones(len(y))

    return  (fourier_series(x,params) - y)/y_err


def lorentz(X,w,zero,tau,gamma):
    
    return w + ( zero / ( 1 + (2 * np.pi * tau * X)**gamma))

def pg_model_params(mod_pg, mod_freq=None):
    
    if mod_freq == None:
        mod_freq = mod_pg.frequency_at_max_power
        
    p = mod_pg._LS_object.model_parameters(frequency=mod_freq)
    if hasattr(p,'value'): 
        p = p.value

    theta = np.full(len(p),np.nan)
    theta[0] = p[0]
    for i in range(1,len(p),2):
        theta[i] = np.sqrt((p[i]**2)+(p[i+1]**2))
        theta[i+1] = np.arctan2(p[i+1],p[i])
        
    return theta   

def get_lc_from_filename(f,**kwargs):
    
    if isinstance(f, lk.LightCurve):
        lc = f.copy()
        lc = lc.remove_nans()        
    else:
        if isinstance(f, fits.BinTableHDU):
            lc = f.copy()
            flux_key = kwargs.get('flux_key','dmag')
            time = lc.data['time']
            flux = lc.data[flux_key]
            
            try:
                flux_err = lc.data[flux_key+'_err']                
            except:
                flux_err = np.zeros(len(time))
                
        else:
            if isinstance(f,str) and f.endswith('txt'):
                lc = np.genfromtxt(f,unpack=True)
            elif (isinstance(f, np.ndarray)) or (isinstance(lc, list)):
                lc = f.copy()                
            else:
                raise TypeError('Object light curve does not have supportive format!')
      
            time = lc[0]
            flux = lc[1]
            try:
                flux_err = lc[2]
            except:
                flux_err = np.zeros(len(time))
                
            
        flux = flux[20:-20]; time = time[20:-20]; flux_err = flux_err[20:-20]
        lc = lk.LightCurve(time=time,flux=flux,flux_err=flux_err).remove_nans()

    return lc    


class Extract(GridTemplate):
    
    def __init__(self, 
                 data = [], 
                 plot_key = 'flux', 
                 **kwargs):

     self.data = data.copy()
     self.plot_key = plot_key
     self._validate()     
     
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
            
     fSPOC = os.listdir(path_to_spoc_files)
     fTSPOC = os.listdir(path_to_tess_spoc_files)
     fFFI = os.listdir(path_to_tesscut_files)
                 
     stars = self.data['STAR']
     tics = self.data['TIC']
     spcs = self.data['SpC']
     ras = self.data['RA']
     decs = self.data['DEC']
     
     for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
         
         spoc_files = [f for f in fSPOC if str(tic) in f]
         spoc_files = remove_slow_lcs(spoc_files)
         sec_spoc = [int(f.split('-')[1][1:]) for f in spoc_files]
         
         tspoc_files = [f for f in fTSPOC if str(tic) in f]
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
                 path_to_input_file = os.path.join(path_to_spoc_files, spoc_f[0], spoc_f[0]+'_lc.fits')
                 ax_lc = self._extract_lc(path_to_input_file, time_bin, **kwargs)
                 print('{}: extracted Sector {} of TIC {} (SPOC)'.format(star,sect,tic))
             elif (len(tspoc_f) > 0) and (mode == 'ALL'):
                 path_to_input_file = os.path.join(path_to_tess_spoc_files, tspoc_f[0], tspoc_f[0][:-3]+'_lc.fits')
                 ax_lc = self._extract_lc(path_to_input_file, time_bin, **kwargs)
                 print('{}: extracted Sector {} of TIC {} (TESS-SPOC)'.format(star,sect,tic))
                 
             add_plot_features(ax_lc, mode = self.plot_key, upper_left=star,
                                      lower_left=spc,lower_right='{} ({})'.format(tic,sect))
                 
         if hasattr(self,'ff'):
             
             # TO DO
             if kwargs.get('extract_field'):
                 
                 # Adds the field from the latest sector LC file
                 self.ff.add_aperture_from_spoc(lc)
                 
                 ffi_files = [j for j in fFFI if (str(ra)[:5] in j) and (str(dec)[:5] in j)]
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
     
    def _extract_lc(
            self,
            path_to_input_file,
            time_bin,
            **kwargs):
        
        filename = path_to_input_file.split('/')[-1]
        if filename.startswith('tess'):
            pipeline = 'spoc'
        elif filename.startswith('hlsp'):
            pipeline = 'tess-spoc'
        else:
            pipeline = 'unknown'
             
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
        
        fSPOC = os.listdir(path_to_spoc_files)
        fTSPOC = os.listdir(path_to_tess_spoc_files)
               
        stars = self.data['STAR']
        tics = self.data['TIC']
        spcs = self.data['SpC']
        ras = self.data['RA']
        decs = self.data['DEC']
        key_vect = []
        
        for star, spc, tic, ra, dec in zip(stars,spcs,tics,ras,decs):
            
            spoc_files = [f for f in fSPOC if str(tic) in f]
            spoc_files = remove_slow_lcs(spoc_files)
            sec_spoc = [int(f.split('-')[1][1:]) for f in spoc_files]
            
            tspoc_files = [f for f in fTSPOC if str(tic) in f]
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
                     bin_size = '10m',
                     stitched = False,
                     nterms = 3,
                     prew = True,
                     maximum_frequency = None,                     
                     **kwargs):
        
        if bin_size not in ['raw','10m','30m']:
            raise Exception('Set bin_size among raw, 10m, 30m. Aborting..')
            
        if bin_size == '10m':
            bin_size = 0.00694
        elif bin_size == '30m':
            bin_size = 0.02083

        
        stars = self.data['STAR']
        tics = self.data['TIC']
        
        for star, tic in zip(stars,tics):
            
            fits_files = [j for j in os.listdir(path_to_output_fits) if ((star in j) and (str(tic) in j))]
            
            if len(fits_files) > 0:
                
                ff = fits.open(get_fits_name(star,tic))
                print('Extracting LS data for {} - {}'.format(star,tic))
              
                if stitched:                    
                    # TO DO - calculate LS for stitched, 
                    # need to decide how it can be saved on the fits files                   
                    pass
                
                else:
                    hdus = get_hdu_from_keys(ff[1:], TTYPE1 = 'time', BINSIZE = str(bin_size))
                    for hdu in hdus:
                        t0, freq_tab, pg_tab, rn_tab = self.lombscargle(hdu, prew=prew, nterms=nterms, maximum_frequency = maximum_frequency)
                        FitsObject.append_lombscargle(ff, pg_tab, rn_tab, header_source = hdu.header)
                        FitsObject.append_frequencies(ff, freq_tab, t0, header_source = hdu.header)
                        
                ff.writeto(get_fits_name(star,tic), overwrite=True)
                ff.close()
                
        return
    
    @staticmethod 
    def from_tesscut(
            f,   
            mask = None,
            thres_ape = 0.2,
            thres_bkg = 1e-4,
            pca = 2,
            lc_bin = None,
            lc_fit_deg = 2,
            gaia_overlay = False,  
            ax_tpf = None,
            ax_lc = None,
            show_extract = True,
            save_output = False,
            **kwargs): 
        
        tpf = lk.TessTargetPixelFile(f)
        ra = tpf.ra
        dec = tpf.dec
        sect = tpf.sector
        cdn = float(tpf.hdu[1].header['TIMEDEL'])
        
        if type(thres_bkg) == float:            
            bkg_mask = ~tpf.create_threshold_mask(thres_bkg, reference_pixel=None)
        elif thres_bkg.startswith('q'):
            bkg_mask = percentile_mask(tpf, float(thres_bkg[1:]))        
        
        bkg_vec = tpf.flux[:, bkg_mask]
        bkg_vec_n = [x/np.nanmedian(x) for x in bkg_vec.T]
        
        g_size = kwargs.pop('gaia_size',64)        
        G = Gaia(tpf,ra,dec,cat='Gaia3')
        mask = getmask(tpf,thres=thres_ape)
         
        fig, ax = plt.subplots(2,2,figsize=(16,10))
        fig.suptitle(f)
        ax = ax.flatten()        
            
        tpf.plot(ax=ax[0],frame=200,aperture_mask=mask,mask_color='#FD110D',show_colorbar=False,title='')
        plot_bkg_aperture(ax[0],tpf,bkg_mask)
        if gaia_overlay:
            G.plot(ax[0])
                            
        ax[0].text(0.05,0.90,'%s' % thres_ape,color='c',size=12,transform=ax[0].transAxes)
        ax[0].set_xlim(tpf.column,tpf.column+tpf.shape[2]-1)
        ax[0].set_ylim(tpf.row,tpf.row+tpf.shape[1]-1)
        
        ax[1].plot(tpf.time.value,bkg_vec[:,:30],'.')
        ax[1].set_ylabel(r'F$_{bkg}$ [e$^{-}$/s]')
        ax[1].set_xlabel(f"Time - {tpf.header['BJDREFI']} [BTJD d]")
        
        off = -0.03
        for p in [1,2,3]:
            lcc = lccor(tpf, mask, bkg_mask, pca_num = p, **kwargs)
            ax[2].plot(lcc.time.value,
                       off + (lcc.flux.value/np.nanmedian(lcc.flux.value)),'.',
                       label=f'PCA {p}')
            off += 0.03

        ax[2].set_ylabel(r'F/F$_{med}$ + const.')
        ax[2].set_xlabel(f"Time - {tpf.header['BJDREFI']} [BTJD d]")
        ax[2].legend(loc=4)    
            
        lcc = lccor(tpf, mask, bkg_mask, pca_num = pca, **kwargs)
        if isinstance(lc_bin,float):
            lcc = lcc.bin(time_bin_size = lc_bin)
        if isinstance(lc_fit_deg,int):
            lcc = normalize(lcc,flux_key ="flux")[0]
            
       # ax[2].plot(lcc.time.value,lcc.flux.value,'b.',label = f'PCA {pca}')
       # ax[2].plot(lcc.time.value,lcc.fitmodel.value,'r--')
       # ax[2].legend()
       # ax[2].set_ylabel(r'Flux (e$^{-}$/s)')
       # ax[2].set_xlabel(f"Time - {tpf.header['BJDREFI']} [BTJD d]")

        ax[3].plot(lcc.time.value,lcc.dmag.value,'k.')
        ax[3].set_ylabel(r'$\Delta$m [mag]')
        ax[3].set_xlabel(f"Time - {tpf.header['BJDREFI']} [BTJD d]")
        ax[3].invert_yaxis()
        
        
        filename = f"{ra:.5f}_{dec:.5f}_s{sect}_cdn{cdn:.3f}_pca{pca}_ffi{thres_ape}"
        
        if save_output:
            plt.savefig(filename + '_lc.png')
            save_three_col(lcc.time.value,lcc.dmag.value,lcc.dmag_err.value,
                           filename+'_lc.txt')   
            
        tpf.hdu.close()
        
        return fig

    @staticmethod
    def lombscargle(
            f,
            prew = False,
            flux_key = 'dmag',
            nterms = 3,
            maximum_frequency=None,
            red_noise = True,
            plot_lc = False,
            plot_ls = False,
            save_output = None,
            **kwargs):        
       
        lc = f.copy()
        meta = {}
        sector = ''
        
        if isinstance(lc, lk.LightCurve):
            lc = lc.remove_nans()
        else: 
            if isinstance(lc, fits.BinTableHDU):
                time = lc.data['time']
                flux = lc.data[flux_key]
               # try:
                hdr = lc.header
                sector = '_s{}'.format(hdr['SECTOR'])
                meta = {x:hdr[x] for x in meta_keys_ls if x in hdr}
                #except:
                 #   pass                
            elif (isinstance(lc, np.ndarray)) or (isinstance(lc, list)):
                time = lc[0]
                flux = lc[1]
            else:
               raise TypeError('Object light curve does not have supportive format!') 
                    
            flux = flux[20:-20]; time = time[20:-20]            
            lc = lk.LightCurve(time=time,flux=flux).remove_nans()            

        if isinstance(save_output,str):
            save_output += '{}_n{}'.format(sector,nterms)
            
        ls_method = 'fast'
        if nterms > 1: ls_method = 'fastchi2'
        sn_window = kwargs.pop('sn_window',1.)
        freq_mask = kwargs.pop('freq_mask',1/27.)
            
        lc_ini = lc.copy()
        pg = lc_ini.to_periodogram(ls_method=ls_method,nterms=nterms,maximum_frequency=maximum_frequency)
        
        pg_tab = [pg.frequency.value,pg.power.value]
        model = pg.model(time=lc_ini.time,frequency=pg.frequency_at_max_power)  
        theta = pg_model_params(pg)
        
        sn_val = max_freq_sn(pg, sn_window, freq_mask)   

        fit_model = [[pg.frequency_at_max_power.value,sn_val.value,*theta]]        
           #sinf=sinusoidal(m_freqs[-1], m_offsets[-1], m_amplitudes[-1], m_phases[-1])
           #siny = sinf(lc.time.value)
        
        if prew:            
            prw_sn = kwargs.pop('prw_sn',4.)
            prw_maxn = kwargs.pop('prw_maxn',30.)
            
            ind = 0
            while (ind < prw_maxn) and (sn_val >= prw_sn):                
          
                prw_res = lc.copy()
                prw_res.flux = prw_res.flux - model.flux
                
                #fig = plt.figure()
                #plt.plot(pg.frequency.value,np.log10(pg.power.value))
                #fig.show()                
                pg = prw_res.to_periodogram(ls_method=ls_method,nterms=nterms,maximum_frequency=maximum_frequency)                
                sn_val = max_freq_sn(pg, sn_window, freq_mask)
                if sn_val < prw_sn:
                    break
                
                model.flux = model.flux + pg.model(time=prw_res.time, frequency=pg.frequency_at_max_power).flux    
                theta = pg_model_params(pg)
                fit_model.append([pg.frequency_at_max_power.value,sn_val.value,*theta])   
                #sinf=sinusoidal(m_freqs[-1], m_offsets[-1], m_amplitudes[-1], m_phases[-1])
                #siny += sinf(lc.time.value)              

                ind += 1
                
            pg_tab.append(pg.power.value)
                
        rn_model = [] 
        if red_noise:
            for p in pg_tab[1:]:
                popt, perr = Extract.rednoise(x=pg_tab[0], y=p, fit_scale='log')
                rn_model.append(np.hstack([popt,perr]))
                
        if plot_lc:
            plt.rcParams.update(PLOT_PARAMS['lc'])
            
            if 'ax_lc' in kwargs:
                ax_lc = kwargs['ax_lc']
            else:
                _, ax_lc = plt.subplots(figsize=(10, 5))
            
            ax_lc.plot(lc.time.value,lc.flux.value,'b',label='observed')
            ax_lc.plot(model.time.value,model.flux.value,'r',label='model')
            ax_lc.set_xlabel(PLOT_XLABEL['lc'])
            ax_lc.set_ylabel(PLOT_YLABEL['dmag'])
            ax_lc.legend()
           
            if isinstance(save_output,str):
                plt.savefig(save_output + '_lc.png')
            #ax_lc.plot(lc.time.value,siny,'g--')
              
        if plot_ls:
            if 'ax_ls' in kwargs:
                ax_ls = kwargs['ax_ls']
            else:
                _, ax_ls = plt.subplots(figsize=(10, 7))            
              
            alpha = np.linspace(1,0.4,len(pg_tab)-1)  
            for i in range(1,len(pg_tab)):
                x=pg_tab[0]                
                ax_ls.plot(x,np.log10(pg_tab[i]),'k',alpha=alpha[i-1])
                
            for i in range(len(rn_model)):                
                rn_prop = rn_model[i][:4]
                ax_ls.plot(x,np.log10(lorentz(x,*rn_prop)),'r',alpha=alpha[i])
                    
            ax_ls.set_xscale('log')
            ax_ls.set_xlim([0.05,x[-1]])                
            ax_ls.set_ylim([-6.4,None])
            ax_ls.set_xlabel(PLOT_XLABEL['ls'])
            ax_ls.set_ylabel(PLOT_YLABEL['ls'])
            
            if isinstance(save_output,str):
                plt.savefig(save_output + '_ls.png')
            
        if isinstance(save_output,str):
            ls_params(freq_tab=fit_model,rn_tab=rn_model,
                      output=save_output + '_params', meta = meta, **kwargs)            
        
            
        return lc.time[0].value, fit_model, pg_tab, rn_model
 
    @staticmethod
    def lombscargle_global(
            f,
            nprew = 30,
            term_sn = 4.,
            term_criterion = 'window',
            opt_step = 3,
            opt_range = 0.1,
            nterms = 3,
            maximum_frequency=None,
            red_noise = True,
            plot_lc = False,
            plot_ls = False,
            save_output = None,
            **kwargs):
        
        lc = get_lc_from_filename(f,**kwargs)
        
        if nprew < 1:
            raise Exception('Define at least one step for frequency extraction (nprew). Aborting..')
        
        meta = {}
        sector = ''         
        if isinstance(save_output,str):
            save_output += '{}_n{}'.format(sector,nterms)
            
        ls_method = 'fast'
        if nterms > 1: ls_method = 'fastchi2'
        sn_window = kwargs.pop('sn_window',1.)
        freq_mask = kwargs.pop('freq_mask',1/27.)        
 
        fit_model = [] 
        pg_tab = []
        prw_res = lc.copy()
        model = lk.LightCurve(time = prw_res.time, flux = np.zeros(len(prw_res.flux)))

        step = 0
        params = Parameters()
        
        params.add('nterms',nterms,vary=False)
        params.add('nmod',step,vary=False)
        save_rn=[]

        rel_change = []
        while step < nprew:
            
            pg = prw_res.to_periodogram(ls_method=ls_method,
                                        nterms=nterms,
                                        maximum_frequency=maximum_frequency)
            theta = pg_model_params(pg)
            
            
           # w_ini =  np.nanmedian(pg.power.value[(10 < pg.frequency.value) & (pg.frequency.value < 20)])
           # z_ini =  np.nanmedian(pg.power.value[(0.07 < pg.frequency.value) & (pg.frequency.value < 0.12)])
           # print(w_ini,z_ini)

            popt_rn, perr_rn, bic = Extract.rednoise(
                    x=pg.frequency.value, 
                    y=pg.power.value,
                    fit_scale='log',
                    npar = step * (2*nterms + 1) + 1,
                    nt = len(prw_res.time)
                    )
            print(popt_rn)
            
            if step == 0: 
                pg_tab = [pg.frequency.value,pg.power.value]
            
            if step > 0:
                rel_change.append(abs(a0_old-popt_rn[0]-popt_rn[1])/a0_old)            
            a0_old =  popt_rn[0]  + popt_rn[1] 
            save_rn.append(a0_old)
            
            n_m, n_std = noise_stats(pg,pg.frequency_at_max_power.value,sn_window,freq_mask)
            if np.nanmean(rel_change[-3:]) < 0.1 :
                print(f'Step {step} converged')
                sn_val = max(pg.power)/lorentz(pg.frequency_at_max_power.value,*popt_rn)
            else:
                sn_val = max(pg.power)/n_std
                
            if sn_val < term_sn:
                break
            
            params['nmod'].set(value=step + 1)            
            params.add(f'f_{step}', value=pg.frequency_at_max_power.value,
                       min = pg.frequency_at_max_power.value - freq_mask,
                       max = pg.frequency_at_max_power.value + freq_mask)            
            params.add(f'offset_{step}', value=theta[0], 
                       min = theta[0]*(1.-opt_range),
                       max = theta[0]*(1.+opt_range),
                       )
            h=1
            for j in range(1,len(theta),2):
                    params.add(f'ampl_{step}_{h}', value = theta[j],
                               min = theta[j]*(1.-opt_range),
                               max = theta[j]*(1.+opt_range),
                               )
                    params.add(f'phase_{step}_{h}', value = theta[j+1],
                               min = theta[j+1]*(1.-opt_range),
                               max = theta[j+1]*(1.+opt_range),
                               )                    
                    h += 1                
            
            if (step+1) % opt_step == 0:
                
                result = minimize(fit_residuals,params,
                              args=(lc.time.value,lc.flux.value,lc.flux_err.value))
                params = result.params
                for p in params:
                    params[p].vary = False
            
            fit_model.append([pg.frequency_at_max_power.value,sn_val.value,*theta])
            model.flux = model.flux + pg.model(time=prw_res.time, frequency=pg.frequency_at_max_power).flux
            prw_res.flux = lc.flux - model.flux
            
            step += 1
            
        plt.plot(save_rn,'o')
        plt.show()
            
       # print(params,fit_model)   
        pg_tab.append(pg.power.value)          

    #    result = minimize(fit_residuals,params,
    #                     args=(lc.time.value,lc.flux.value,lc.flux_err.value))
    #    params = result.params
    #    print(params)
       # print(fit_model)

        sinf1 = fourier_series(lc.time.value,params)
                
        rn_model = [] 
        if red_noise:
            for p in pg_tab[1:]:
                popt, perr,_ = Extract.rednoise(x=pg_tab[0], y=p, fit_scale='log')
                rn_model.append(np.hstack([popt,perr]))
                
        if plot_lc:
            plt.rcParams.update(PLOT_PARAMS['lc'])
            
            if 'ax_lc' in kwargs:
                ax_lc = kwargs['ax_lc']
            else:
                _, ax_lc = plt.subplots(figsize=(10, 5))
            
            ax_lc.plot(lc.time.value,lc.flux.value,'k',label='observed')
            ax_lc.plot(lc.time.value,sinf1,'r',label=f'opt{opt_step} model')
         #   ax_lc.plot(lc.time.value,sinf2,'r',label='model2')
            ax_lc.plot(model.time.value,model.flux.value,'c',label='opt1 model')
           

            ax_lc.set_xlabel(PLOT_XLABEL['lc'])
            ax_lc.set_ylabel(PLOT_YLABEL['dmag'])
            
            ax_lc.invert_yaxis() 
            ax_lc.legend()
           
            if isinstance(save_output,str):
                plt.savefig(save_output + '_lc.png')
            #ax_lc.plot(lc.time.value,siny,'g--')
              
        if plot_ls:
            if 'ax_ls' in kwargs:
                ax_ls = kwargs['ax_ls']
            else:
                _, ax_ls = plt.subplots(figsize=(10, 7))            
              
            alpha = np.linspace(1,0.4,len(pg_tab)-1)  
            for i in range(1,len(pg_tab)):
                x=pg_tab[0]                
                ax_ls.plot(x,np.log10(pg_tab[i]),'k',alpha=alpha[i-1])
                
            for i in range(len(rn_model)):                
                rn_prop = rn_model[i][:4]
                ax_ls.plot(x,np.log10(lorentz(x,*rn_prop)),'r',alpha=alpha[i])
            #ax_ls.plot(x,np.log10(lorentz(x,3.55209592e-05,3.55209592e-05,6.43796300e-02,2.88226612e+00)),'b--')
               

            ax_ls.set_xscale('log')
            ax_ls.set_xlim([0.05,x[-1]])                
            ax_ls.set_ylim([-6.4,None])
            ax_ls.set_xlabel(PLOT_XLABEL['ls'])
            ax_ls.set_ylabel(PLOT_YLABEL['ls'])
            
            if isinstance(save_output,str):
                plt.savefig(save_output + '_ls.png')
            
        if isinstance(save_output,str):
            ls_params(freq_tab=fit_model,rn_tab=rn_model,
                      output=save_output + '_params', meta = meta, **kwargs)            
        
            
        return lc.time[0].value, fit_model, pg_tab, rn_model

      
    @staticmethod
    def rednoise(x,
                 y,
                 fit_scale = 'log',
                 low_lim = 2/27., 
                 up_lim = 30.,
                 npar = 0.,
                 nt = 2000,
                 **kwargs):       
        
        df = npar + 4
        
        nan_ind = np.isnan(y)
        x = x[~nan_ind]
        y = y[~nan_ind]
        
        fit_func = lorentz        
        if fit_scale == 'log':
            y = np.log10(y)
            fit_func = lambda x,w,z,t,g: np.log10(lorentz(x,w,z,t,g))

       
        mask = (np.array(x) > low_lim) & (np.array(x) < up_lim)
        try:
            popt, pconv = curve_fit(fit_func,x[mask],y[mask], p0=[3e-5, 1e-3, 0.1, 3], method= 'lm', maxfev=300)#,bounds=bounds)
            perr = np.sqrt(np.diag(pconv))
                
            mse = sum((y - np.log10(lorentz(x,*popt)))**2)/(len(y)*(0.434**2))  
            bic = len(y)*np.log(mse) + df * np.log(nt)           
            
        except RuntimeError:
            #popt, pconv = curve_fit(fit_func,x[mask],y[mask], p0=[3e-5,0.,np.inf,np.inf], method= 'lm', maxfev=300)#,bounds=bounds)
            #perr = np.sqrt(np.diag(pconv))

            popt = np.nan * np.empty(4)
            perr = np.nan * np.empty(4)            
            bic = np.nan
        
        return popt, perr, bic
    
       
                  
        
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



