from astropy.io import fits
import lightkurve as lk
from methods.functions import *
import matplotlib.pyplot as plt
from constants.paths import *
from tess.measures import *
from tess.extract import *
from tess.contamination import *
from tables.match import *
import os

rawf = os.listdir(path_to_spoc_files)
rawf = remove_slow_lcs(rawf)
rawfs = os.listdir(path_to_tess_spoc_files)
procf = os.listdir(path_to_output_fits)

def fits_raw(tic,sector,flux_column = 'sap_flux'):
    
    sector = 's%04d' % sector    
    for f in rawf:        
        if (str(tic) in f) and (sector in f):
            
            if 'fast' in f:
                lc_file = f'{path_to_spoc_files}{f}/{f}-lc.fits'
            else:
                lc_file = f'{path_to_spoc_files}{f}/{f}_lc.fits'  
                
            lc = lk.TessLightCurveFile(lc_file, flux_column = flux_column).remove_outliers()
            lc = lc[lc.quality==0]
           
            return lc
        
    for f in rawfs:
        if (str(tic) in f) and (sector in f):
            
            lc_file = f'{path_to_tess_spoc_files}{f}/{f[:-3]}_lc.fits'  
        
            lc = lk.TessLightCurveFile(lc_file, flux_column = flux_column).remove_outliers()
            lc = lc[lc.quality==0]
            
            return lc
       
        
    return None

def fits_proc(tic, header = 'all'):
    
    for f in procf:        
        if str(tic) in f:
            
            fl = fits.open(f'{path_to_output_fits}{f}')
            
            if isinstance(header,int):                
                return fl[header]
            else:
                return fl
            
            return 
            
    return None

def eval_polyfit(tic, sector, lc_bin = True, gaps = [], *args,**kwargs):
    
    lc = fits_raw(tic,sector)
    if lc_bin:
        lc = lc.bin(time_bin_size = 0.00694)
        
    for n in gaps:
        mask = (lc.time.value > n[0]) & (lc.time.value < n[1])
        lc.flux[mask] = np.nan
    lc = lc.remove_nans()  

    lcn = normalize_break(lc, *args,**kwargs)
    
    

    s = np.nanstd(lcn.nflux)

    fig,ax=plt.subplots(nrows=2,ncols=1) 
    ax[0].plot(lcn.time.value,lcn.flux) 
    ax[0].plot(lcn.time.value,lcn.trend)
    ax[1].plot(lcn.time.value,lcn.nflux)
    
    ax[1].set_ylim(1.-8*s,1.+8*s)
    
    return fig

def eval_flatten(tic,sector, f_win, lc_bin = True, break_tolerance = 5, 
                 gaps = [], **kwargs):
    
    lc = fits_raw(tic,sector)    
    if lc_bin:
        lc = lc.bin(time_bin_size = 0.00694)   
        
        
    for n in gaps:
        mask = (lc.time.value > n[0]) & (lc.time.value < n[1])
        lc.flux[mask] = np.nan
    lc = lc.remove_nans()        
    
    bin_size = np.nanmedian(lc.time[1:] - lc.time[0:-1]).value
    window_length = int(f_win / bin_size)
    break_tolerance = int(break_tolerance / bin_size)
    
    lcf, trend = lc.flatten(window_length = window_length, return_trend = True,
                            break_tolerance = break_tolerance, **kwargs)
    
    fig,ax=plt.subplots(nrows=2,ncols=1) 
    ax[0].plot(lc.time.value,lc.flux)
    ax[0].plot(trend.time.value,trend.flux)   
    ax[1].plot(lcf.time.value,lcf.flux)
    
    return fig




#xc = PreProcess(args).compile()
#cm = PostProcess(xc).xmatch().append('combined').write_to_csv('input_sample_ems')

    