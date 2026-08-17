from tables import *
from methods import *
import args 
import matplotlib.pyplot as plt
import numpy as np
from astroquery.simbad import Simbad
#from evolution import *
import os
from astropy.coordinates import SkyCoord
from constants import *
import pandas as pd
from tess import *
import numpy as np
from astropy.io import ascii

#xc = PreProcess(args).compile()
#GDW =  [(not any(y in x for y in ('III','V','O9','A0','A1','A2','A3','ON','OC'))) or ('LBV' in x) for x in xc['SpC']]
#xc = xc[GDW]
#cm = PostProcess(xc).xmatch().append('combined').write_to_csv('input_sample_ems')
cm = ascii.read('input_sample_ems')
cm.sort(['RA'])


# FILTERING THE SAMPLE
LOC = [('MW' in x) for x in cm['GAL']]
#LOC = [('MW' in x) or ('LMC' in x) or ('SMC' in x) for x in cm['GAL']]
cm = cm[LOC]

#r = np.where(cm['STAR']=='HD62623')[0][0]; cm = cm[r:r+1]
#cstars = ['HD96918','HR5171','HR8752','6 CAS','RHO CAS']
#cstars = ['HD62623','P Cyg']
#cstars = [
#    'MWC137','HD80077','GG Car','HR5171','WRAY 16-137','[GKF2010] MN48','[B61] 2',
#    'P Cyg','V439 Cyg','MWC349','6 CAS']
#r= [np.where(cm['STAR']==x)[0][0] for x in cstars]; cm = cm[r]

################ QUERYING FROM MAST
# LIGHTUCRVES
#mast_query(cm, download_dir='data/')
# TPFs
#mast_query(cm, product="Targetpixelfile")
#download_tpfs(cm, frame = 0, del_original = True)
# CBVS
#download_cbvs(cm)
# CBV VALIDATION
#CBVs(cm).validate(custom_types='lc_types',time_out = 200)


############# EXTRACTION 
# TIME DOMAIN

'''LCs = Extract(data=cm, 
              plot_key='flux',
              plot_name='x_ems', 
              rows_page = 4, cols_page = 1,
              output_format='png')
LCs.lightcurves(time_bin = [0.00694],
                type_file = 'lc_types', 
                gap_file = 'lc_gaps', 
                save_fits = True, extract_field = False)
#'''
# FREQUENCY DOMAIN
'''PGs = Extract(data=cm, 
              plot_name='xf_ems', 
              fig_xlabel = '', fig_ylabel = '', 
              rows_page = 4, cols_page = 2,
              figsize = (12,22),
              output_format='png')
PGs.periodograms(snr_file = 'ls_snr',
                 maximum_frequency=40)
#'''
############# VISUALIZATION
# LIGHTCURVES
'''
LC = Visualize(data=cm,
                plot_name='vn_ems', 
                plot_key='dmag', 
                figsize = (32,20),
                rows_page=6, 
                cols_page=6, 
                output_format='png').lightcurves(models=True, trend = True)
'''
# PERIODOGRAMS
'''
LS = Visualize(data=cm,
                plot_name='ls_ems', 
                plot_key='ls',
                figsize = (32,20),
                rows_page=6, 
                cols_page=6, 
                output_format='png').periodograms()
'''

############# METRICS
time_metrics = ['EMSE1','EMSE0','MAD','MAD_RAW','ETA']
freq_metrics = ['WFM','WFD']
rn_metrics = ['W0','R0','TAU','GAMMA']

# RESETING - REMOVING
#fl = FitsList(cm); fl.add_header_keys(key_dict={'HDUTYPE':'LIGHTCURVE'})
#fl = FitsList(cm); fl.remove_hdu(hdutypes=['FREQUENCIES','PERIODOGRAMS'])

# TIME DOMAIN
#td = TimeDomain(data = cm, measures = time_metrics).calculate()
# FREQUENCY DOMAIN
#fl = FitsList(cm); fl.remove_header_keys(keys = freq_metrics)
#fd = FrequencyDomain(data = cm, measures = freq_metrics).calculate(min_freq = 0.1)

  


feats = Features()
ftab = feats.get_from_sectors(
    input_cat = cm,
    time_keys = time_metrics + ['SECTOR','CROWDSAP'], 
    freq_keys = freq_metrics,
    rn_keys = rn_metrics,
    calc_keys = ['JH','HK','KW4','W14','W24','W34','Q_JHK'],
    log_convert = ['IQR','PSI','ETA','MAD','MAD_RAW','TOP','W0','R0','MSE0','EMSE0'],
    save_output = None)


#print(ftab[['STAR','SECTOR'] + freq_metrics])
 
#ftab_agg = feats._aggregate(cols = rn_metrics, mode='median',
#                     group_by = ['STAR','SpC','TIC'],
#                     save_output = None) 
#print(ftab_agg[['STAR'] + rn_metrics])

#print(ftab)

#feats.pair_plot(plot_cols = time_metrics, hue = 'SpC',aggregate_type='none',outlier_sigma=5)
#feats.pair_plot(plot_cols = freq_metrics, hue = 'SpC',aggregate_type='none')
#feats.pair_plot(plot_cols = rn_metrics, hue = 'SpC',aggregate_type='none',outlier_sigma=10.)
#feats.pair_plot(plot_cols = time_metrics+freq_metrics+rn_metrics, hue = 'SpC',aggregate_type='none',outlier_sigma=10.)

#feats.pair_plot_single(pair_cols=['Q_JHK','HK'])

#feats.pair_plot(plot_cols = time_metrics + frequency_metrics + rn_metrics,
#                hue = 'SpC',aggregate_type='none')


pca = 7
min_dist = 0.09
feats.knn_classify(var_cols = time_metrics+rn_metrics+freq_metrics,
                 aggregate_type='median', scaler_type = 'standard',
                 pca_components=pca,n_perm=0)

feats.knn_regress(var_cols = time_metrics+rn_metrics+freq_metrics,
                  regress_col='Tmag', aggregate_type='median', 
                  scaler_type = 'standard', pca_components=pca,n_perm=1000)


#print(feats)
'''

fig, ax = plt.subplots(1,2); ax= ax.flatten()
feats.umap_plot(ax=ax[0], var_cols = time_metrics+rn_metrics+freq_metrics,
                aggregate_type = 'median', scaler_type = 'standard',
                n_neighbors=6, min_dist=min_dist, pca_components=pca)
feats.umap_plot(ax=ax[1], var_cols = time_metrics+rn_metrics+freq_metrics,
                aggregate_type = 'median', scaler_type = 'standard',
                n_neighbors=20, min_dist=min_dist, pca_components=pca)

fig, ax = plt.subplots(1,3); ax= ax.flatten()
feats.umap_plot(ax=ax[0], var_cols = time_metrics+rn_metrics+freq_metrics,
                cbar_col = ['Tmag'], 
                aggregate_type = 'median', scaler_type = 'standard',
                n_neighbors=20, min_dist=min_dist, pca_components=pca)
feats.umap_plot(ax=ax[1], var_cols = time_metrics+rn_metrics+freq_metrics,
                cbar_col = ['CROWDSAP'], 
                aggregate_type = 'median', scaler_type = 'standard',
                n_neighbors=20, min_dist=min_dist, pca_components=pca)
feats.umap_plot(ax=ax[2], var_cols = time_metrics+rn_metrics+freq_metrics,
                cbar_col = ['Q_JHK'], 
                aggregate_type = 'median', scaler_type = 'standard',
                n_neighbors=20, min_dist=min_dist, pca_components=pca)
plt.tight_layout(wspace=0,hspace=0)
'''