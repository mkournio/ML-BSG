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
from astropy.table import Table

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
time_metrics = ['MAD_RAW','MAD','ETA','EMSE1','EMSE0']
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


############# ML METHODS
feats = Features()
feats.get_from_sectors(
    input_cat = cm,
    time_keys = time_metrics + ['SECTOR','CROWDSAP'], 
    freq_keys = freq_metrics,
    rn_keys = rn_metrics,
    calc_keys = ['JH','HK','KW4','W14','W24','W34','Q_JHK'],
    log_convert = ['IQR','ETA','W0','R0','PSI','MAD','MAD_RAW','TOP','MSE0','EMSE0'],
    save_output = None)
feats.merge_cand()

### CORNER PLOTS
#feats.corner_plot(plot_cols = time_metrics, hue = 'SpC',outlier_sigma=5.)
#feats.corner_plot(plot_cols = freq_metrics+rn_metrics, hue = 'SpC',outlier_sigma=5.)

var_cols = time_metrics + rn_metrics + freq_metrics
agg_cols = var_cols + ['Tmag','CROWDSAP','Q_JHK','HK','JH','Jmag','Hmag','Kmag']
ml_kwargs = {
    'agg_type' :'median', 'agg_cols': agg_cols, 'split_cand' : False,
    'var_cols' : var_cols, 'scaler_type' : 'standard', 'pca_components': 6,
    'kn' : 3, 'n_perm': 0, 'umap_min_dist': 0.1, 'min_cluster_size': 3
    }
   
TT = TexTab()
print(feats.df[['STAR','SECTOR'] + freq_metrics])

feats.aggregate(group_by = ['STAR','SpC','TIC','RA','DEC'], **ml_kwargs)

# SAMPLE PRINTING
#TT.TabSample(feats.df[['STAR','RA','DEC','SpC','Tmag','Q_JHK','TIC','CROWDSAP']].sort_values('RA'))

#### PAIR PLOTS
'''feats.pair_plot(['W0','Tmag'], star_labels = True)
fig, ax = plt.subplots(2, 1, sharex =  True, figsize = (8,16)); ax= ax.flatten()
feats.pair_plot(ax = ax[0], pair_c = ['HK','JH'])
feats.pair_plot(ax = ax[1], pair_c = ['HK','Q_JHK'], star_labels = True)
ax[0].set_xticklabels([]); ax[0].set_xlabel('')
fig.subplots_adjust(hspace=0.02); fig.tight_layout()'''


ml_kwargs['var_cols'].remove('W0')
#### CLASSIFIERS - REGRESSORS
#feats.knn_classify(**ml_kwargs)
#feats.knn_regress(regress_col='HK', exclude_labels=[],**ml_kwargs)
#nn_table = feats.nearest_neighbors(**ml_kwargs)#; print(nn_table)
#feats.umap_plot(umap_n = 20, **ml_kwargs)

#feats.hier_clustering(k = 2, **ml_kwargs)
#feats.hdbscan_clustering(**ml_kwargs)
#feats.umap_plot(cbar_col = ['DIST_PCA'], umap_n = 20, **ml_kwargs)
#feats.local_outlier_factor(**ml_kwargs)


#### UMAP PLOTS
fig, ax = plt.subplots(1, 2, figsize = (18,10)); ax= ax.flatten()
feats.umap_plot(ax=ax[0], umap_n = 6, **ml_kwargs)
feats.umap_plot(ax=ax[1], umap_n = 20, **ml_kwargs) 
fig.tight_layout(); fig.savefig('UMAP.eps', format = 'eps')

'''fig, ax = plt.subplots(1,3, figsize = (18,6)); ax= ax.flatten()
feats.umap_plot(ax=ax[0], cbar_col = ['Tmag'], umap_n = 6, **ml_kwargs)
feats.umap_plot(ax=ax[1], cbar_col = ['Q_JHK'], umap_n = 6, **ml_kwargs)
feats.umap_plot(ax=ax[2], cbar_col = ['CROWDSAP'], umap_n = 6, **ml_kwargs)
ax[1].set_yticklabels([]); ax[1].set_ylabel('')
ax[2].set_yticklabels([]); ax[2].set_ylabel('')
fig.subplots_adjust(wspace=0.0); fig.tight_layout()'''

