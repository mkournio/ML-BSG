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

CLASS = [('B[e]SG' in x) or ('LBV' in x) for x in cm['SpC']]
cm = cm[CLASS]

CROBJ = ['Hen 3-519','[GKF2010] MN44','WRAY 16-232','V439 Cyg','MWC349']
crmask = [(x not in CROBJ) for x in cm['STAR']]
cm = cm[crmask]

#cm = cm[:1]

#r = np.where(cm['STAR']=='GG Car')[0][0] 
#cm = cm[r:r+1]
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
                save_fits = True, extract_field = False)'''
#'''
# FREQUENCY DOMAIN
'''PGs = Extract(data=cm, 
              plot_name='xf_ems', 
              fig_xlabel = '', fig_ylabel = '', 
              rows_page = 4, cols_page = 2,
              figsize = (12,22),
              output_format='png')
PGs.periodograms(snr_file = 'ls_snr',
                 maximum_frequency=40)'''
#
############# VISUALIZATION
# LIGHTCURVES

'''LC = Visualize(data=cm,
                plot_name='lc', 
                plot_key='flux', 
                figsize = (25,30),
                rows_page=8, 
                cols_page=3, 
                output_format='pdf').lightcurves(preview = 'preview_fig')'''
                                                 
'''LC = Visualize(data=cm,
                plot_name='nlc', 
                plot_key='dmag', 
                figsize = (25,30),
                rows_page=8, 
                cols_page=3, 
                output_format='pdf').lightcurves(models = True, preview = 'preview_fig')'''

# PERIODOGRAMS

'''LS = Visualize(data=cm,
                plot_name='ls', 
                plot_key='ls',
                figsize = (25,30),
                rows_page=8, 
                cols_page=3, 
                output_format='pdf').periodograms(preview = 'preview_fig')'''


############# METRICS
time_metrics = ['MAD_RAW','MAD','ETA','EMSE0','EMSE1']
freq_metrics = ['WFM','WFD']
rn_metrics = ['R0','TAU','GAMMA']
#tl_metrics = ['TLPW','TLLW','TLPS','TLLS']
# RESETING - REMOVING
#fl = FitsList(cm); fl.add_header_keys(key_dict={'HDUTYPE':'LIGHTCURVE'})
#fl = FitsList(cm); fl.remove_hdu(hdutypes=['FREQUENCIES','PERIODOGRAMS'])

# TIME DOMAIN
#fl = FitsList(cm); fl.remove_header_keys(keys = ['EMSE0','EMSE1','EMSE2'])
#td = TimeDomain(data = cm, measures = time_metrics).calculate()
# FREQUENCY DOMAIN
#fl = FitsList(cm); fl.remove_header_keys(keys = freq_metrics)
#fd = FrequencyDomain(data = cm, measures = freq_metrics).calculate(min_freq = 0.1)

#fl = FitsList(cm)
#fd = FrequencyDomain(data = cm, measures = ['FSTEPS']).calculate(min_freq = 0.1)


# CONTAMINATION - TESS LOCALIZE
#fl = FitsList(cm); fl.remove_header_keys(keys = tl_metrics)
#ct = Contamination(data = cm, measures = tl_metrics).calculate()

############# ML METHODS

feats = Features()
feats_tab= feats.get_from_sectors(
    input_cat = cm,
    time_keys = time_metrics + ['SECTOR','CROWDSAP'], 
    freq_keys = freq_metrics,
    rn_keys = rn_metrics,
    calc_keys = ['JH','HK','KW4','W14','W24','W34','Q_JHK'],
    log_convert = ['ETA','W0','R0','PSI','MAD','MAD_RAW'],
    save_output = None)

#feats_tab[['STAR','TIC','SECTOR','SpC','CROWDSAP','Tmag','Q_JHK'] + time_metrics + freq_metrics +rn_metrics].to_csv('feats_all0.csv',index=False)

#feats.corner_plot(plot_cols = time_metrics, hue = 'SpC', merge_cand=True, outlier_sigma=50.)
### CORNER PLOTS
#feats.corner_plot(plot_cols = time_metrics, hue = 'SpC', merge_cand=True, outlier_sigma=5.,save_output='cp_td.pdf')
#feats.corner_plot(plot_cols = freq_metrics+rn_metrics, hue = 'SpC', merge_cand=True, outlier_sigma=5.,save_output='cp_fd.pdf')

var_cols = time_metrics + freq_metrics + rn_metrics
agg_cols = var_cols + ['Tmag','CROWDSAP','Q_JHK','HK','JH','Jmag','Hmag','Kmag']
ml_kwargs = {
    'agg_type' :'mean', 'agg_cols': agg_cols, 'split_cand' : False,
    'var_cols' : var_cols, 'scaler_type' : 'standard', 'pca_components': 6,
    'kn' : 1, 'n_perm': 1000, 'umap_min_dist': 0.1, 'min_cluster_size': 3
    }
   
#feats.multivar_plot(plot_cols = var_cols, hue='SpC', type_shown = 'all', save_output = 'kde_metrics.pdf')

feats_agg, _ = feats.aggregate(group_by = ['STAR','SpC','TIC','RA','DEC'], **ml_kwargs)
#feats_agg[['STAR','SpC','CROWDSAP','Tmag','Q_JHK'] + time_metrics + freq_metrics +rn_metrics].to_csv('feats_agg.csv',index=False)

# SAMPLE PRINTING
#TT = TexTab()
#TT.TabSample(feats.df[['STAR','RA','DEC','SpC','Q_JHK','Tmag','TIC','CROWDSAP']].sort_values('RA'))

#### PAIR PLOTS
#feats.pair_plot(['W0','Tmag'], star_labels = True) # strong (0.84, p<1e-7)
#feats.pair_plot(['R0','Tmag'], star_labels = True) # weaker (0.44, p=0.03)
#feats.pair_plot(['EMSE0','Tmag'], star_labels = True) #uncorrelated (-0.07, p=0.71)
#feats.pair_plot(['EMSE1','Tmag'], star_labels = True) #uncorrelated (0.14, p=0.50)
#feats.pair_plot(['EMSE1','R0'], star_labels = True)  # strong (0.77, p<1e-5)
#feats.pair_plot(['EMSE0','TAU'], star_labels = True)  # strong (-0.65, p<1e-3)


'''fig, ax = plt.subplots(2, 1, sharex =  True, figsize = (8,16)); ax= ax.flatten()
feats.pair_plot(ax = ax[0], pair_c = ['HK','JH'])
feats.pair_plot(ax = ax[1], pair_c = ['HK','Q_JHK'], star_labels = True)
ax[0].set_xticklabels([]); ax[0].set_xlabel('')
fig.subplots_adjust(hspace=0.02); fig.tight_layout()'''


try:
    ml_kwargs['var_cols'].remove('W0')
#ml_kwargs['var_cols'].remove('FSTEPS')
#ml_kwargs['var_cols'].remove('R0')
except:
    pass

#### CLASSIFIERS - REGRESSORS
#feats.knn_classify(merge_cand=True, **ml_kwargs)
#feats.knn_regress(regress_col='Tmag', exclude_labels=[],**ml_kwargs)
#nn_table = feats.nearest_neighbors(**ml_kwargs)#; print(nn_table)
#
#feats.hier_clustering(k = None, **ml_kwargs)
#feats.hdbscan_clustering(**ml_kwargs)
#feats.local_outlier_factor(**ml_kwargs)

#### UMAP PLOTS
fig, ax = plt.subplots()#1, 2, figsize = (18,10)); ax= ax.flatten()
feats.umap_plot(ax=ax, umap_n = 10, show_mutual=True,**ml_kwargs)
#feats.umap_plot(ax=ax[1], umap_n = 15, **ml_kwargs) 
fig.tight_layout();#fig.savefig('UMAP.pdf', format = 'pdf')'''

ml_kwargs['kn'] = 3
knn_tab = feats.nearest_neighbors(show_matrix = True, **ml_kwargs)
#print(knn_tab)

'''fig, ax = plt.subplots(1,3, figsize = (18,6)); ax= ax.flatten()
feats.umap_plot(ax=ax[0], cbar_col = ['Tmag'], umap_n = 6, **ml_kwargs)
feats.umap_plot(ax=ax[1], cbar_col = ['Q_JHK'], umap_n = 6, **ml_kwargs)
feats.umap_plot(ax=ax[2], cbar_col = ['CROWDSAP'], umap_n = 6, **ml_kwargs)
ax[1].set_yticklabels([]); ax[1].set_ylabel('')
ax[2].set_yticklabels([]); ax[2].set_ylabel('')
fig.subplots_adjust(wspace=0.0); fig.tight_layout()'''

