# scatter plot
plot_mcs = {'s' : 25, 'marker' : 's', 'ls': 'None'}
plot_mw = {'s' : 30, 'marker' : '.', 'ls': 'None'}

#plot
plot_LBV  = {'ms' : 13, 'c' : 'k', 'marker' : 'd', 'mew': 0.7, 'ls': 'None', 'mfc' : 'None'}
plot_BREs  = {'ms' : 12, 'c' : 'k', 'marker' : 'p', 'mew': 0.7, 'ls': 'None', 'mfc' : 'None'}
plot_cLBV = {'ms' : 16, 'c' : 'k', 'marker' :'$\u25cc$', 'ls': 'None', 'mfc' : 'None'}
plot_cBRC  = {'ms' : 15, 'c' : 'k', 'marker' : '$\u2b1a$', 'ls': 'None', 'mfc' : 'None'}

PLOT_XLC_NCOL = 1  #2
PLOT_XLC_NROW = 4  #5

PLOT_XLABEL =   { 
		'flux' : r'Time $-$ 2457000 [BTJD d]',
		'nflux' : r'Time $-$ 2457000 [BTJD d]',
        'dmag' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]',
		'sed': r'Wavelength (A)'
		}

PLOT_YLABEL =   {
    'flux' : 'Flux [e-/s]',
    'nflux': 'Normalized flux',
    'dmag' : r'$\Delta$m [mag]',
    'ls' : r'Amplitude (mag)'
    }
        
GAIA_UPMARK = 64

SIZE_FONT_SUB = 12
SIZE_XLABEL_FIG = 22
SIZE_YLABEL_FIG = 22
SIZE_XLABEL_SUB = 11
SIZE_YLABEL_SUB = 11

SIZE_GRID = (26,20)#(16,20)

PLOT_PARAMS =	{
'lc'		:
		{'legend.fontsize': 8,
	 	'font.size':  SIZE_FONT_SUB,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': SIZE_XLABEL_SUB,
         	'ytick.labelsize': SIZE_YLABEL_SUB},
'ls'		:
		{'legend.fontsize': 12,
	 	'font.size':  12,
         	'axes.labelsize': 16,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 17,
         	'ytick.labelsize': 17},
'prew'		:
		{'legend.fontsize': 8,
	 	'font.size':  SIZE_FONT_SUB+6,
         	'axes.labelsize': 24,
         	'axes.titlesize': 13,
         	'xtick.labelsize': SIZE_XLABEL_SUB+5,
         	'ytick.labelsize': SIZE_YLABEL_SUB+5},
'sed'		:
		{'legend.fontsize': 10,
	 	'font.size':  8,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 18,
         	'ytick.labelsize': 11},
'cr'		:
		{'legend.fontsize': 11,
		#'text.usetex' : True,
	 	'font.size':  11,
         	'axes.labelsize': 16,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 15,
         	'ytick.labelsize': 15},
'panel'		:
		{'legend.fontsize': 12,
	 	'font.size':  14,
         	'axes.labelsize': 22,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 23,
         	'ytick.labelsize': 23}
		}
    
CBAR_TITLE_SIZE = 11
CBAR_TICK_SIZE = 10


def styled_label(key):
    
    if 'log' in key:
        
        return r'log$_{10}$(%s)' % STY_LB[key.replace('log','')]
    
    elif 'e_' in key:
        
        return r'var(%s)' % STY_LB[key.replace('e_','')]
    
    else:
        return STY_LB[key]

STY_LB = 	{
		'VSINI' : r'$v$sin $i$','MDOT' : r'$\dot{M}$','LOGLM': r'log(L/M)', 
        'LOGQ' : r'log$_{10}Q$', 'LOGG' : r'log$g$', 'MASS' : '$M_{evol}$', 'VMAC': r'$v_{mac}$ [km s$^{-1}$]', 
        'VMIC' : r'$v_{mic}$ [km s$^{-1}$]', 'NABUN' : r'N/H', 'LOGD' : r'log$_{10}D$', 'S_MASS' : r'$M_{Rg}$',
		'W0' : r'log$_{10}(W_{0})$', 'R0' : r'log$_{10}(R_{0})$', 'TAU' : r'$\tau$', 'GAMMA' : r'$\gamma$' ,
		'TEFF' :  r'T$_{\rm eff}$ [K]', 'SpCt' : 'B(*)I/II', 'TESS_time' : r'Time $-$ 2457000 [BTJD d]','TESS_freq' : r'Frequency [d$^{-1}$]',
        'SLOGL' : r'log$_{10}(\mathcal{L}/\mathcal{L}_{\odot})$', 'LOGL' : r'log$_{10}$($L$/L$_{\odot})$',        
		'MAD' : r'log$_{10}(MAD$)', 'MAD_RAW': r'log$_{10}(MAD_{0})$', 'STD': r'$\sigma$ [mag]', 'ZCROSS' : r'$D_{0}$', 'PSI': r'log$_{10}(\psi^2)$', 'IQR': r'IQR',
		'ETA' : r'$\eta$', 'SKW': r'skw', 'A_V' : r'$A_{V}$', 'EDD' : r'$\Gamma_{e}$', 'KRT': r'kurt', 'ZCR': r'Zcr',
        'EMSE1': r'$m_{E}$', 'EMSE0': r'log$_{10}(\bar{E}^2)$',
        'WFM': r'$\bar{f}_{w}$', 'WFD': r'$\sigma_{f,w}$',
        'HK': r'$H-K_{s}$', 'Q_JHK': r'$Q_{JHK}$', 'JH': r'$J-H$',
        'MSM' : r'MSM', 'MSP' : r'$\bar{{E_s}^2}$', 'MSD' : r'$\sigma_{E_s}$', 'MSC' : r'$\kappa_{E_s}$', 'MSS' : r'$m_{E_s}$',
        'AVECROWD': r'CROWDSAP', 'MINCROWD': r'CROWDSAP', 'Tmag': r'T [mag]', 'RUWE':r'RUWE',
		'MG' : r'$M_{G}$ [mag]', 'MJ' : r'$M_{J}$ [mag]', 'MH' : r'$M_{H}$ [mag]', 'MK': r'$M_{K}$ [mag]', 
        'JK': r'$J-K_{s}$', 'VCHAR' : r'log($\nu_{char}$ [d$^{-1}$])', 'BR': r'$B_{p}-R_{p}$',
		'FF' : r'$f_{i}$ [d$^{-1}$]', 'A_FF' : r'$A_{i}$ [mag]', 'HF' : r'$jf$', 'FFR' : r'$f_{1}/f_{2}$', 
		'A_FFR' : r'$A_{f_{1}}/A_{f_{2}}$', 'BETA' : r'$\beta$', 'VINF' : r'$v_{inf}$ [km s$^{-1}$]',
		'INDFF' : r'#$f_{i}$','INDFFS' : r'#$f_{i,sec}$'}

LC_COLOR = 	{
		'spoc'	 : 'c',
        'tess-spoc' : 'lime',
        'model' : 'r',
        'binned': 'k',
		'tesscut'	 : 'r',
        'fit': 'r',
        'raw': 'pink',
		'any': 'k'
		}

TESS_AP_C = {
    'spoc' : 'r',
    'thres': 'c',
    }


#displayed ID

#binarity status
#merger product
#circumstellar/dust envelope
#circumstellar/molecular disk
#colliding winds & X-ray emission, 
#documented eruption/outburst


STAR_IDS = {
    '6 CAS': [r'6$\,$Cas',1,0,0,0,0,0],
    'AG Car': [r'AG$\,$Car',0,1,1,1,0,1],
    'CD-42 11721': [r'CD-42$\,$11721',1,0,1,1,1,1],
    'CPD-52 9243': [r'CPD-52$\,$9243'],
    'Cyg OB2 12': [r'Cyg$\,$OB2-12'],
    'GG Car': [r'GG$\,$Car'],
    'HD152236': [r'$\zeta1\,$Sco'],
    'HD148937': [r'HD$\,$148937'],
    'HD316285': [r'HD$\,$316285'],
    'HD326823': [r'HD$\,$326823'],
    'HD327083': [r'HD$\,$327083'],
    'HD62623': [r'HD$\,$62623'],
    'HD80077': [r'HD$\,$80077'],
    'HD87643': [r'HD$\,$87643'],
    'HD96918': [r'V382$\,$Car'],
    'HR Car': [r'HR$\,$Car'],
    'HR5171': [r'HR$\,$5171A'],
    'HR8752': [r'V509$\,$Cas'],
    'Hen 3-1383': [r'Hen$\,$3-1383'],
    'Hen 3-519': [r'Hen$\,$3-519'],
    'Hen3-298': [r'Hen$\,$3-298'],
    'IRC+10420': [r'IRC$\,$+10420'],
    'MWC 930': [r'MWC$\,$930'],
    'MWC137': [r'MWC$\,$137'],
    'MWC300': [r'MWC$\,$300'],
    'MWC342': [r'MWC$\,$342'],
    'MWC349': [r'MWC$\,$349'],
    'P Cyg': [r'P$\,$Cyg'],    
    'RHO CAS' : [r'$\rho$ Cas'],
    'V1429 Aql': [r'MWC$\,$314'],
    'V432 Car': [r'WRA$\,$751'],
    'V439 Cyg': [r'V439$\,$Cyg'],
    'WRAY 16-137': [r'WRAY$\,$16-137'],
    'WRAY 16-232': [r'WRAY$\,$16-232'],
    '[B61] 2': [r'WRAY$\,$17-96'],
    '[GKF2010] MN44': [r'MN44'],
    '[GKF2010] MN48': [r'MN48']
    }

cornerplot_kwargs={
    "axes.labelsize" : 14, 
    "legend.title_fontsize": 14,
    "legend.markerscale": 1.2,
    "legend.fontsize" : 12
    }

pairplot_kwargs={
    "axes.labelsize" : 14, 
    "font.size": 6,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,

    "lines.markersize" : 9,
    "legend.title_fontsize": 14,
    "legend.markerscale": 1.2,
    "legend.fontsize" : 12
    }

dendro_kwargs={
    "axes.labelsize" : 14, 
    "axes.titlesize" : 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    }


CLASS_M = {
    'B[e]SG' : 's',
    'LBV' : '^',
    'YHG' : 'o'
    }

CLASS_C = {
    'B[e]SG' : 'b',
    'LBV' : 'g',
    'YHG' : 'r'
    }