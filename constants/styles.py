# scatter plot
plot_mcs = {'s' : 25, 'marker' : 's', 'ls': 'None'}
plot_mw = {'s' : 30, 'marker' : '.', 'ls': 'None'}

#plot
plot_LBV  = {'ms' : 13, 'c' : 'k', 'marker' : 'd', 'mew': 0.7, 'ls': 'None', 'mfc' : 'None'}
plot_BREs  = {'ms' : 12, 'c' : 'k', 'marker' : 'p', 'mew': 0.7, 'ls': 'None', 'mfc' : 'None'}
plot_cLBV = {'ms' : 16, 'c' : 'k', 'marker' :'$\u25cc$', 'ls': 'None', 'mfc' : 'None'}
plot_cBRC  = {'ms' : 15, 'c' : 'k', 'marker' : '$\u2b1a$', 'ls': 'None', 'mfc' : 'None'}

PLOT_XLC_NCOL = 2
PLOT_XLC_NROW = 5

PLOT_XLABEL =   { 
		'lc' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]',
		'sed': r'Wavelength (A)'
		}

PLOT_YLABEL =   {
    'flux' : 'Flux [e-/s]',
    'nflux': 'Normalized flux',
    'ls' : r'Amplitude (mag)',
    'dmag' : r'$\Delta$m [mag]'
    }
        
        
        

SIZE_FONT_SUB = 12
SIZE_XLABEL_FIG = 22
SIZE_YLABEL_FIG = 22
SIZE_XLABEL_SUB = 11
SIZE_YLABEL_SUB = 11

SIZE_GRID = (16,20)

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
		'LOGW' : r'log$_{10}$($W$ [mag])', 'LOGR0' : r'log$_{10}$($R_{0}$ [mag])', 'TAU' : r'$\tau$ [d]', 'GAMMA' : r'$\gamma$' ,
		'TEFF' :  r'T$_{\rm eff}$ [K]', 'SpCt' : 'B(*)I/II', 'TESS_time' : r'Time $-$ 2457000 [BTJD d]','TESS_freq' : r'Frequency [d$^{-1}$]',
        'SLOGL' : r'log$_{10}(\mathcal{L}/\mathcal{L}_{\odot})$', 'LOGL' : r'log$_{10}$($L$/L$_{\odot})$',        
		'MAD' : r'MAD', 'STD': r'$\sigma$ [mag]', 'ZCROSS' : r'$D_{0}$', 'PSI': r'$\psi^2$', 'IQR': r'IQR',
		'ETA' : r'$\eta$', 'SKW': r'skw', 'A_V' : r'$A_{V}$', 'EDD' : r'$\Gamma_{e}$', 'KRT': r'kurt', 'ZCR': r'Zcr',
        'MSM' : r'MSM', 'MSP' : r'$\bar{{E_s}^2}$', 'MSD' : r'$\sigma_{E_s}$', 'MSC' : r'$\kappa_{E_s}$', 'MSS' : r'$m_{E_s}$',
        'AVECROWD': r'CROWDSAP', 'MINCROWD': r'CROWDSAP', 'Tmag': r'T [mag]',
		'MG' : r'$M_{G}$ [mag]', 'MJ' : r'$M_{J}$ [mag]', 'MH' : r'$M_{H}$ [mag]', 'MK': r'$M_{K}$ [mag]', 
        'JK': r'$J-K_{s}$', 'VCHAR' : r'log($\nu_{char}$ [d$^{-1}$])',
		'FF' : r'$f_{i}$ [d$^{-1}$]', 'A_FF' : r'$A_{i}$ [mag]', 'HF' : r'$jf$', 'FFR' : r'$f_{1}/f_{2}$', 
		'A_FFR' : r'$A_{f_{1}}/A_{f_{2}}$', 'BETA' : r'$\beta$', 'VINF' : r'$v_{inf}$ [km s$^{-1}$]',
		'INDFF' : r'#$f_{i}$','INDFFS' : r'#$f_{i,sec}$'}

LC_COLOR = 	{
		'spoc'	 : 'c',
        'tess-spoc' : 'lime',
        'binned': 'k',
		'tesscut'	 : 'r',
        'fit': 'k',
        'raw': 'pink',
		'any': 'k'
		}

TESS_AP_C = {
    'spoc' : 'r',
    'thres': 'c',
    }