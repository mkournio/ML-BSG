
plot_all = {'ms' : 6, 'c' : 'c', 'marker' : '.', 'ls': 'None'}
plot_LBV  = {'ms' : 14, 'c' : 'k', 'marker' : 'o', 'mew': 1.5, 'ls': 'None', 'mfc' : 'None'}
plot_cLBV = {'ms' : 16, 'c' : 'k', 'marker' :'$\u25cc$', 'ls': 'None', 'mfc' : 'None'}
plot_BRC  = {'ms' : 13, 'c' : 'k', 'marker' : 's', 'mew': 1.5, 'ls': 'None', 'mfc' : 'None'}
plot_cBRC  = {'ms' : 15, 'c' : 'k', 'marker' : '$\u2b1a$', 'ls': 'None', 'mfc' : 'None'}

PLOT_XLC_NCOL = 2

PLOT_XLC_NROW = 5

PLOT_XLABEL =   { 
		'lc' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]',
		'sed': r'Wavelength (A)'
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

STY_LB = 	{
		'VSINI' : r'$v$sin $i$','MDOT' : r'$\dot{M}$',
		'LOGLM': r'log(L/M)', 'LOGL' : r'log(L/L$_{\odot}$)', 'LOGQ' : r'log$Q$',
		'LOGG' : r'log$g$', 'MASS' : '$M_{evol}$', 'VMAC': r'$v_{mac}$ [km s$^{-1}$]', 'VMIC' : r'$v_{mic}$ [km s$^{-1}$]',
		'NABUN' : r'N/H', 'LOGD' : r'log$D$', 'S_MASS' : r'$M_{Rg}$',
		'LOGW' : r'log($W$ [mag])', 'LOGR0' : r'log($R_{0}$ [mag])', 'TAU' : r'$\tau$ [d]', 'GAMMA' : r'$\gamma$' ,
		'TEFF' :  r'log(T$_{\rm eff}$ [K])', 'TESS_time' : r'Time $-$ 2457000 [BTJD d]',
		'TESS_freq' : r'Frequency [d$^{-1}$]','S_LOGL' : r'log($L$/L$_{\odot})$',
		'MAD' : r'MAD', 'SVAR': r'$\sigma$ [mag]', 'ZCROSS' : r'$D_{0}$', 'PSI': r'$\psi^2$',
		'ETA' : r'log($\eta$)', 'SKEW': r'skw', 'A_V' : r'$A_{V}$', 'EDD' : r'$\Gamma_{e}$',
		'GABS' : r'$M_{G}$ [mag]', 'VCHAR' : r'log($\nu_{char}$ [d$^{-1}$])',
		'FF' : r'$f_{i}$ [d$^{-1}$]', 'A_FF' : r'$A_{i}$ [mag]', 'HF' : r'$jf$', 'FFR' : r'$f_{1}/f_{2}$', 
		'A_FFR' : r'$A_{f_{1}}/A_{f_{2}}$', 'BETA' : r'$\beta$', 'VINF' : r'$v_{inf}$ [km s$^{-1}$]',
		'INDFF' : r'#$f_{i}$','INDFFS' : r'#$f_{i,sec}$'}

LC_COLOR = 	{
		'spoc'	 : 'c',
        'spoc_binned': 'b',
		'tesscut'	 : 'r',
        'fit': 'k',
		'any': 'k'
		}

TESS_AP_C = {
    'spoc' : 'r',
    'thres': 'c',
    }