from plot_methods import GridTemplate
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from functions import *
from prewhitening import PreWhitening
import lightkurve as lk
import numpy as np
import os 

class ExtractLS(GridTemplate):
	
	def __init__(self, data, prewhiten = False, **kwargs):

		self.data = data
		self.prewhiten = prewhiten
		
		super(ExtractLS,self).__init__(fig_xlabel=PLOT_XLABEL['ls'], fig_ylabel=PLOT_YLABEL['ls'],**kwargs)		
		self._proc()
		self.GridClose()

	def _proc(self):	

		stars_id = self.data['STAR']
		for star in stars_id:
			lc_files = [j for j in os.listdir(TESS_LC_PATH) if (star in j)]
			if len(lc_files) > 0: self._proc_single(star,lc_files)
		return		

	def _proc_single(self, star, lc_files):

		ax = self.GridAx()
		divider = make_axes_locatable(ax)

		g_lc_files = group_consec_lcs(lc_files)	
		len_gls = len(g_lc_files)

		k = 0
		for gl in g_lc_files:

			g_times = np.empty(0); g_fluxes = np.empty(0); g_sects = []
			for l in gl:
				time, flux = lc_read(TESS_LC_PATH + l)
				g_sects.append(sect_from_lc(l))	
				offset = np.nanmedian(flux)
				g_times = np.append(g_times,time)
				g_fluxes = np.append(g_fluxes,flux - offset)
			g_sect_tx = '+'.join([str(x) for x in g_sects])

			# Building light curve from the merged (stitched) sectors
			lc = lk.LightCurve(time=g_times,flux=g_fluxes).remove_nans()

			# Pre-whitening
			if self.prewhiten :
				obj_prw = PreWhitening(lc, star, g_sect_tx, save_files = False, plot_rn = True, rows_page = 6, cols_page = 2, output_format = 'pdf')

			pg = lc.to_periodogram()
			x = np.log10(np.array(pg.frequency))
			y = np.log10(np.array(pg.power))

			ax.plot(x,y,'b')
			ax.set_ylim(-5.9,-2.1)	


			ax.xaxis.set_tick_params(direction="in",labelsize=6)
			ax.yaxis.set_tick_params(labelsize=9)
			if k == 0: ax.text(0.05,0.85,star,color='r',fontsize=10,transform=ax.transAxes)
			ax.text(0.4,0.05,g_sect_tx,color='b',fontsize=9,transform=ax.transAxes	)

			if (len_gls > 1) and (k < len_gls - 1):
				d=.02
				kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
				ax.plot((1,1),(-d,+d), **kwargs) 
				ax.plot((1,1),(1-d,1+d), **kwargs) 
				ax.spines['right'].set_visible(False)

				ax = divider.append_axes("right", size="100%", pad=0.12)

				kwargs.update(transform=ax.transAxes) 
				ax.plot((0,0),(-d,+d), **kwargs) 
				ax.plot((0,0),(1-d,1+d), **kwargs) 
				ax.spines['left'].set_visible(False)
				ax.tick_params(labelleft = False) 
					

			k += 1

		return	
