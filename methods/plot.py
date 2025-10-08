from pypdf import PdfWriter
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib import patches
from constants.styles import *
from methods.tools import check_header_key
from methods.functions import round_to
import os
from astropy.visualization import PercentileInterval, ImageNormalize, LinearStretch
import numpy as np
from astropy.io import fits
import warnings
import lightkurve as lk

def plot_lc_single(ax, 
                   lc,
                   m = '',
                   flux_key ="flux",
                   lc_type = 'any',
                   **kwargs):
    
    if lc is None:
        return
    
    if ax is None:
     _, ax = plt.subplots()    
     
    if isinstance(lc, lk.LightCurve):
        time = lc.time.value
        flux = lc[flux_key].value
        try:
            flux_err = lc[flux_key+"_err"].value
        except:
            flux_err = np.zeros(len(flux)) 
    elif isinstance(lc, fits.BinTableHDU):
        time = lc.data['time']
        flux = lc.data[flux_key]
        try:
            flux_err = lc.data[flux_key+"_err"]
        except:
            flux_err = np.zeros(len(flux)) 
    elif isinstance(lc, np.ndarray):
        time = lc[0]
        flux = lc[1]
        if lc.shape[0] > 2:
            flux_err = lc[2]
    else:
        raise TypeError('Object light curve does not have supportive format!')

    ax.plot(time,flux,m,c = LC_COLOR[lc_type])
    
    return ax

def plot_lc_multi(axes,
                  grouped_hdus,
                  m = '',
                  flux_key ="flux",
                  lc_type = 'any',
                  **kwargs):
    
    if not isinstance(axes,list) and len(axes) != len(grouped_hdus):
        raise TypeError('Size of grouped HDUs does not match that of axes.')
   
    for g_hdu, ax in zip(grouped_hdus,axes):
        
        g_times = np.empty(0); g_fluxes = np.empty(0); g_sects = []
        for hdu in g_hdu:
            
            if flux_key == 'dmag':
                offset = np.nanmedian(hdu.data[flux_key])
            else:
                offset = 0.
            g_times = np.append(g_times,hdu.data['time'])
            hdu.data[flux_key][-1] = np.nan
            g_fluxes = np.append(g_fluxes,hdu.data[flux_key]-offset)
            g_sects.append(hdu.header['SECTOR'])
 
        if len(g_sects) > 1: 
            g_sect_tx = r'{}$-${}'.format(g_sects[0],g_sects[-1])
        else:
            g_sect_tx = r'{}'.format(g_sects[0])
        ax.plot(g_times, g_fluxes,m,c = LC_COLOR[lc_type])
        ax.text(0.48,0.05,g_sect_tx,color='k',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)
        
        x1_p = (2*g_times[0]+g_times[-1])/3.
        x2_p = (2*g_times[-1]+g_times[0])/3.
        ax.set_xticks([round_to(x1_p,5)[0], round_to(x2_p,5)[1]])
        
    return axes   

def plot_tess_field(field, ax = None, spoc_aperture = None, thr_aperture = None, **kwargs):
    
    # Plots TESS field from ImageHDU or ndarray
    
    frow, fcol = 0, 0
    if isinstance(field,fits.ImageHDU):
     frow = field.header['2CRV5P']
     fcol = field.header['1CRV5P']
     field = field.data
    elif isinstance(field,np.ndarray):
     #TODO
     pass 
    else:
     raise TypeError('Object field should be either ImageHDU or ndarray!')
 
    extent = (fcol - 0.5, fcol + field.shape[1] - 0.5, 
              frow - 0.5, frow + field.shape[0] - 0.5)

    # Based on the LightKurve plot_image function  
    if ax is None:
     _, ax = plt.subplots()    
    
    mask = np.nan_to_num(field) > 0
    if mask.any() > 0:
      vmin, vmax = PercentileInterval(95.0).get_limits(field[mask])
    else:
      vmin, vmax = 0, 0
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch(), clip=False)

    ax.imshow(field, origin='lower', norm=norm, extent = extent, **kwargs)
    
    if spoc_aperture is not None:
        if isinstance(spoc_aperture,fits.ImageHDU):
            if check_header_key(spoc_aperture,'EXTNAME','APERTURE'):
                ref_row, ref_col = spoc_aperture.header['CRVAL2P'], spoc_aperture.header['CRVAL1P']
                ap_data = spoc_aperture.data
                
                # Convert TESS flags to boolean mask
                ap_data = [[np.binary_repr(x, width = 8)[-2] == '1' for x in row] for row in ap_data]
                
                in_aperture = np.where(ap_data)
                ap_row = in_aperture[0] + ref_row - 0.5
                ap_col = in_aperture[1] + ref_col - 0.5
                for ii in range(len(ap_row)):
                    rect=patches.Rectangle((ap_col[ii],ap_row[ii]),1,1, fill=False, hatch="//", color=TESS_AP_C['spoc'])
                    ax.add_patch(rect)
            else:
                warnings.warn('Invalid SPOC aperture provided.')
        else:
            warnings.warn('SPOC aperture has to be of fits.ImageHDU type.')
                
     
                
    #ax.set_xlabel(xlabel)
    #ax.set_ylabel(ylabel)
    #ax.set_title(title)
   # if show_colorbar:
    #    cbar = plt.colorbar(cax, ax=ax, label=clabel)
     #   cbar.ax.yaxis.set_tick_params(tick1On=False, tick2On=False)
      #  cbar.ax.minorticks_off()
        
    return ax

def add_plot_features(ax,mode = 'flux',upper_left='',lower_left='',lower_right='', upper_right='', y_min_max = None):
     
    if not isinstance(ax,list):
        ax=[ax]
                     
    if y_min_max is not None:
        ymin, ymax = y_min_max
        r = 0.4 * (ymax - ymin)
        for ax_d in ax:
            ax_d.set_ylim(ymin-r, ymax+r)            
  
    if mode == 'dmag':
        for ax_d in ax:
                     ax_d.invert_yaxis()            
            
    ax[0].text(0.05,0.85,upper_left,color='r',fontsize=SIZE_FONT_SUB,transform=ax[0].transAxes)
    ax[0].text(0.05,0.05,lower_left,color='b',fontsize=SIZE_FONT_SUB,transform=ax[0].transAxes)
    ax[-1].text(0.6,0.05,lower_right,color='b',fontsize=SIZE_FONT_SUB,transform=ax[-1].transAxes)
    ax[-1].text(0.6,0.85,upper_right,color='k',fontsize=SIZE_FONT_SUB,transform=ax[-1].transAxes)  
    
    return ax

def get_filename(fname,fformat):
        
        fl_id = 0 
        while os.path.exists('%s_%s.%s' % (fname,fl_id,fformat)): 
            fl_id += 1
            
        return '%s_%s.%s' % (fname,fl_id,fformat)
    
class GridTemplate(object):

	# Class for the creation and management of plotting grid
    def __init__(self,
                 rows_page = 3,
                 cols_page = 1,
                 output_format = 'pdf',
                 join_pages = False,
                 params = PLOT_PARAMS['lc'],
                 figsize = SIZE_GRID,
                 inter = False,
                 **kwargs):

        plt.rcParams.update(params)
        
        self.__dict__.update(**kwargs)
        self.rows_page = rows_page
        self.cols_page = cols_page
        self.output_format = output_format
        self.inter = inter
        self.figsize = figsize
        self.join_pages = join_pages
        
        self.filename = kwargs.pop('plot_name','GridPlot')
        self.fig_xlabel  = kwargs.pop('fig_xlabel','X LABEL')
        self.fig_ylabel  = kwargs.pop('fig_ylabel','Y LABEL')
        self.coll_x = kwargs.pop('coll_x',False)
        self.coll_y = kwargs.pop('coll_y',False)
        
		#Initiating pdf book
        if self.join_pages:
            self.output_format = 'pdf'
            self.pdf_list = []

        self.ind = 0
        
    def GridAx(self,
               divide=False,
               ax_scaling=[]):
        
        # Returns the position axis on the grid when called from the child class

        plot_row = int((self.ind / self.cols_page) % self.rows_page)
        plot_col = int(self.ind % self.cols_page)
        plot_tot = (self.cols_page*self.rows_page)
        
        if (self.ind % plot_tot  == 0):
            
            if self.ind > 0:
                self._page_close()	# Close/Save the existing grid
                
            self._grid_new()

        ax = self.fig.add_subplot(self.gs[plot_row,plot_col])
        self.ax_pos = plot_row
        
        if plot_row < self.rows_page - 1:
            if self.coll_x:
                ax.tick_params(labelbottom=False)
                self.gs.update(hspace=0.05)
        else:
            if 'col_labels' in self.__dict__:
                ax.set_xlabel(STY_LB[self.col_labels[pl_ot_col]])
                
        if plot_row == 0 and 'sup_xlabels' in self.__dict__:
            ax.set_title(self.sup_xlabels[plot_col])
            
        if plot_col > 0:
            if self.coll_y:
                ax.tick_params(labelleft=False)
                self.gs.update(wspace=0.05)
        else:
            if 'row_labels' in self.__dict__:
                ax.set_ylabel(STY_LB[self.row_labels[plot_row]])	
                
        if 'xlim' in self.__dict__: ax.set_xlim(self.xlim)
        if 'ylim' in self.__dict__: ax.set_ylim(self.ylim)
        
        if divide:
            ax = self._axis_divide(ax,ax_scaling)
        
        self.ind += 1
        
        return ax
    
    def _axis_divide(self,ax,ax_scaling,**kwargs):
        
        ax_v = [ax]        
        if len(ax_scaling) == 0:            
            return [ax]
        
        d  = kwargs.pop('d', 0.02)
        divider = make_axes_locatable(ax)
        for xs in ax_scaling:                      

            axargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            ax.plot((1,1),(-d,+d), **axargs)
            ax.plot((1,1),(1-d,1+d), **axargs)
            ax.spines["right"].set_visible(False)
            
            ax = divider.append_axes("right", size="%s%%" % xs, pad=0.14)
            
            axargs.update(transform=ax.transAxes) 
            ax.plot((0,0),(-d,+d), **axargs)
            ax.plot((0,0),(1-d,1+d), **axargs)
            ax.spines["left"].set_visible(False)
            ax.tick_params(labelleft = False)
            
            ax_v.append(ax)            
        
        return ax_v
  
    def _grid_new(self):
        
        #Create grid on new page
        
        self.fig = plt.figure(figsize=self.figsize)
        self.gs = GridSpec(self.rows_page, self.cols_page, figure=self.fig)
        
        return
    
    def _page_close(self):
        
        #Close/Save grid in page
        
        self._save_output()
        if self.inter: 
            plt.show()
        
        self.fig.clf()
        plt.close(self.fig)
        
        return
    
    def _save_output(self):
        
        self.glob_ax = self.fig.add_subplot(self.gs[:self.ax_pos+1, :], frameon=False)
        self.glob_ax.tick_params(labelleft = False, labelbottom = False, bottom = False, left = False)
        self.glob_ax.set_xlabel(self.fig_xlabel, fontsize = SIZE_XLABEL_FIG, labelpad=30)
        self.glob_ax.set_ylabel(self.fig_ylabel, fontsize = SIZE_YLABEL_FIG, labelpad=55)
        
        if self.output_format != None:
            
            file = get_filename(self.filename,self.output_format)
            self.fig.savefig(file, format = self.output_format, bbox_inches='tight')
            if self.join_pages:
                self.pdf_list.append(file)
                 
        return
    
    def close_plot(self):
        
        if hasattr(self,'fig'):
            self._page_close()
            
        if self.join_pages:
            merger = PdfWriter()
            for f in self.pdf_list:
                merger.append(f)
            merger.write(get_filename(self.filename+'_m','pdf'))
            merger.close()
            for r in self.pdf_list: 
                os.remove(r)
            
        return
        
        
        
        