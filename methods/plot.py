from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib import patches
from constants.styles import *
from methods.tools import check_header_key
import os
from astropy.visualization import PercentileInterval, ImageNormalize, LinearStretch
import numpy as np
from astropy.io import fits
import warnings
import lightkurve as lk

def plot_lc_single(lc, 
                   ax = None,
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

def plot_lc_multi(lc, ax = None, **kwargs):
    
    return

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

def add_plot_features(ax,mode = 'flux',upper_left='',lower_left='',lower_right=''):
    
    ax.text(0.05,0.85,upper_left,color='r',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)
    ax.text(0.05,0.05,lower_left,color='b',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)
    ax.text(0.6,0.05,lower_right,color='b',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)
  
    if mode == 'dmag': ax.invert_yaxis()

    return ax
   
    
class GridTemplate(object):

	# Class for the creation and management of plotting grid
    def __init__(self, rows_page = 3, cols_page = 1, output_format = 'pdf', params = PLOT_PARAMS['lc'], inter = False, figsize = SIZE_GRID, **kwargs):

        plt.rcParams.update(params)
        
        self.__dict__.update(**kwargs)
        self.rows_page = rows_page
        self.cols_page = cols_page
        self.output_format = output_format
        self.inter = inter
        self.figsize = figsize
        
        self.filename = kwargs.pop('plot_name','GridPlot')
        self.fig_xlabel  = kwargs.pop('fig_xlabel','X LABEL')
        self.fig_ylabel  = kwargs.pop('fig_ylabel','Y LABEL')
        self.coll_x = kwargs.pop('coll_x',False)
        self.coll_y = kwargs.pop('coll_y',False)

        self.ind = 0
        
    def GridAx(self):
        
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
                ax.set_xlabel(STY_LB[self.col_labels[plot_col]])
                
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
        
        self.ind += 1
        
        return ax
   
    def _get_filename(self):
        
        fl_id = 0 
        while os.path.exists('%s_%s.%s' % (self.filename,fl_id,self.output_format)): 
            fl_id += 1
            
        return '%s_%s' % (self.filename,fl_id)
    
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
            filename = '%s.%s' % (self._get_filename(),self.output_format)
            self.fig.savefig('%s.%s' % (self._get_filename(),self.output_format), format = self.output_format,bbox_inches='tight')
            
        return