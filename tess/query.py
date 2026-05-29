import lightkurve as lk
from lightkurve.correctors import download_tess_cbvs
from constants.paths import path_to_output_fits
from astropy.io import fits
from methods.tools import get_hdu_from_keys
import shutil
import os

def mast_query(table, 
               product = "Lightcurve", 
               download = True, 
               download_dir='data/'):

     if not os.path.isdir(download_dir):
         os.mkdir(download_dir)
         
     if not set(['STAR','TIC']).issubset(table.columns):
         raise KeyError("Table does not contain STAR and/or TIC columns")
         
     for row in table:
       
       print('Querying from MAST star %s with TIC %s' % (row['STAR'],row['TIC']))
       
       if product == "Lightcurve":
           try:
               q = lk.search_lightcurve('TIC %s' % row['TIC'], mission="TESS")
               q = q[(q.author == 'SPOC') | (q.author == 'TESS-SPOC')]
               if download:
                   q.download_all(download_dir=download_dir)
                   print('Lightcurve: Fetched {} product(s)'.format(len(q)))
                   
           except Exception as e:
               print('Lightcurve: {}'.format(e))
               
       elif product == "Targetpixelfile":
           try:
               q = lk.search_targetpixelfile('TIC %s' % row['TIC'], mission="TESS")
               q = q[(q.author == 'SPOC') | (q.author == 'TESS-SPOC')]
               print(q)               
           except Exception as e:
               print('Targetpixelfile: {}'.format(e))
     
     return
 
def download_tpfs(table,
              frame = 0,
              download_dir = 'data/TPFs/',
              del_original = False):
    
    if not os.path.exists(download_dir+'cache/'):
        os.mkdir(download_dir+'cache/')
    
    for l in table:
        
        star = l['STAR']
        tic = l['TIC']

        filename = [f for f in os.listdir(path_to_output_fits) if star in f]
        if len(filename) > 0:
            
            print(f'Downloading TPF for {star} - TIC {tic}')
            
            hdulist = fits.open(os.path.join(path_to_output_fits,filename[0]))
            hdu_raw = get_hdu_from_keys(hdulist, HDUTYPE = 'LIGHTCURVE', BINNING = 'F')
            
            for h in hdu_raw:
                
                sector = h.header['SECTOR']
                author = h.header['PIPELINE'].upper()

                print(f'- SECTOR {sector}')
                tab = lk.search_targetpixelfile(f'TIC {tic}', 
                                                sector = sector, author=author)
                t = tab[0].download(download_dir=download_dir+'cache/')
                if isinstance(frame,int):                    
                    fits_name = f'{tic}_{sector}_{frame}.fits'
                    if not os.path.exists(download_dir+fits_name):
                        t[frame].to_fits(download_dir+fits_name)         
    if del_original:
        shutil.rmtree(download_dir+'cache/')
    
    return
    

def download_cbvs(table, 
                  cbv_type = "SingleScale", 
                  download_dir = 'data/CBVs/', 
                  band = None):
    
    for l in table:
        
        star = l['STAR']
        filename = [f for f in os.listdir(path_to_output_fits) if star in f]
        if len(filename) > 0:
            
            print(f'Downloading CBVs for {star}')
            
            hdulist = fits.open(os.path.join(path_to_output_fits,filename[0]))
            hdu_raw = get_hdu_from_keys(hdulist, HDUTYPE = 'LIGHTCURVE', BINNING = 'F')

            for h in hdu_raw:
                
                    sector = h.header['SECTOR']
                    camera = h.header['CAMERA']
                    ccd = h.header['CCD']
                    
                    print(f'- SECTOR {sector} CAMERA {camera} CCD {ccd}')                   
                    cbvs = download_tess_cbvs(sector=sector, camera=camera,
                                              ccd=ccd, cbv_type=cbv_type,
                                              band=band, save_dir=download_dir)
    
    
    return