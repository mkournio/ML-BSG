import lightkurve as lk
import os

def mast_query(table,download_dir='data/',product="Lightcurve"):

     if not os.path.isdir(download_dir):
         os.mkdir(download_dir)
         
     if product not in ["ffi","Lightcurve"]:
         raise KeyError("Unrecognized filetype")
         
     if not set(['STAR','TIC']).issubset(table.columns):
         raise KeyError("Table does not contain STAR and/or TIC columns")
         
     t_ind = 0
     for row in table:
       
       print('Querying from MAST star %s with TIC %s' % (row['STAR'],row['TIC']))
       
       if product == "Lightcurve":       
          try:
           q = lk.search_lightcurve('TIC %s' % row['TIC'], mission="TESS")
           q = q[(q.author == 'SPOC') | (q.author == 'TESS-SPOC')]
           q.download_all(download_dir=download_dir)
           print('Lightcurve: Fetched {} product(s)'.format(len(q)))
           if len(q)>0: t_ind +=1

          except Exception as e:
           print('Lightcurve: {}'.format(e))
           
     print('Data available for {} objects'.format(t_ind))
     
     return
