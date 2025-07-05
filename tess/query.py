import lightkurve as lk
import os

def mast_query(table,download_dir='data/',products=["ffi","Lightcurve"]):

     if not os.path.isdir(download_dir):
         os.mkdir(download_dir)
         
     if not set(products).issubset(["ffi","Lightcurve"]):
         raise KeyError("Unrecognized filetype")
         
     if not set(['STAR','TIC']).issubset(table.columns):
         raise KeyError("Table does not contain STAR and/or TIC columns")
    
     for row in table:
       
       print('Querying from MAST star %s with TIC %s' % (row['STAR'],row['TIC']))
       
       for prod in products:
       
          try:
           q = lk._search_products('TIC %s' % row['TIC'], filetype = prod, mission="TESS")
           if prod == "Lightcurve":
             q = q[(q.author == 'SPOC') | (q.author == 'TESS-SPOC')]
           q.download_all(download_dir=download_dir)
           print('{}: Fetched {} product(s)'.format(prod,len(q)))
          except Exception as e:
           print('{}: {}'.format(prod,e))
     
     return
