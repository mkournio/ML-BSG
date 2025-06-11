import lightkurve as lk
import os

#class MastQuery(object):

#    def __init__(self, tab, download_dir):
#     pass

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

          q = lk._search_products('TIC %s' % row['TIC'], filetype = prod, mission="TESS")
          try:
           q.download_all(download_dir=download_dir)
           print("..fetched %s %s products" % (len(q),prod))
          except SearchError: 
           print('..no %s products were found' % prod)
           pass
           #print('..no %s products were found' % prod) 
     
     return
