from args import *
import numpy as np
from astropy.table import unique

def update_table(func):
         return combine_procs(append_col, rename_col)(func)
 
def combine_procs(*procs):

         def cproc(func):
         
          for proc in reversed(procs):
            func = proc(func)
          return func
          
         return cproc            
         
def append_col(ftab):

         def wrapper(*args):

           tab = ftab(*args).copy()
           if APPEND_COLS == {} or len(args) < 2:
            return tab
           else:
            for cat, col_dir in APPEND_COLS.items():
             if cat == args[1]:
              for ncol, nval in col_dir.items():              
               if callable(nval):
                try:
                 tab[ncol] = nval(tab)
                except KeyError:
                 tab[ncol] = np.nan
               else:
                tab[ncol] = nval   

            return tab
                  
         return wrapper

def rename_col(ftab):

         def wrapper(*args):

           tab = ftab(*args).copy()
           if RENAME_COLS == {}:
            return tab
           else:
            for ncol, pcols in RENAME_COLS.items():
             for pc in pcols:
              if pc in tab.columns: 
               tab[pc].name = ncol 
                      
            return tab
                  
         return wrapper
         
def remove_doubles(ftab):

         def wrapper(*args):
         
          tab = ftab(*args).copy() 
          untab = unique(tab, keys=DOUBLE_KEYS)
          
          return untab
          
         return wrapper
         

