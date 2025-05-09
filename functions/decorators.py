from args import *
import numpy as np

def update_table(func):
         return combine_procs(new_col, trans_col)(func)
         
def combine_procs(*procs):
    def cproc(func):
        for proc in reversed(procs):
            func = proc(func)
        return func
    return cproc
            
         
def new_col(ftab):
        
         def wrapper(self,fcat):
    
           tab = ftab(self,fcat).copy()
           if NEW_COLS == {}:
            return tab
           else:
            for cat, col_dir in NEW_COLS.items():
             if cat == fcat:
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

def trans_col(ftab):
        
         def wrapper(self,*args):
    
           tab = ftab(self,*args).copy()
           if TR_COLS == {}:
            return tab
           else:
            for ncol, pcols in TR_COLS.items():
             for pc in pcols:
              if pc in tab.columns: 
               tab[pc].name = ncol        
            return tab
                  
         return wrapper
         
