from astropy import units as u
from astroquery.xmatch import XMatch
from astropy.table import Table,  hstack, vstack, unique, join
import numpy as np
from astropy.io import ascii
#from constants import *
from methods import *
from args import *
import pandas as pd

def read_viztab(cat):

        from astroquery.vizier import Vizier
        
        vizier = Vizier(columns=["_RA","_DE"])
        vizier.ROW_LIMIT = -1
        
        return vizier.get_catalogs(cat.split(':')[1])[0]



class PostProcess(Table):

        def __init__(self, ini_tab = None,  xrad = 2, **kwargs):

                        self.xcats = XM_CATS
                        self.db_keys = DOUBLE_KEYS
                        self.xrad = xrad                        
                        super().__init__(ini_tab,**kwargs)

             
        @update_table                     
        @remove_doubles                                                
        def get_xmtab(self,cat,cols):

                        xm_tab = XMatch.query(cat1 = self, cat2 = cat, max_distance = self.xrad * u.arcsec, colRA1 = 'RA', colDec1 = 'DEC')     

                        return xm_tab[cols]

        def xmatch(self):
                        jtab = self.copy()
                        
                        for xcat in jtab.xcats:  
                         xcol = jtab.db_keys + jtab.xcats[xcat]                
                         tab = jtab.get_xmtab(xcat,xcol)
                         jtab = join(jtab, tab, join_type='left', keys= jtab.db_keys)
                        
                        return jtab
                        
        @append_col                
        def append(self,cat='combined'):
                        app_tab = self.copy()
                        return app_tab
                         

class PreProcess(object):

        def __init__(self, args):
                        self.cats = args.INI_CATS
                        
        @update_table
        def read_xtab(self,cat):
                        return ascii.read(cat)
        @update_table
        def read_vtab(self,cat): 
                        return read_viztab(cat)
                       
        @remove_doubles                  
        def compile(self):

                        tabs = []	
                        for cat in self.cats:
                         if cat.startswith('vizier:'): 
                          tab = self.read_vtab(cat)
                         else:
                          tab = self.read_xtab(cat)

                         tabs.append(tab)
                         
                        vstab = vstack(tabs)
                        if len(MERGE_KEYS) > 0:
                          for mk in MERGE_KEYS:
                           vstab[mk] = merged_col(vstab[DOUBLE_KEYS],vstab[mk])
                        
                        return vstab

