import args 
from tables import *
from astropy.io import ascii
import numpy as np
from astropy.table import Table, Column

'''class SelfTab(Table):
        def __init__(self):
         super().__init__()
         
        def __getitem__(self,col):
         try:
          return super().__getitem__(col)
         except KeyError:
          return Column()
          
         
a = SelfTab()
print(a['c'])
col_c = Column(name='c', data=['x', 'y', 'z'])
col_d = Column(name='d', data=['u', 'v', 'w'])
t2=Table()
t2.add_columns([col_c, col_d])
a.add_columns(list(t2.columns.values()))
print(a['f'])
#print(a['TEFF'])
'''
XC = Compiler(args)
print(XC.compiled['STAR','RA','DEC','BETA','NABUN','TEFF','e_TEFF','E_TEFF','LOGG','e_LOGG','E_LOGG'].pprint(max_lines=-1))

