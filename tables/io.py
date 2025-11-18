#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:48:29 2025

@author: michalis
"""

from astropy.io import ascii
import numpy as np

def tab_to_csv(table,output,**kwargs):
    
    ctab = table.copy().to_pandas()
    ctab.to_csv(output)
    
    return
