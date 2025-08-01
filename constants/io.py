# -*- coding: utf-8 -*-

tess_header_keys_p = ['TSTART','TSTOP','DATE-OBS','DATE-END','PROCVER',\
                      'OBJECT','TICID','SECTOR','CAMERA','CCD','RA_OBJ',\
                      'DEC_OBJ','PMRA','PMDEC','PMTOTAL','INSTRUME','PXTABLE','TESSMAG']
tess_header_keys_b = ['BJDREFI','TIMEDEL','CROWDSAP','FLFRCSAP','PDCVAR']



flux_units = {
    'time': 'd',
    'flux': 'e-/s',
    'flux_err': 'e-/s',
    'fitmodel': 'e-/s',
    'nflux': '',
    'nflux_err': '',
    'dmag': 'mag',
    'dmag_err': 'mag',
    }