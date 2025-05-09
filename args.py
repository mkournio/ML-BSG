INI_CATS = ['vizier:J/A+A/687/A228/tabled1', 'FRASER+10', 'HAUCKE+19'] 
          
TR_COLS = {
          'TEFF' : ['Teff'],
          'e_TEFF' : ['e_Teff'],
          'E_TEFF' : ['E_Teff'], 
          'LOGG' : ['logg'],
          'e_LOGG' : ['e_logg'],  
          'E_LOGG' : ['E_logg'],           
          'RA'   : ['_RAJ2000'],
          'DEC'  : ['_DEJ2000'],
          'STAR' : ['Name'],
          'VSINI': ['vsini'],
          'SpC'  : ['SPTYPE']
          }
          
NEW_COLS = {
          'FRASER+10' : {'e_TEFF': 1000, 'E_TEFF': 1000, 'e_LOGG': 0.1, 'E_LOGG': 0.1},
          'HAUCKE+19' : {'E_TEFF': lambda x: x['e_TEFF'], 'E_LOGG': lambda x: x['e_LOGG']}
	  }


xm_cols = {
          'vizier:IV/38/tic' : ['TIC'],
          }
        


