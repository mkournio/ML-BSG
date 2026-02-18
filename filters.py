

#### FILTER FOR MATIAS
#MAT= (cm['RA'] > 80) & (cm['RA'] < 180) & (cm['DEC'] < 10) & (cm['DEC'] > -70)  & (cm['GAL'] == 'MW')  & (cm['Gmag'] < 10) 
#cm = cm[MAT]
#print(cm['STAR','SpC','GAL','RA','DEC','Gmag',].pprint(max_lines=-1,max_width=-1))
#ra_c,dec_c = np.genfromtxt('casleo_mat2',delimiter=',',unpack=True,usecols=(0,1),comments='#')
#print(ra_c)
#mask = [float('%.4f'% x) in ra_c for x in cm['RA']] 
#cm = cm[mask]
#print(cm['STAR','RA','DEC','SpC','GAL','Gmag','CROWDSAP'].pprint(max_lines=-1,max_width=-1))
#Visualize(data=cm, plot_name='MATIAS', plot_key='dmag', rows_page=7, cols_page=1,join_pages=False, output_format='png').lightcurves(stitched=True, bin_size = '10m')
#################

#### FILTER FOR CASLEO 2026A PROPOSAL
#LOC = [('MW' in x) for x in cm['GAL']]
#cm = cm[LOC]
#CASLEO = (100<cm['RA']) & \
#         (cm['RA']<300) & \
#         (cm['DEC']<0) & \
#         (cm['BPmag']<11) & (cm['Gmag']<11)
#cm = cm[CASLEO]
#ra_c,dec_c = np.genfromtxt('casleo_pv',unpack=True,usecols=(0,1),comments='#')
#mask = [float('%.4f'% x) in ra_c for x in cm['RA']] 
#cm = cm[mask]
#####################################

class Filters(object):      
           
        @staticmethod           
        def EMS(cm):
           
          
           EMS = [('LBV' in x) or ('B[e]SG' in x) or ('YHG' in x) for x in cm['SpC']]           
           cm = cm[EMS]
           
           return 
       
       

