from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
import time as time
#import concurrent.futures
from astropy.table import Table

####opening the catalog I want to obtain more information
print("Opening your main catalog")
data_my = Table.read('your_cool_main_catalog.fits'
                        ,format='fits')

#####selecting the columns you want to select from your main catalog
print("selecting columns from the main catalog")
####add your columns from the main catalog here
QUICK_OBJECT_ID = np.array(data_my['QUICK_OBJECT_ID'])
RA = np.array(data_my['RA'])
DEC = np.array(data_my['DEC'])


#######opening the catalog that has additional information
print("Opening the catalog with additional information")
data_add = Table.read('you_cool_ancillary_data.fits', format='fits')
df_add = data_add.to_pandas()
df_add.info()
print(len(data_add))

print("selecting columns from the catalog with additional information")
RA_big = np.array(data_add['ra'])
DEC_big = np.array(data_add['dec'])

print("Starting the crossmatch")
inicio = time.time()
print("running the match")
c_my = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
c_add = SkyCoord(ra=RA_big*u.degree, dec=DEC_big*u.degree)
idx, d2d, d3d = c_my.match_to_catalog_sky(c_add)
max_separation = 1.0*u.arcsec
sep_constraint = d2d < max_separation
final = time.time()
print("time to peform the match:",final - inicio, "seconds")
print("time to perform the match:", (final - inicio)/60, "minutes")
print("Size of the match cataog:",len(idx))

print("selecting objects that are trully matched")
if len(idx[sep_constraint]) > 0:
    ####these are columns from the main catalog
    QUICK_OBJECT_ID_match = QUICK_OBJECT_ID[sep_constraint]
    QUICK_OBJECT_ID_match_f = QUICK_OBJECT_ID_match.reshape(-1,1)
    RA_match = RA[sep_constraint]
    RA_match_f = RA_match.reshape(-1,1)
    DEC_match = DEC[sep_constraint]
    DEC_match_f = DEC_match.reshape(-1,1)
    
    #####these are columns from the addtional information catalog
    RA_big_match = RA_big[idx[sep_constraint]]
    RA_big_match_f = RA_big_match.reshape(-1,1)
    DEC_big_match = DEC_big[idx[sep_constraint]]
    DEC_big_match_f = DEC_big_match.reshape(-1,1)
    
    ####saving all the columns into a single array
    con = np.concatenate((QUICK_OBJECT_ID_match_f,RA_match_f,DEC_match_f,RA_big_match_f,DEC_big_match_f,), axis=1)
         
    ###saving the con array  into a table
    df_out = pd.DataFrame(con, columns= ['QUICK_OBJECT_ID','RA','DEC','ra_add','dec_add'])
    df_out.info()
    t = Table.from_pandas(df_out)
    t.write('your_match_catalog.fits', overwrite=True) 
else:
    print("Not enough objects inside the crossmatch Radius")    
print("it works")
