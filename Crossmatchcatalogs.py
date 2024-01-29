import os
from multiprocessing import Pool, cpu_count
import matplotlib.colors as mcolors
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview, newprojplot
import time as time

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import astropy_healpix as ah
import ligo.skymap.plot
import re
from ligo.skymap.postprocess import find_greedy_credible_levels
from astropy.table import Table, QTable
from astropy import units as u
from astropy.io import fits
from ligo.skymap.io.fits import read_sky_map
#import pandas as pd
from astropy.coordinates import Angle, SkyCoord
import matplotlib.colors as mcolors
from tqdm.notebook import tqdm

###load the catalog you want to obtain more information
data = Table.read('/luidhy_docker/projeto/GW_dark_sirens/GW_DKsirens_O4/S230922g/S230922g.multiorder_table.fits', format='fits')
#data = Table.read('/luidhy_docker/projeto/GW_dark_sirens/GW_DKsirens_O4/S230922g/Dominguez_RADEC.fits', format='fits')
df_my = data.to_pandas()
df_my.info()

###select just the collumns needed for this file
#objid = np.array(data['objID'])
ra_my = np.array(data['RA'])
dec_my = np.array(data['DEC'])
#TType = np.array(data['TType'])

###load the general catalog you want, for instande DELVE, LEGACY
def process_file(filename):
    var = Table.read(filename, format='fits')
    var.keep_columns(['RA', 'DEC'])    ###collumns you want to keep 'Z_PHOT_PEAK'
    return var


if __name__ == "__main__":
    # Loading the fits files
    #directory = '/luidhy_docker/astrodados/DELVE_DR2_cleaned' ###path yo your general ancilllary data
    directory = '/luidhy_docker/astrodados/LEGACY_DR10_phz'
    filenames = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".fits")]

    # Specify the number of processes (adjust as needed)
    num_processes = min(cpu_count(), len(filenames))  # Use the minimum of CPU cores and files
    print(f"Using {num_processes} processes.")

    total_files = len(filenames)
    processed_files = 0

    # Using multiprocessing to parallelize file processing
    with Pool(processes=num_processes) as pool:
        results = []
        for result in pool.imap_unordered(process_file, filenames):
            results.append(result)
            processed_files += 1
            progress = processed_files / total_files * 100
            print(f"Progress: {progress:.2f}%", end='\r')

    # Concatenating the processed tables
    stack = np.hstack(results)

##collumns from the stacked data
ra_add = stack['RA']
dec_add = stack['DEC']


print(len(ra_add))

#####doing the crossmatch
print("Starting the crossmatch")
inicio = time.time()
print("running the match")
c_my = SkyCoord(ra=ra_my*u.degree, dec=dec_my*u.degree)
c_add = SkyCoord(ra=ra_add*u.degree, dec=dec_add*u.degree)
idx, d2d, d3d = c_my.match_to_catalog_sky(c_add)
#max_separation = 1.0*u.arcsec
#sep_constraint = d2d < max_separation
final = time.time()
print("time to peform the match:",final - inicio, "seconds")
print("time to perform the match:", (final - inicio)/60, "minutes")
print("Size of the match cataog:",len(idx))


max_separation = 1.0 * u.arcsec
filtered_indices = [i for i, d in enumerate(d2d) if d < max_separation]

print(len(filtered_indices))

# Create new arrays for the objects variables that match the separation condition
filtered_ra_my = ra_my[filtered_indices]
filtered_dec_my = dec_my[filtered_indices]
filtered_ra_add = ra_add[idx[filtered_indices]]
filtered_dec_add = dec_add[idx[filtered_indices]]

filtered_table = Table([filtered_ra_my, filtered_dec_my, filtered_ra_add, filtered_dec_add],
                       names=['RA_my', 'DEC_my', 'RA_add', 'DEC_add'])

# Saving the table to a FITS file
filtered_table.write('filtered_catalog_teste3.fits', overwrite=True)
# Print the number of objects in the filtered catalog
print("Number of objects with separation < 1 arcsec:", len(filtered_table))










