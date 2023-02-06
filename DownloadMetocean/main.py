"""
Author: Stian Brurås
Date: 26.01.2023

"""

# Import statements
import concurrent.futures
import os
from itertools import repeat
from download_functions import download_year_NORA3, download_year_NORA3_wave, download_year_NORA10
from misc import *
from inputs import *

"""
Code for downloading environmental data from MET Norway
Supported datasets:
    - NORA10EI: https://thredds.met.no/thredds/catalog/nora10ei/catalog.html
        -ERA-Interim
        -10 km spatial and 1 hr temporal resolution
        -Wind and waves
    - NORA3: https://thredds.met.no/thredds/catalog/nora3/catalog.html
        -ERA-5
        -3 km spatial and 1 hr temporal resolution
        -Wind
    - NORA3 - Windsurfer: https://thredds.met.no/thredds/catalog/windsurfer/mywavewam3km_files/catalog.html
        -Wave model based on NORA3 dataset
        -3 km spatial and 1 hr temporal resolution
"""

years = np.arange(startyear, endyear+1)


# Create folder if it does not exist
try:
    os.mkdir(savefolder)
except FileExistsError:
    pass

if dataset=="NORA3":
    func = download_year_NORA3
elif dataset=="NORA3 - Windsurfer":
    func = download_year_NORA3_wave
elif dataset=="NORA10EI":
    func = download_year_NORA10
else:
    raise Exception("Choose an existing data set")


if __name__ == '__main__':
    """
    Download data using multiprocessing.
    Be careful with using many processes, might get banned from website if pulling to many requests.
    """

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_parallel) as executor:
        executor.map(func, years, repeat(variables), repeat(lat), repeat(lon), repeat(savefolder))