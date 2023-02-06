lon=3.666664  # Desired longitude
lat = 60.66664  # Desired latitude

dataset = "NORA10EI"  # Desired dataset [NORA10EI;NORA3;NORA3 - Windsurfer]

startyear = 1981  # First included year
endyear = 2021  # Last included year

n_parallel = 6  # Number of parallel processes
"""
Variables to extract from dataset, 
Names must match exactly the names in the dataset.

Definition of variables may be found in the links above, by entering a folder, pressing on a .nc-file,
and selecting OPENDAP
"""
variables = ['time',
             'x_wind_10m',
             'y_wind_10m',
             'wind_speed',
             'wind_direction']  # Variables to download

savefolder = "NORA3_Wind"  # Folder where data is stored