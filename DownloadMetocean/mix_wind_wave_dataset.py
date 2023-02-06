import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

datafolder_wave = "C:\\Users\\stian\\OneDrive - NTNU\\5. klasse\\Masteroppgave\\MasterThesisCode\\NORA3_Wave\\Combined_data\\"

data_wave = pd.read_hdf(datafolder_wave+"combined_dataset.h5", key="df")

datafolder_wind = "C:\\Users\\stian\\OneDrive - NTNU\\5. klasse\\Masteroppgave\\MasterThesisCode\\NORA3_Wind\\Combined_data\\"

data_wind = pd.read_hdf(datafolder_wind+"combined_dataset.h5", key="df")

data_wave.drop(data_wave.head(3).index, axis=0, inplace=True)
data_wave = data_wave.reset_index()
data_wave.drop("index", axis=1, inplace=True)


data_wind.drop(data_wind.tail(3).index, axis=0, inplace=True)
data_wind = data_wind.reset_index()
data_wind.drop("index", axis=1, inplace=True)

wind_speed = data_wind["wind_speed"]
wind_dir = data_wind["wind_direction"]

data_wave = data_wave.join(wind_speed)
data_wave = data_wave.join(wind_dir)

data_wave.to_hdf("Troll_1981_2021_complete.h5", key="df")

data_wave.drop(["ff", "dd", "hs_sea", "tp_sea", "tmp_sea", "tm1_sea", "thq_sea", "hs_swell", "tp_swell", "thq_swell"], axis=1, inplace=True)
data_wave.to_hdf("Troll_1981_2021_selected.h5", key="df")
debug = True