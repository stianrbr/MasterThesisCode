import numpy as np
import pandas as pd
from read_and_correct_data import bin_data, tp_correction
from distribution_fitting import fit_weibull_wind, fit_weibull_waves, fit_lognormal_period


#Path to the (uncorrected) input dataset
data_path = ""
data_name = "Troll_1981_2021_selected.h5"
data = pd.read_hdf(data_path+data_name, key="df")

# Correct Tp and save, or read a previously correct dataset
correct_tp = False
if correct_tp:
    data["tp"] = tp_correction(data["tp"], data["hs"], plot=False, show_plot=False, save_plot=False)
    data.to_hdf(data_path + data_name.removesuffix(".h5") + "_corr_tp.h5", key="df")
else:
    data_name = "Troll_1981_2021_selected_corr_tp.h5"
    data = pd.read_hdf(data_path + data_name, key="df")

# Bin width for grouping data
wind_bin_width = 0.5
wave_bin_width = 0.5

# Generate bins for wind and wave for conditional distribution
wind_bin_height, wind_bin_edges, wind_mid_bin = bin_data(data["wind_speed"], wind_bin_width)
wave_bin_height, wave_bin_edges, wave_mid_bin = bin_data(data["hs"], wave_bin_width)




alpha_wind, beta_wind = fit_weibull_wind(data["wind_speed"], plot=False, show=False, save=False)
wind_param = pd.DataFrame(np.array([[alpha_wind, beta_wind]]),
                          columns=["alpha", "beta"])
wind_param.to_hdf("Parameters\\"+"wind_param.h5", key="df")


a, b = fit_weibull_waves(data, wind_bin_edges, wind_bin_height, plot_separate=False, plot_fitting=False, show=False, save=False)
wave_height_param = pd.DataFrame(np.stack((a,b)).T, columns=["a", "b"])
wave_height_param.to_hdf("Parameters\\"+"hs_param.h5", key="df")

i, j, k, l = fit_lognormal_period(data, wind_mid_bin, wind_bin_edges, wave_mid_bin, wave_bin_edges, plot=True, show=True, save=True)
wave_period_param = pd.DataFrame(np.stack((i,j, k, l)).T, columns=["i", "j", "k", "l"])
wave_period_param.to_hdf("Parameters\\"+"tp_param.h5", key="df")

