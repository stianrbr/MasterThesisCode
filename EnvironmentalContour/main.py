import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from misc import joint_contour_surface, contour_line_for_windspeed
matplotlib.use('tkagg')
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)

"""
Reading parameters from hdf-files created by FitDistributions project
"""

wind_params = pd.read_hdf("Parameters\\wind_param.h5", key="df")
wind_alpha = wind_params["alpha"][0]
wind_beta = wind_params["beta"][0]

hs_params = pd.read_hdf("Parameters\\hs_param.h5", key="df")
a = tuple([value for index, value in hs_params["a"].items()])
b = tuple([value for index, value in hs_params["b"].items()])

tp_params = pd.read_hdf("Parameters\\tp_param.h5", key="df")
i = tuple([value for index, value in tp_params["i"].items()])
j = tuple([value for index, value in tp_params["j"].items()])
k = tuple([value for index, value in tp_params["k"].items()])
l = tuple([value for index, value in tp_params["l"].items()])



"""
Plotting the environmental contour surface
"""
hs, tp, uw = joint_contour_surface(return_period=100, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)

fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.set_xlabel("Hs")
ax.set_ylabel("Tp")
ax.set_zlabel("Uw")
ax.plot_wireframe(hs, tp, uw, color="g", alpha=0.3)
ax.set_box_aspect((1,1,1))
plt.show()

"""
Plotting a environmental contour line
"""

windspeed = 24

hs100, tp100 = contour_line_for_windspeed(windspeed=windspeed, return_period=100, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)
hs50, tp50 = contour_line_for_windspeed(windspeed=windspeed, return_period=50, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)
hs500, tp500 = contour_line_for_windspeed(windspeed=windspeed, return_period=500, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)
hs10000, tp10000 = contour_line_for_windspeed(windspeed=windspeed, return_period=10000, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)
hs1, tp1 = contour_line_for_windspeed(windspeed=windspeed, return_period=1, seastate_dur=1, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                   a=a, b=b, i=i, j=j, k=k, l=l)
fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("Hs")
ax.set_ylabel("Tp")
plt.plot(hs10000,tp10000, color="k", label=f"R={10000} y")
plt.plot(hs500,tp500, color="b", label=f"R={500} y")
plt.plot(hs100,tp100, color="r", label=f"R={100} y")
plt.plot(hs50,tp50, color="g", label=f"R={50} y")
plt.plot(hs1,tp1, color="orange", label=f"R={1} y")
plt.title(f"Environmental contour\nUw={windspeed} m/s")
plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.tight_layout()
plt.show()
debug = True