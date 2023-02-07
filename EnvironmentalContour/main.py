import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
from misc import joint_contour_surface, contour_line_for_windspeed, u_space_surface
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
Creating folder for storing results
"""

savefolder = "Contour_results\\"

try:
    os.mkdir(savefolder)
except FileExistsError:
    pass


"""
Plotting the environmental contour surface
"""

plot_surf = True
plot_surf_uspace = True
show_surf = True
save_surf = True
if plot_surf:
    return_period = 100
    seastate_dur = 1
    if plot_surf_uspace:
        x, y, z = u_space_surface(return_period, seastate_dur)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel("u2")
        ax.set_ylabel("u3")
        ax.set_zlabel("u1")
        ax.plot_wireframe(x, y, z, color="g", alpha=0.3)
        ax.set_box_aspect((1, 1, 1))
        if save_surf:
            plt.savefig(savefolder + f"u_surf_{return_period}y.png")
        if show_surf:
            plt.show()
        else:
            plt.close(fig)
    hs, tp, uw = joint_contour_surface(return_period=return_period, seastate_dur=seastate_dur, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                       a=a, b=b, i=i, j=j, k=k, l=l)
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    ax.set_xlabel(r"$H_s$ [m]")
    ax.set_ylabel(r"$T_p$ [s]")
    ax.set_zlabel(r"$U_w$ [m/s]")
    ax.plot_wireframe(hs, tp, uw, color="g", alpha=0.3)
    ax.set_box_aspect((1,1,1))
    if save_surf:
        plt.savefig(savefolder+f"contour_surf_{return_period}y.png")
    if show_surf:
        plt.show()
    else:
        plt.close(fig)

"""
Plotting a environmental contour line
"""
plot_single_line = True
save_single_line = True
show_single_line = True
if plot_single_line:
    windspeed = 24
    return_period = 100
    seastate_dur = 1

    hs, tp = contour_line_for_windspeed(windspeed=windspeed, return_period=return_period, seastate_dur=seastate_dur, alpha_wind=wind_alpha, beta_wind=wind_beta,
                                       a=a, b=b, i=i, j=j, k=k, l=l)
    fig = plt.figure()
    ax = plt.axes()
    ax.set_xlabel(r"$H_s$ [m]")
    ax.set_ylabel(r"$T_p$ [s]")
    plt.plot(hs,tp, color="orange", label=f"R={return_period} y")
    plt.title(f"Environmental contour\nUw={windspeed} m/s")
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
    plt.tight_layout()
    if save_single_line:
        plt.savefig(savefolder + f"contour_line_ws_{windspeed}_r_{return_period}y.png")
    if show_single_line:
        plt.show()
    else:
        plt.close(fig)

"""
Plotting multiple return periods and wind speeds
"""
plot_multiple_line = True
save_multiple_line = True
show_multiple_line = False

windspeeds = np.arange(10, 35, step=2)
return_periods = [10000, 100, 50, 5, 1]
seastate_dur = 1
if plot_multiple_line:
    for wind in windspeeds:
        fig = plt.figure()
        ax = plt.axes()
        ax.set_xlabel(r"$H_s$ [m]")
        ax.set_ylabel(r"$T_p$ [s]")
        clrs = sns.color_palette("muted", len(return_periods))
        for c, r in enumerate(return_periods):
            try:
                hs, tp = contour_line_for_windspeed(windspeed=wind, return_period=r, seastate_dur=seastate_dur,
                                                alpha_wind=wind_alpha, beta_wind=wind_beta, a=a, b=b, i=i, j=j, k=k, l=l)
            except:
                print(f"Failed calculating return period: {r} for wind speed {wind}")
                pass
            else:
                plt.plot(hs, tp, color=clrs[c], label=f"R={r} y")
        plt.title(f"Environmental contour\nUw={wind} m/s")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.tight_layout()
        if save_multiple_line:
            plt.savefig(savefolder + f"contour_line_w_{wind}_{'_'.join(str(r) for r in return_periods)}y.png")
        if show_multiple_line:
            plt.show()
        else:
            plt.close(fig)
