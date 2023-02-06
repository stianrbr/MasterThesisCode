import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
from scipy.optimize import curve_fit
import nmmn.plots
import astropy
import seaborn as sns
sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)

from finished_models import weibull_wind
from curve_fitting import power_function, linear_function, exponential


def fit_weibull_wind(wind_data, plot=False, show=False, save=False):
    """
    Fit a 2-parameter Weibull distribution to wind data using MLE.
    Plotting functionality included.
    :param wind_data: (ndarray) Wind data
    :param plot: (bool) Flag to decide if plotting
    :param show: (bool) Flag to decide showing of plot
    :param save: (bool) Flag to decide saving of plot
    :return: shape: (float) Weibull shape parameter
             scale: (float) Weibull scale parameter
    """
    shape, loc, scale = stats.weibull_min.fit(wind_data, floc=0, method="MLE")  # Scipy fitting

    if plot:
        f_k = np.arange(1, len(wind_data) + 1) / (len(wind_data) + 1)  # Empirical CDF
        F = weibull_wind(np.sort(wind_data), shape, scale)  # Fitted CDF
        fig = plt.figure()
        plt.scatter(np.log(np.sort(wind_data)), np.log(-np.log(1 - f_k)), c="k", s=2, label="Sample")
        plt.plot(np.log(np.sort(wind_data)), np.log(-np.log(1 - F)), c="r", label="Fitted\ndistribution")
        # Label and title uses latex formatting.
        plt.xlabel(r'$\ln(U_{w})$')
        plt.ylabel(r'$\ln(-\ln(1-F_{U_w}))$')
        plt.title(
            r'Marginal distribution of mean wind speed ($U_{w}$)' + "\n" + '2-parameter Weibull model - Weibull scale',
            fontsize=14)
        lgnd = plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        lgnd.legendHandles[0].set_sizes([20.0])  # Setting point size in legend.
        plt.tight_layout()
        if save:
            plt.savefig("Weibull_plot_wind.png")
        if show:
            plt.show()
        else:
            plt.close(fig)
    return shape, scale

def fit_weibull_waves(data, wind_bin_edges, wind_bin_height, plot_separate=False, plot_fitting=False, show=False, save=False):
    """
    Fit a 2-parameter Weibull distribution to sign. wave height inside each wind bin.
    Approximate Weibull parameters as a function of wind speed by curve fitting.
    :param data: (DataFrame) Hindcast data
    :param wind_bin_edges: (ndarray) Edges of wind bins, including leftmost and rightmost edge
    :param wind_bin_height: (ndarray) Number of values inside each wind bin
    :param plot_separate: (bool) Flag for plotting distribution inside each bin, default=False
    :param plot_fitting: (bool) Flag for plotting curve fit approximations, default=False
    :param show: (bool) Flag for showing plots, default=False
    :param save: (bool) Flag for saving plots, default=False
    :return: a_param: (ndarray) Curve fit results (power function) for Weibull shape parameter
             b_param: (ndarray) Curve fit results (power function) for Weibull scale parameter
    """

    try:  # Creating directory for results
        os.mkdir("H_s_Weibull_fit\\")
    except FileExistsError:
        pass

    df = pd.DataFrame(columns=["Bin left", "Bin mid", "Bin right", "Shape", "Scale"])  # Allocating empty dataframe
    """
    Loop over all wind bins and fit separate Weibull distributions
    """
    for i in range(len(wind_bin_edges)-1):
        if wind_bin_height[i]<50:  # Only fit distribution if sample size larger than 50
            break
        left_edge, right_edge = wind_bin_edges[i], wind_bin_edges[i+1]
        filtered_data = data[(data["wind_speed"]>=left_edge) & (data["wind_speed"]<right_edge)]
        hs_bin = filtered_data["hs"].to_numpy()
        f_k_hs_bin = np.arange(1, len(hs_bin) + 1) / (len(hs_bin) + 1)
        shape, loc, scale = stats.weibull_min.fit(hs_bin, floc=0, method="MLE")  # Scipy fitting
        temp = {"Bin left": [left_edge], "Bin mid": [(right_edge - left_edge) / 2 + left_edge],
                "Bin right": [right_edge], "Shape": [shape], "Scale": [scale]}
        df = pd.concat([df, pd.DataFrame(temp)])  # Append results
        if plot_separate:
            F_bin = stats.weibull_min.cdf(np.sort(hs_bin), c=shape, loc=loc, scale=scale)
            fig = plt.figure()
            plt.scatter(np.log(np.sort(hs_bin)), np.log(-np.log(1 - f_k_hs_bin)), c="k", s=2, label="Raw data")
            plt.plot(np.log(np.sort(hs_bin)), np.log(-np.log(1 - F_bin)), c="r", label="Fitted\ndistribution")
            plt.xlabel(r'$\ln(U_{w})$')
            plt.ylabel(r'$\ln(-\ln(1-F_{U_w}))$')
            plt.title(r'Distribution of significant wave height given' +"\n"+'mean wind speed ($H_{s}|U_{W}\in $'+f"[{left_edge} m/s ,{right_edge} m/s))"+ "\n" + 'Weibull model - Weibull scale', fontsize=14)
            plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
            plt.tight_layout()
            if save:
                plt.savefig(f"H_s_Weibull_fit\\Hs_U_{left_edge}_{right_edge}.png")
            if show:
                plt.show()
            else:
                plt.close(fig)
        if plot_separate:
            F_bin = stats.weibull_min.cdf(np.sort(hs_bin), c=shape, loc=loc, scale=scale)
            fig = plt.figure()
            plt.scatter(np.sort(hs_bin), f_k_hs_bin, c="k", s=2, label="Raw data")
            plt.plot(np.sort(hs_bin), F_bin, c="r", label="Fitted\ndistribution")
            plt.xlabel(r'$U_{w}$')
            plt.ylabel(r'$F_{U_w}$')
            plt.title(r'Distribution of significant wave height given' +"\n"+'mean wind speed ($H_{s}|U_{W}\in $'+f"[{left_edge} m/s ,{right_edge} m/s))"+ "\n" + 'Weibull model', fontsize=14)
            plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
            plt.tight_layout()
            if save:
                plt.savefig(f"H_s_Weibull_fit\\Hs_U_{left_edge}_{right_edge}_dist.png")
            if show:
                plt.show()
            else:
                plt.close(fig)
    df.to_hdf("H_s_Weibull_fit\\fitted_parameters.h5", key="df")

    """
    Approximate parameters by curve fitting
    """

    mid_bin = df["Bin mid"].to_numpy()
    shapelist = df["Shape"].to_numpy()
    scalelist = df["Scale"].to_numpy()
    # Initial guess must be changed manually
    a_param, _ = curve_fit(power_function, mid_bin, shapelist, (1.75, 0.07, 1.24))
    b_param, _ = curve_fit(power_function, mid_bin, scalelist, (1.34, 0.01, 1.97))
    if plot_fitting:
        fig = plt.figure()
        plt.scatter(mid_bin, shapelist)
        plt.plot(mid_bin, power_function(mid_bin, *a_param),
                 label=r"$\alpha_{HC}(u)$ (Shape)")
        plt.scatter(mid_bin, scalelist)
        plt.plot(mid_bin, power_function(mid_bin, *b_param),
                 label=r"$\beta_{HC}(u)$ (Scale)")
        plt.xlabel(r"$U_w$ [m/s]")
        plt.ylabel("Value [-]")
        plt.title(r"Curve fitting of conditional $H_s$ distribution parameters")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.tight_layout()
        if save:
            plt.savefig("H_s_Weibull_fit\\shape_scale.png")
        if show:
            plt.show()
        else:
            plt.close(fig)
    return a_param, b_param

def fit_lognormal_period(data, wind_mid_bin, wind_bin_edges, wave_mid_bin, wave_bin_edges,
                         plot=False, show=False, save=False):
    """
    :param data: (DataFrame) Hindcast data
    :param wind_mid_bin: (ndarray) Midpoint of wind bins
    :param wind_bin_edges: (ndarray) Edges of wind bins, including leftmost and rightmost edge
    :param wave_mid_bin: (ndarray) Midpoint of wave bins
    :param wave_bin_edges: (ndarray) Edges of wave bins, including leftmost and rightmost edge
    :param plot: (bool) Flag for plotting, default=False
    :param show: (bool) Flag for showing plots, default=False
    :param save: (bool) Flag for saving plots, default=False
    :return: i_param: (ndarray) Curve fit results (power function) for expected spectral period
             j_param: (ndarray) Curve fit results (power function) for expected mean wind speed
             k_param: (ndarray) Curve fit results (exponential) for Tp coefficient of variance
             l_param: (ndarray) Curve fit results (power function) for theta (slope)
    """
    try:
        os.mkdir("T_p_LogNormal_fit\\")
    except FileExistsError:
        pass
    res_mu = []
    res_sigma = []
    res_nu = []
    for i in range(len(wind_bin_edges) - 1):
        wave_res_mu = np.zeros_like(wave_mid_bin)
        wave_res_sigma = np.zeros_like(wave_mid_bin)
        wave_res_nu = np.zeros_like(wave_mid_bin)
        for j in range(len(wave_bin_edges) - 1):
            left_wind_edge = wind_bin_edges[i]
            right_wind_edge = wind_bin_edges[i + 1]
            left_wave_edge = wave_bin_edges[j]
            right_wave_edge = wave_bin_edges[j + 1]
            filtered_wind = data.loc[(data["wind_speed"] >= left_wind_edge) & (data["wind_speed"] < right_wind_edge)]
            filtered = filtered_wind.loc[
                (filtered_wind["hs"] >= left_wave_edge) & (filtered_wind["hs"] < right_wave_edge)]
            tp_bin = filtered["tp"].to_numpy()

            if len(tp_bin) < 3:
                mu_tp = np.nan
                sigma_tp = np.nan
                nu_tp = np.nan
            else:
                mu_tp = np.mean(tp_bin)
                sigma_tp = np.std(tp_bin, ddof=1)
                nu_tp = sigma_tp / mu_tp
            wave_res_mu[j] = mu_tp
            wave_res_sigma[j] = sigma_tp
            wave_res_nu[j] = nu_tp
        res_mu.append(wave_res_mu)
        res_sigma.append(wave_res_sigma)
        res_nu.append(wave_res_nu)
    res_mu = np.array(res_mu)
    res_sigma = np.array(res_sigma)
    res_nu = np.array(res_nu)
    if plot:
        ##
        my_cmap = nmmn.plots.turbocmap()

        fig, ax = plt.subplots()
        im = ax.pcolormesh(wave_bin_edges, wind_bin_edges, res_mu, edgecolor="whitesmoke", linewidth=0.005,
                           cmap=my_cmap)
        fig.colorbar(im, ax=ax)
        plt.xlabel(r"$H_{s}$")
        plt.ylabel(r"$U_{w}$")
        plt.title(r"Mean of peak spectral period ($\mu_{T_p}$)")
        plt.tight_layout()
        if save:
            plt.savefig("T_p_LogNormal_fit\\mu_tp_colormesh.png")
        if show:
            plt.show()
        else:
            plt.close(fig)

        fig, ax = plt.subplots()
        im = plt.pcolormesh(wave_bin_edges, wind_bin_edges, res_sigma, edgecolors="whitesmoke", linewidth=.003,
                            cmap=my_cmap)
        fig.colorbar(im, ax=ax)
        plt.xlabel(r"$H_{s}$")
        plt.ylabel(r"$U_{w}$")
        plt.title(r"Standard deviation of peak spectral period ($\sigma_{T_p}$)")
        plt.tight_layout()
        if save:
            plt.savefig("T_p_LogNormal_fit\\sigma_tp_colormesh.png")
        if show:
            plt.show()
        else:
            plt.close(fig)

        fig, ax = plt.subplots()
        im = plt.pcolormesh(wave_bin_edges, wind_bin_edges, res_nu, edgecolors="whitesmoke", linewidth=.003,
                            cmap=my_cmap)
        fig.colorbar(im, ax=ax)
        plt.xlabel(r"$H_{s}$")
        plt.ylabel(r"$U_{w}$")
        plt.title(r"Coefficient of variance of peak spectral period ($\nu_{T_p}$)")
        plt.tight_layout()
        if save:
            plt.savefig("T_p_LogNormal_fit\\nu_tp_colormesh.png")
        if show:
            plt.show()
        else:
            plt.close(fig)
    t_bar_h = np.zeros_like(wave_mid_bin)
    u_bar_h = np.zeros_like(wave_mid_bin)
    for j in range(len(wave_bin_edges) - 1):
        left_wave_edge = wave_bin_edges[j]
        right_wave_edge = wave_bin_edges[j + 1]
        filtered_wave = data.loc[(data["hs"] >= left_wave_edge) & (data["hs"] < right_wave_edge)]
        t_bar_h[j] = np.mean(filtered_wave["tp"].to_numpy())
        u_bar_h[j] = np.mean(filtered_wave["wind_speed"].to_numpy())
    sigma = np.ones(len(wave_mid_bin) + 1)
    sigma[[0]] = 0.01
    i_param, _ = curve_fit(power_function, wave_mid_bin, t_bar_h, (8, 0.6, 0.9))
    j_param, _ = curve_fit(power_function, np.concatenate((np.array([0]), wave_mid_bin)),
                           np.concatenate((np.array([0]), u_bar_h)), (4, 2, 0.9), sigma=sigma)
    if plot:
        fig = plt.figure()
        plt.scatter(wave_mid_bin, t_bar_h, label=r"$\bar{T}_p(h)$")
        plt.plot(wave_mid_bin, power_function(wave_mid_bin, *i_param))
        plt.scatter(wave_mid_bin, u_bar_h, label=r"$\bar{u}(h)$")
        plt.plot(wave_mid_bin, power_function(wave_mid_bin, *j_param))
        plt.title(r"Approximation of expected spectral peak period $\bar{T}_p(h)$"+" and \n"+ r"expected mean wind speed $\bar{u}(h)$ as a function of $H_{s}$")
        plt.xlabel(r"$H_{s} [m]$")
        plt.ylabel(r"$\bar{T}_p(h) [s]$, $\bar{U}(h) [m/s]$")
        lgnd = plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        for handle in lgnd.legendHandles:
            handle.set_sizes([20.0])
        plt.tight_layout()
        if save:
            plt.savefig("T_p_LogNormal_fit\\tbar_ubar_fit.png")
        if show:
            plt.show()
        else:
            plt.close(fig)

    inter_res, theta_res, hs_res = [], [], []
    for i in range(len(wave_mid_bin)):
        hs = wave_mid_bin[i]
        mu = res_mu[:, i]
        mw = wind_mid_bin
        mu_dropnan = mu[~np.isnan(mu)]
        mw_dropnan = mw[~np.isnan(mu)]
        if len(mu_dropnan) > 5:
            hs_res.append(hs)
            t_bar = power_function(hs, *i_param)
            u_bar = power_function(hs, *j_param)
            norm_tp = (mu_dropnan - t_bar) / t_bar
            norm_u = (mw_dropnan - u_bar) / u_bar
            param, _ = curve_fit(linear_function, norm_u, norm_tp, (0.05, -0.06))
            theta_res.append(param[1])
            inter_res.append(param[0])
            plot = 1
            show = 0
            if plot:
                fig = plt.figure()
                plt.scatter(norm_u, norm_tp, label="Raw data")
                plt.plot(norm_u, linear_function(norm_u, *param), label="Curve fit")
                plt.xlabel(r"$\frac{u-\bar{u}(h)}{\bar{u}(h)}$")
                plt.ylabel(r"$\frac{\bar{T}_p(u,h)-\bar{T}_p(h)}{\bar{T}_p(h)}$")
                plt.title(
                    r"Approximation of $\theta$ for $H_s=$" + f"{hs} m\n Linear approximation " + r"($\gamma = 1)$")
                plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                plt.tight_layout()
                if save:
                    plt.savefig(f"T_p_LogNormal_fit\\theta_for_hs_{hs}.png")
                if show:
                    plt.show()
                else:
                    plt.close(fig)
    theta_res = np.concatenate((np.array([0]), np.asarray(theta_res)))
    hs_res = np.concatenate((np.array([0]), np.asarray(hs_res)))
    theta_mean = np.mean(theta_res)
    sigma = np.ones(len(hs_res))
    sigma[[0]] = 0.01
    l_param, _ = curve_fit(power_function, hs_res, theta_res, (0.2, -0.32, 0.4), sigma=sigma)

    if plot:
        fig = plt.figure()
        plt.scatter(hs_res, theta_res, label="Estimated slopes")
        plt.plot(hs_res, power_function(hs_res, *l_param), label="Power function")
        plt.plot(hs_res, np.repeat(theta_mean, len(hs_res)), label="Mean")
        plt.title(r"Approximation of $\theta$")
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        plt.tight_layout()
        if save:
            plt.savefig(f"T_p_LogNormal_fit\\theta_fitting.png")
        if show:
            plt.show()
        else:
            plt.close(fig)
    nu_hs_bins = []
    corr_hs = []
    for i in range(len(wave_bin_edges) - 1):
        left_wave_edge = wave_bin_edges[i]
        right_wave_edge = wave_bin_edges[i + 1]
        filtered_wave = data.loc[(data["hs"] >= left_wave_edge) & (data["hs"] < right_wave_edge)]
        tp_bin = filtered_wave["tp"]
        if len(tp_bin) > 2:
            mu_tp = np.mean(tp_bin)
            sigma_tp = np.std(tp_bin, ddof=1)
            nu_tp = sigma_tp / mu_tp
            nu_hs_bins.append(nu_tp)
            corr_hs.append(wave_mid_bin[i])
    nu_hs_bins = np.array(nu_hs_bins)
    corr_hs = np.array(corr_hs)
    k_param, _ = curve_fit(exponential, corr_hs, nu_hs_bins, (-0.05, 0.3, -0.2))
    if plot:
        fig = plt.figure()
        plt.scatter(corr_hs, nu_hs_bins, )
        plt.plot(corr_hs, exponential(corr_hs, *k_param))
        plt.xlabel(r"$H_{s}$")
        plt.ylabel(r"$\nu_{T_p}$")
        plt.title(r"Approximation of $\nu_{T_p}(h)$")
        if save:
            plt.savefig(f"T_p_LogNormal_fit\\nu_fitting.png")
        if show:
            plt.show()
        else:
            plt.close(fig)
    return i_param, j_param, k_param, l_param