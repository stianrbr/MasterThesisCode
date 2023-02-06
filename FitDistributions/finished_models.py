import numpy as np
from curve_fitting import power_function, sigmoid, exponential
from scipy.stats import norm

def weibull_wind(u, alpha, beta):
    """
    Cumulative marginal distribution of mean wind speed at 10 meter altitude
    :param u: (ndarray, float) Mean wind speed
    :param alpha: (float) Weibull shape parameter
    :param beta: (float) Weibull scale parameter
    :return: (ndarray, float) cumulative probability P(U<u)
    """
    return 1-np.exp(-(u/beta)**alpha)


def weibull_wave(h, u, a, b):
    """
    Cumulative conditional distribution of significant wave height,
    formulated as a 2-parameter Weibull model
    :param h: (ndarray, float) Significant wave height
    :param u: (float) Mean wind speed
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :return: (ndarray, float) cumulative probability P(U<u)
    """
    # Calculate parameterized values
    alpha_h = power_function(u, *a)  # Shape parameter
    beta_h = power_function(u, *b)  # Scale parameter
    return 1-np.exp(-(h/beta_h)**alpha_h)

def lonowe_wave(h, u, c, d, e, f, g):
    """
    Cumulative conditional distribution of significant wave height,
    formulated as a hybrid log-normal-Weibull (LoNoWe) model
    :param h: (ndarray, float) Significant wave height
    :param u: (float) Mean wind speed
    :param c:
    :param d:
    :param e:
    :param f:
    :param g:
    :return: (ndarray, float) cumulative probability P(U<u)
    """
    eta = power_function(u, *c)
    sigma = sigmoid(u, *d)
    mu = power_function(u, *e)
    beta = power_function(u, *f)
    alpha = power_function(u, *g)

    h = np.asarray(h)
    h_below_eta = h[np.where(h <= eta)]
    h_above_eta = h[np.where(h > eta)]
    F_below_eta = norm.cdf((np.log(h_below_eta)-mu)/sigma)
    F_above_eta = 1-np.exp(-(h_above_eta/beta)**alpha)
    return np.concatenate((F_below_eta, F_above_eta))


def log_normal_wave_period(t, u, h, i, j, k, l):
    t_bar = power_function(h, *i)
    u_bar = power_function(h, *j)
    nu = exponential(h, *k)
    theta = power_function(h, *l)
    mu_tp = t_bar*(1+theta*((u-u_bar)/u_bar))
    sigma_ln_tp = np.sqrt(np.log(nu**2+1))
    mu_ln_tp = np.log(mu_tp/(np.sqrt(1+nu**2)))
    return norm.cdf((np.log(t) - mu_ln_tp) / sigma_ln_tp)


