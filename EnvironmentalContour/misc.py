import numpy as np
from scipy.stats import norm
from curve_fitting import power_function, exponential
def weibull_wind(u, alpha, beta):
    """
    Cumulative marginal distribution of mean wind speed at 10 meter altitude
    :param u: (ndarray, float) Mean wind speed
    :param alpha: (float) Weibull shape parameter
    :param beta: (float) Weibull scale parameter
    :return: (ndarray, float) cumulative probability P(U<u)
    """
    return 1-np.exp(-(u/beta)**alpha)

def weibull_wind_pdf(u, alpha, beta):
    """
    Probability density distribution of mean wind speed at 10 meter altitude
    :param u:
    :param alpha:
    :param beta:
    :return:
    """
    return alpha/beta * (u/beta)**(alpha-1)*np.exp(-(u/beta)**alpha)

def weibull_wave(h, u, a, b):
    """
    Cumulative conditional distribution of significant wave height,
    formulated as a 2-parameter Weibull model
    :param h: (ndarray, float) Significant wave height
    :param u: (float) Mean wind speed
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :return: (ndarray, float) cumulative probability P(H<h|u)
    """
    # Calculate parameterized values
    alpha_h = power_function(u, *a)  # Shape parameter
    beta_h = power_function(u, *b)  # Scale parameter
    return 1-np.exp(-(h/beta_h)**alpha_h)

def weibull_wave_pdf(h, u, a, b):
    """
    Cumulative conditional distribution of significant wave height,
    formulated as a 2-parameter Weibull model
    :param h: (ndarray, float) Significant wave height
    :param u: (float) Mean wind speed
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :return: (ndarray, float) cumulative probability P(H<h|u)
    """
    # Calculate parameterized values
    alpha_h = power_function(u, *a)  # Shape parameter
    beta_h = power_function(u, *b)  # Scale parameter
    return alpha_h/beta_h * (h/beta_h)**(alpha_h-1)*np.exp(-(h/beta_h)**alpha_h)

def joint_hs_tp(u, h, t, a, b, i, j, k, l):
    temp = np.nan_to_num(weibull_wave_pdf(h, u, a, b)*log_normal_wave_period_pdf(t, u, h, i, j, k, l))
    return np.trapz(temp, u)

def log_normal_wave_period(t, u, h, i, j, k, l):
    """
    Cumulative conditional distribution of spectral peak period,
    formulated as a log-normal model
    :param t: (ndarray, float) Spectral peak period
    :param u: (float) Mean wind speed
    :param h: (float) Significant wave height
    :param i: (tuple) Coeff. for power func. parameterization of expected spectral peak period
    :param j: (tuple) Coeff. for power func. parameterization of expected mean wind speed
    :param k: (tuple) Coeff. for exponential func. parameterization of coeff. of variance
    :param l: (tuple) Coeff. for power func. parameterization of theta (slope)
    :return: (ndarray, float) cumulative probability P(T<t|u,h)
    """
    t_bar = power_function(h, *i)
    u_bar = power_function(h, *j)
    nu = exponential(h, *k)
    theta = power_function(h, *l)
    mu_tp = t_bar*(1+theta*((u-u_bar)/u_bar))
    sigma_ln_tp = np.sqrt(np.log(nu**2+1))
    mu_ln_tp = np.log(mu_tp/(np.sqrt(1+nu**2)))
    return norm.cdf((np.log(t) - mu_ln_tp) / sigma_ln_tp)

def log_normal_wave_period_pdf(t, u, h, i, j, k, l):
    """
    Cumulative conditional distribution of spectral peak period,
    formulated as a log-normal model
    :param t: (ndarray, float) Spectral peak period
    :param u: (float) Mean wind speed
    :param h: (float) Significant wave height
    :param i: (tuple) Coeff. for power func. parameterization of expected spectral peak period
    :param j: (tuple) Coeff. for power func. parameterization of expected mean wind speed
    :param k: (tuple) Coeff. for exponential func. parameterization of coeff. of variance
    :param l: (tuple) Coeff. for power func. parameterization of theta (slope)
    :return: (ndarray, float) cumulative probability P(T<t|u,h)
    """
    t_bar = power_function(h, *i)
    u_bar = power_function(h, *j)
    nu = exponential(h, *k)
    theta = power_function(h, *l)
    mu_tp = t_bar*(1+theta*((u-u_bar)/u_bar))
    sigma_ln_tp = np.sqrt(np.log(nu**2+1))
    mu_ln_tp = np.log(mu_tp/(np.sqrt(1+nu**2)))
    return 1/(t*sigma_ln_tp*np.sqrt(2*np.pi)) * np.exp(-((np.log(t)-mu_ln_tp)**2/(2*sigma_ln_tp**2)))

def transform_u1(u1, alpha, beta):
    """
    Transformation of first variable in u-space to mean wind speed in physical space
    :param u1: (ndarray, float) u-space variable, normal distributed variable
    :param alpha: (float) Weibull shape parameter
    :param beta: (float) Weibull scale parameter
    :return: (ndarray, float) Mean wind speed
    """
    return beta*(-np.log(1-norm.cdf(u1)))**(1/alpha)

def tranform_u2(u2, uw, a, b):
    """
    Transformation of second variable in u-space to sign. wave height in physical space
    :param u2: (ndarray, float) u-space variable, normal distributed variable
    :param uw: (float) mean wind speed, for conditional distribution
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :return: (ndarray, float) Sign. wave height
    """
    alpha_h = power_function(uw, *a)  # Shape parameter
    beta_h = power_function(uw, *b)  # Scale parameter
    return beta_h*(-np.log(1-norm.cdf(u2)))**(1/alpha_h)

def transform_u3(u3, uw, h, i, j, k, l):
    """
    Transformation of third variable in u-space to spectral peak period in physical space
    :param u3: (ndarray, float) u-space variable, normal distributed variable
    :param uw: (float) mean wind speed, for conditional distribution
    :param h: (float) sign. wave height, for conditional distribution
    :param i: (tuple) Coeff. for power func. parameterization of expected spectral peak period
    :param j: (tuple) Coeff. for power func. parameterization of expected mean wind speed
    :param k: (tuple) Coeff. for exponential func. parameterization of coeff. of variance
    :param l: (tuple) Coeff. for power func. parameterization of theta (slope)
    :return: (ndarray, float) Spectral peak period
    """
    t_bar = power_function(h, *i)
    u_bar = power_function(h, *j)
    nu = exponential(h, *k)
    theta = power_function(h, *l)
    mu_tp = t_bar*(1+theta*((uw-u_bar)/u_bar))
    sigma_ln_tp = np.sqrt(np.log(nu**2+1))
    mu_ln_tp = np.log(mu_tp/(np.sqrt(1+nu**2)))
    return np.exp(sigma_ln_tp*u3+mu_ln_tp)

def u_from_uw(uw, alpha, beta):
    """
    Transform mean wind speed from physical space to u-space
    :param uw: (ndarray, float) Mean wind speed
    :param alpha: (float) Weibull shape parameter
    :param beta: (float) Weibull scale parameter
    :return: (ndarray, float) Corresponding value in u-space
    """
    F = weibull_wind(uw, alpha, beta)
    return norm.ppf(F)

def joint_contour_surface(return_period, seastate_dur, alpha_wind, beta_wind, a, b, i, j, k, l):
    """
    Create a sphere in u-space for a given return period and transform it to physical space
    :param return_period: (float) Return period for exceedance probability
    :param seastate_dur: (float) Duration of considered sea states
    :param alpha_wind: (float) Weibull shape parameter
    :param beta_wind: (float) Weibull sscale parameter
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :param i: (tuple) Coeff. for power func. parameterization of expected spectral peak period
    :param j: (tuple) Coeff. for power func. parameterization of expected mean wind speed
    :param k: (tuple) Coeff. for exponential func. parameterization of coeff. of variance
    :param l: (tuple) Coeff. for power func. parameterization of theta (slope)
    :return: hs: (ndarray) Sign. wave height points of sphere in physical space
             tp: (ndarray) Spectral peak period points of sphere in physical space
             uw: (ndarray) Mean wind speed points of sphere
    """
    q = 1 / return_period
    n = 24 / seastate_dur * 365

    beta = -1 * norm.ppf(q / n)
    n_grid, m_grid = np.mgrid[0:2 * np.pi:500j, 0:np.pi / 2:500j]
    u1 = beta * np.cos(m_grid)
    u2 = beta*np.cos(n_grid)*np.sin(m_grid)
    u3 = beta*np.sin(n_grid)*np.sin(m_grid)
    uw = transform_u1(u1, alpha_wind, beta_wind)
    hs = tranform_u2(u2, uw, a, b)
    tp = transform_u3(u3, uw, hs, i, j, k, l)
    return hs, tp, uw

def u_space_surface(return_period, seastate_dur):
    """
    Create a sphere in u-space for a given return period
    :param return_period: (float) Return period for exceedance probability
    :param seastate_dur: (float) Duration of considered sea states
    """
    q = 1 / return_period
    n = 24 / seastate_dur * 365

    beta = -1 * norm.ppf(q / n)
    n_grid, m_grid = np.mgrid[0:2 * np.pi:500j, 0:np.pi / 2:500j]
    z = beta * np.cos(m_grid)
    x = beta*np.cos(n_grid)*np.sin(m_grid)
    y = beta*np.sin(n_grid)*np.sin(m_grid)
    return x, y, z

def contour_line_for_windspeed(windspeed, return_period, seastate_dur, alpha_wind, beta_wind, a, b, i, j, k, l):
    """
    Create a contour line for a given wind speed, return period and transform it to physical space
    :param windspeed: (float) Wind speed
    :param return_period: (float) Return period for exceedance probability
    :param seastate_dur: (float) Duration of considered sea states
    :param alpha_wind: (float) Weibull shape parameter
    :param beta_wind: (float) Weibull sscale parameter
    :param a: (tuple) Intercept, slope and exponent of shape parameterization
    :param b: (tuple) Intercept, slope and exponent of scale parameterization
    :param i: (tuple) Coeff. for power func. parameterization of expected spectral peak period
    :param j: (tuple) Coeff. for power func. parameterization of expected mean wind speed
    :param k: (tuple) Coeff. for exponential func. parameterization of coeff. of variance
    :param l: (tuple) Coeff. for power func. parameterization of theta (slope)
    :return: hs: (ndarray) Sign. wave height points of circle in physical space
             tp: (ndarray) Spectral peak period points of circle in physical space
    """
    q = 1 / return_period
    n = 24 / seastate_dur * 365

    beta = -1 * norm.ppf(q / n)
    v_z = np.arccos(u_from_uw(windspeed, alpha_wind, beta_wind) / beta)
    n_grid = np.linspace(0, 2 * np.pi, 500)
    u2 = beta * np.cos(n_grid) * np.sin(v_z)
    u3 = beta * np.sin(n_grid) * np.sin(v_z)
    hs = tranform_u2(u2, windspeed, a, b)
    tp = transform_u3(u3, windspeed, hs, i, j, k, l)
    return hs, tp