import numpy as np
import matplotlib.pyplot as plt

def power_function(x, a1, a2, a3):
    return a1+a2*(x**a3)

def linear_function(x, a1, a2):
    return a1+a2*x

def exponential(x, a1, a2, a3):
    return a1+a2*np.exp(x*a3)

def sigmoid(x, a1, a2, a3, a4):
    return a1+a2*(1/(1+np.exp(-a4*x+a3)))

