import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

def plot_QTF_contour(ax, qtf, period1, period2):
    P1, P2 = np.meshgrid(period1, period2)
    cp = ax.contourf(P1, P2, qtf)
    return cp

def plot_QTF_surface(ax, qtf, period1, period2):
    P1, P2 = np.meshgrid(period1, period2)
    cp = ax.plot_surface(P1, P2, qtf, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    return cp