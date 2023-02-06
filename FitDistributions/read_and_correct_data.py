import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(rc={'text.usetex' : True})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.5)


def bin_data(data, width):
    """
    Bins data in a numpy arrays into bins of specified width. Bins start at zero and goes until max of data.
    :param data: (ndarray) Input data
    :param width: (float) Width of bins
    :return: bin_height: (ndarray, size N) Number of values in bin.
             bin_edges: (ndarray, size N+1) Edges of bins, including leftmost and rightmost edge.
             mid_bin: (ndarray, size N) Point between left and right edge in each bin.
    """
    bin_height, bin_edge = np.histogram(data, bins=np.arange(0, np.max(data)+width, step=width))
    mid_bin = (bin_edge[1:]-bin_edge[:-1])/2 + bin_edge[:-1]
    return bin_height, bin_edge, mid_bin


def tp_correction(Tp, Hs, plot=False, save_plot=False, show_plot=False):
    """
    Function to correct log-spacing in NORA-dataset
    :param Tp: (ndarray, size N) Uncorrected Tp
    :param Hs: (ndarray, size N) Sign. wave height, for scatter plot
    :param plot: (bool) Flag to decide plotting
    :param save_plot: (bool) Flag to decide saving of plot
    :param show_plot: (bool) Flag to decide showing of plot
    :return: Tp_modified (ndarray, size N) Corrected Tp
    """
    Tp_modified = np.zeros_like(Tp)  # Allocating
    sorted_unique_Tp = np.sort(np.unique(Tp))  # Allocating

    """
    Removing instances of infinite frequencies (= Periods of zero)
    """
    if sorted_unique_Tp[0] == 0:
        min_non_zero_Tp = sorted_unique_Tp[1]
    else:
        min_non_zero_Tp = sorted_unique_Tp[0]

    for j in range(len(Tp)):
        """
        Modifying Tp as in App. D (Haver)
        """
        if Tp[j] == 0:
            i = np.around(1+np.log(min_non_zero_Tp/3.244)/0.09525)
        else:
            i = np.around(1+np.log(Tp[j]/3.244)/0.09525)
        Tp_modified[j] = 3.244 * np.exp(0.09525*(i-0.5-np.random.uniform()))

    if plot:
        fig = plt.figure()
        plt.scatter(Tp_modified, Hs, c="k", s=1, label="Corrected")
        plt.scatter(Tp, Hs, c="r", s=1, label="Raw data")
        lgnd = plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
        for handle in lgnd.legendHandles:
            handle.set_sizes([20.0])
        plt.xlabel(r"$T_{p} [s]$")
        plt.ylabel(r"$H_{s} [m]$")
        plt.title("Correction of peak spectral period")
        plt.tight_layout()
        if save_plot:
            plt.savefig("tp_correction.png")
        if show_plot:
            plt.show()
        else:
            plt.close(fig)
    return Tp_modified