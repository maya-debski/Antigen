import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
sns.set_context('talk')
sns.set_style('ticks')


def plot_wavelength(lines, W, wavelength, outfolder=None):
    """
    Purpose: Plots the residuals of the wavelength solution using a violin plot.

    Args:
        lines : 1D ndarray, Expected wavelengths of arc lamp lines.
        W : 2D ndarray, Measured positions of the arc lines for different fibers.
        wavelength : 2D ndarray, Wavelength solution for the fibers.

    Returns:
        None
    """

    # Prepare data for seaborn violin plot
    data = []
    num_fibers = W.shape[0]
    for i, line in enumerate(lines):
        residuals = []
        for fiber in range(2, num_fibers): # skip the first two broken fibers
            X = np.arange(wavelength.shape[1])
            pred = np.interp(W[fiber, i], X, wavelength[fiber])
            if np.all(wavelength[fiber] != 0.):
                residuals.append(line - pred)
        for res in residuals:
            data.append({'Wavelength': line, 'Residual': res})

    df = pd.DataFrame(data)

    # Create the violin plot
    plt.figure(figsize=(14, 7))
    plt.gca().set_position([0.15, 0.19, 0.75, 0.71])
    sns.violinplot(x='Wavelength', y='Residual', data=df, palette='coolwarm', inner=None, saturation=0.8)
    # Customize plot appearance
    plt.xlabel(r'Arc Line Wavelength ($\mathrm{\AA}$)')
    plt.ylabel(r'Measured - Expected ($\mathrm{\AA}$)')
    plt.title('Wavelength Solution Residuals')
    plt.xticks(rotation=45)

    # Save the plot as a PNG file with the given name
    plt.savefig(os.path.join(outfolder, 'wavelength_measures.png'))
    return None


def plot_trace(full_trace, trace, x, orders=[5, 130, 230], outfolder=None):
    """
    Purpose: Plots the residuals of the trace correction and saves the figure.

    Args:
        trace_cor (np.ndarray): 2D ndarray, The array of trace correction residuals to be plotted.
        full_trace (?): undocumented
        trace (?): undocumented
        x (?): undocumented
        orders (list): undocumented
        outfolder (str): directory to save PNG file

    Returns:
        None
    """
    # TODO: update docstring to match input args, function signature
    # TODO: replace 'x' and 'X' with better identifiers, for search, for human readers

    X = np.arange(full_trace.shape[1])
    # Create a figure with specified size
    plt.figure(figsize=(8, 7))
    colors = plt.get_cmap('Set2')(np.linspace(0, 1, len(orders)))
    for order, color in zip(orders, colors):
        mean_trace = np.mean(full_trace[order])
        plt.scatter(x, trace[order] - mean_trace, color='k', edgecolor='k', s=30,)
        plt.scatter(x, trace[order] - mean_trace, color=color, edgecolor='k', s=20, alpha=0.5)
        plt.plot(X, full_trace[order] - mean_trace, color=color, lw=1, label='Fiber: %i' % (order+1))

    plt.legend()

    # Adjust the appearance of the ticks on both axes
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', zorder=3)
    ax.tick_params(axis='y', which='both', left=True, right=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=True)
    ax.tick_params(axis='both', which='major', length=8, width=2)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.minorticks_on()

    # Label the axes
    plt.xlabel('Column')
    plt.ylabel('Trace - Mean(Trace)')

    # Save the plot as a PNG file with the given name
    plt.savefig(os.path.join(outfolder, 'trace_measures.png'))
    return None
