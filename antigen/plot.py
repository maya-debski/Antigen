import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.cm as mpl_cm

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


def plot_trace(full_trace, chunk_trace, chunk_column, orders=[5, 130, 230], outfolder=None, ylims=(-10,10)):
    """
    Purpose: Plots the residuals of the trace correction and saves the figure.

    Args:
        full_trace (np.ndarray): 2D ndarray, The array of trace correction residuals to be plotted.
        chunk_trace (np.ndarray?): undocumented
        chunk_column (np.ndarray?) undocumented
        orders (list): 'orders' of what? fibers? undocumented
        outfolder (str): directory to save PNG file

    Returns:
        fig (matplotlib.figure.Figure)
        ax (matplotlib.axes.Axes)
    """
    # TODO: update docstring to match input args, function signature

    full_trace_X = np.arange(full_trace.shape[1])
    # Create a figure with specified size
    plt.figure(figsize=(8, 7))
    fig = plt.gcf()  # Get Current Figure

    colors = plt.get_cmap('Set2')(np.linspace(0, 1, len(orders)))
    for order, color in zip(orders, colors):
        fiber_label = order+1
        mean_trace = np.mean(full_trace[order])
        plt.scatter(chunk_column, chunk_trace[order] - mean_trace, color='k', edgecolor='k', s=30,)
        plt.scatter(chunk_column, chunk_trace[order] - mean_trace, color=color, edgecolor='k', s=20, alpha=0.5)
        plt.plot(full_trace_X, full_trace[order] - mean_trace, color=color, lw=1, label=f'Fiber: {fiber_label}')

    plt.legend()

    # Adjust the appearance of the ticks on both axes
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', zorder=3)
    ax.tick_params(axis='y', which='both', left=True, right=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=True)
    ax.tick_params(axis='both', which='major', length=8, width=2)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.minorticks_on()
    ax.set_ylim(ylims)

    # Label the axes
    plt.xlabel('Column')
    plt.ylabel('Trace - Mean(Trace)')

    # Save the plot as a PNG file with the given name
    plt.savefig(os.path.join(outfolder, 'trace_measures.png'))

    return fig, ax


def plot_frame(data_array, save_file=None, title=None, alter_method=None):
    """
    Purpose: Take a single ndarray, and plot a single heatmap
    Note: Input data_array could be the Archon frame buffer (1x4, all 4 CCDs)
          or just a single channel_array (1x1, single CCD read-out)

    Args:
        data_array (np.ndarray): frame data as numpy.ndarray, read from fits.data
        save_file (str): default=None, file path for save, will not save if name is None or file already exists
        title (str): default=None, for matplotlib ax.set_title()
        alter_method (str): default=None
            option "log", "ln" rescale all data_array values
            option "clip" cuts off data_array values outside 5th-95th percentiles
    Returns:
        fig (matplotlib.figure.Figure)
        ax (matplotlib.axes.Axes)
    """
    if alter_method:
        data_array, alter_desc = alter_array_data(data_array, alter_method)
        if title:
            title += f' \n {alter_desc}'
        else:
            title = f' Exposure Data \n {alter_desc} '

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    cmap_name = 'viridis'
    img_obj = ax.imshow(data_array, cmap=cmap_name, origin='lower', aspect='auto')
    cbar_ax = make_colomap_axes(fig, ax)
    _ = fig.colorbar(img_obj, cax=cbar_ax)

    if title:
        ax.set_title(title)

    if save_file:
        if os.path.isfile(save_file):
            raise FileExistsError(f'ERROR: File name already exists. Will not over-write: {save_file}')
        else:
            fig.savefig(save_file)
    return fig, ax


def alter_array_data(data_array, alter_method=None):
    """
    Purpose: Alter the 2D array data of a frame for purposes of generating debug plots for human visual inspection
    Methods provided are focused on mitigating outliers washing out the colormap so that the data trends can be seen

    Args:
        data_array (np.ndarray): 2D numpy array of frame data numeric values
        alter_method (str): default=None
            option "log", "ln" rescale all data_array values
            option "clip" cuts off data_array values outside 5th-95th percentiles
    """
    if alter_method == 'log':
        altered_data = np.log10(data_array)
        description = 'Altered: LOG10(DATA)'  # np.log() is base-e, not base-10
    elif alter_method == 'ln':
        altered_data = np.log(data_array)
        description = 'Altered: LN(DATA)'
    elif alter_method == 'clip':
        lower_percentile = 5
        upper_percentile = 95
        lower_bound = np.percentile(data_array, lower_percentile)
        upper_bound = np.percentile(data_array, upper_percentile)
        altered_data = np.clip(data_array, lower_bound, upper_bound)
        description = 'Altered: CLIP at [5th, 95th] Percentiles'
    else:
        error_msg = ('Supported alter_method values = "log", "ln", "clip", '
                     f'Received unsupported alter_method={alter_method}')
        raise ValueError(error_msg)
    return altered_data, description


def make_colomap_axes(fig, ax):
    """
    Purpose: Fix the awful sizing/positioning of default colorbar in matplotlib

    Args:
        fig (matplotlib.figure.Figure)
        ax (matplotlib.axes.Axes)
    Returns:
        cbar_ax (): matplotlib colorbar axes object
    """
    cbar_pad = 0.01
    cbar_x0 = ax.get_position().x1 + cbar_pad
    cbar_y0 = ax.get_position().y0
    cbar_dx = 2*cbar_pad
    cbar_dy = ax.get_position().height
    cbar_coordinates = (cbar_x0, cbar_y0, cbar_dx, cbar_dy)
    cbar_ax = fig.add_axes(cbar_coordinates)
    return cbar_ax


def get_array_data_summary(data_array, verbose=False):
    """
    Purpose: get summary stats of data_array values, analogous to pandas.DataFrame.describe()

    Args:
        data_array (np.ndarray): 2D numpy array read from FITS file
        verbose (bool): default=Fault, if True, print(data_summary)
    Returns
        array_summary (dict): summary stats from numpy methods stored as a dict with
                              keys = ['shape', 'dtype', 'min', 'max', 'median', 'std']
    """
    # convert from numpy types to python primatives as prep for JSON serialization
    ptiles_ints = [5, 25, 50, 75, 95]
    percentile_vals = np.percentile(data_array, ptiles_ints)
    ptiles_keys = [str(i)+'th' for i in ptiles_ints]
    percentiles_dict = dict(zip(ptiles_keys, percentile_vals))

    data_summary = {
        'shape': data_array.shape,
        'dtype': str(data_array.dtype),
        'min': float(np.min(data_array)),
        'max': float(np.max(data_array)),
        'mean': float(np.mean(data_array)),
        'std': float(np.std(data_array)),
        **percentiles_dict
    }
    if verbose:
        print(data_summary)

    return data_summary
