import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from astropy.visualization import (ImageNormalize, LogStretch, MinMaxInterval,
                                   AsinhStretch, SqrtStretch, PercentileInterval)


plt.rcParams["font.family"] = "Times New Roman"
sns.set_context('talk')
sns.set_style('ticks')

logger = logging.getLogger('antigen.plot')

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
    sns.violinplot(
        data=df,
        x='Wavelength',
        y='Residual',
        hue='Wavelength',  # Add hue parameter
        legend=False,      # Hide the legend since it's redundant
        inner=None,
        palette='coolwarm',
        saturation=0.8
    )

    # Customize plot appearance
    plt.xlabel(r'Arc Line Wavelength ($\mathrm{\AA}$)')
    plt.ylabel(r'Measured - Expected ($\mathrm{\AA}$)')
    plt.title('Wavelength Solution Residuals')
    plt.xticks(rotation=45)

    # Save the plot as a PNG file with the given name
    plt.savefig(os.path.join(outfolder, 'wavelength_measures.png'))
    return None


def plot_spectrum_with_standard(spectrum, spectrum_error, wavelength, standard_star_wavelength,
                                standard_star_flux, response, outfolder=None, ylims=None, xlims=None):
    """
    Purpose:
        Plot a science spectrum with uncertainties, overlaid with a standard star flux.

    Args:
        spectrum (ndarray): 1D science spectrum (flux units).
        spectrum_error (ndarray): 1D uncertainties for the science spectrum.
        wavelength (ndarray): 1D wavelength array for the science spectrum.
        standard_star_wavelength (ndarray): 1D wavelength array for the standard star.
        standard_star_flux (ndarray): 1D flux array for the standard star.
        response (ndarray): 1D throughput array
        outfolder (str, optional): Directory to save PNG file. If None, plot is not saved.
        ylims (tuple, optional): y-axis limits for the flux plot.
        xlims (tuple, optional): x-axis limits for the flux plot.
    Returns:
        fig (matplotlib.figure.Figure): The matplotlib figure object.
        ax (matplotlib.axes.Axes): The matplotlib axis object.
    """
    # Interpolate standard star flux onto science wavelength grid for visual comparison
    interp_flux = interp1d(
        standard_star_wavelength,
        standard_star_flux,
        bounds_error=False,
        fill_value="extrapolate"
    )
    standard_on_science_grid = interp_flux(wavelength)
    throughput = 1./ response
    # Create figure
    fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                            gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.1})

    ax = axs[0]
    # Plot science spectrum with error bars
    ax.plot(wavelength, spectrum, color="tab:blue", lw=1.5, label="Science Spectrum")
    ax.fill_between(
        wavelength,
        spectrum - spectrum_error,
        spectrum + spectrum_error,
        color="tab:blue",
        alpha=0.3,
        label="Uncertainty"
    )

    # Plot standard star interpolated onto science grid
    ax.plot(
        wavelength,
        standard_on_science_grid,
        color="tab:red",
        lw=1.2,
        linestyle="--",
        label="Standard Star (interpolated)"
    )

    # Axis formatting
    ax.set_ylabel(r"F$_{\lambda}$ (erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}$$^{-1}$)")
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)

    ax = axs[1]
    ax.plot(
        wavelength,
        throughput,
        color="black",
        lw=1.2,
        linestyle="-",
        label="Throughput"
    )
    ax.set_xlabel("Wavelength")
    ax.set_ylabel("Throughput")
    for ax in axs:
        ax.legend()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='major', length=8, width=2)
        ax.tick_params(axis='both', which='minor', length=4, width=1)
    # Save file
    if outfolder is not None:
        Path(outfolder).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outfolder, "spectrum_with_standard.png")
        plt.savefig(outfile, dpi=200)

    return fig, ax


def plot_trace(full_trace, chunk_trace, chunk_column, fiber_indices=[5, 130, 230], outfolder=None, ylims=(-10,10)):
    """
    Purpose: Plots the residuals of the trace correction and saves the figure.

    Args:
        full_trace (np.ndarray): 2D ndarray, The array of trace correction residuals to be plotted.
        chunk_trace (np.ndarray?): undocumented
        chunk_column (np.ndarray?) undocumented
        fiber_indices (list): fiber indices to investigate the trace
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

    colors = plt.get_cmap('Set2')(np.linspace(0, 1, len(fiber_indices)))
    for fiber_index, color in zip(fiber_indices, colors):
        fiber_label = fiber_index+1
        mean_trace = np.mean(full_trace[fiber_index])
        plt.scatter(chunk_column, chunk_trace[fiber_index] - mean_trace, color='k', edgecolor='k', s=30,)
        plt.scatter(chunk_column, chunk_trace[fiber_index] - mean_trace, color=color, edgecolor='k', s=20, alpha=0.5)
        plt.plot(full_trace_X, full_trace[fiber_index] - mean_trace, color=color, lw=1, label=f'Fiber: {fiber_label}')

    plt.legend()

    # Adjust the appearance of the ticks on both axes
    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in', zorder=3)
    ax.tick_params(axis='y', which='both', left=True, right=True)
    ax.tick_params(axis='x', which='both', bottom=True, top=True)
    ax.tick_params(axis='both', which='major', length=8, width=2)
    ax.tick_params(axis='both', which='minor', length=4, width=1)
    ax.minorticks_on()
    #ax.set_ylim(ylims)

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
        logger.info(data_summary)

    return data_summary

def plot_scattered_light_fit(x_norm, y_norm, group_values, residuals, poly_order, 
                           rms_residual, n_groups, outfolder=None):
    """
    Create diagnostic plot for scattered light 2D polynomial fitting showing group positions,
    values, and fit residuals.
    
    Args:
        x_norm (np.ndarray): Normalized x coordinates of pixel groups (0-1)
        y_norm (np.ndarray): Normalized y coordinates of pixel groups (0-1)  
        group_values (np.ndarray): Background values for each pixel group
        residuals (np.ndarray): Fit residuals (observed - model) for each group
        poly_order (int): Order of 2D polynomial that was fitted
        rms_residual (float): RMS of fit residuals
        n_groups (int): Total number of pixel groups used
        outfolder (str, optional): Directory to save PNG file. If None, plot is not saved.
        
    Returns:
        fig (matplotlib.figure.Figure): The matplotlib figure object
        ax (matplotlib.axes.Axes): The matplotlib axis object
    """
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # Calculate point sizes based on absolute residuals
    # Scale residuals to reasonable point sizes (10-100 pixels)
    abs_residuals = np.abs(residuals)
    if abs_residuals.max() > 0:
        # Normalize to [0,1] then scale to [10, 100]
        size_norm = abs_residuals / abs_residuals.max()
        point_sizes = 10 + 90 * size_norm
    else:
        point_sizes = np.full_like(residuals, 30)
    
    # Create scatter plot with group values as color and residuals as size
    scatter = ax.scatter(x_norm, y_norm, 
                        c=group_values,
                        s=point_sizes,
                        alpha=0.7,
                        cmap='viridis',
                        edgecolors='black',
                        linewidths=0.5)
    
    # Add colorbar for group values
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label('Background Value (ADU)', fontsize=12)
    
    # Create legend for residual sizes
    # Show a few representative sizes
    residual_percentiles = [10, 50, 90]
    residual_values = np.percentile(abs_residuals, residual_percentiles)
    legend_sizes = []
    legend_labels = []
    
    for perc, val in zip(residual_percentiles, residual_values):
        if val > 0:
            # Calculate corresponding point size
            size_norm_val = val / abs_residuals.max() if abs_residuals.max() > 0 else 0
            legend_size = 10 + 90 * size_norm_val
            legend_sizes.append(legend_size)
            legend_labels.append(f'{perc}th percentile: {val:.1f}')
    
    # Create dummy scatter points for legend
    legend_elements = []
    for size, label in zip(legend_sizes, legend_labels):
        legend_elements.append(
            plt.scatter([], [], s=size, c='gray', alpha=0.7, 
                       edgecolors='black', linewidths=0.5, label=label)
        )
    
    # Add legend for residual sizes
    if legend_elements:
        legend = ax.legend(handles=legend_elements, 
                          title='|Residual| Size Scale',
                          loc='upper left', 
                          bbox_to_anchor=(0.02, 0.98),
                          framealpha=0.9)
        legend.get_title().set_fontsize(10)
        legend.get_title().set_fontweight('bold')
    
    # Formatting
    ax.set_xlabel('Normalized X Position', fontsize=12)
    ax.set_ylabel('Normalized Y Position', fontsize=12)
    ax.set_title(f'Scattered Light 2D Polynomial Fit Diagnostics\n'
                f'Poly Order: {poly_order}, Groups: {n_groups:,}, '
                f'RMS Residual: {rms_residual:.2f} ADU',
                fontsize=14, pad=20)
    
    # Set equal aspect ratio to show true spatial relationships
    ax.set_aspect('equal')
    
    # Add grid for reference
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Format ticks
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    ax.tick_params(axis='both', which='major', length=6, width=1.5)
    ax.tick_params(axis='both', which='minor', length=3, width=1)
    ax.minorticks_on()
    
    # Add statistics text box
    stats_text = (f'Statistics:\n'
                 f'Mean Residual: {np.mean(residuals):.2f}\n'
                 f'Std Residual: {np.std(residuals):.2f}\n'
                 f'Max |Residual|: {abs_residuals.max():.2f}\n'
                 f'Groups Used: {n_groups:,}')
    
    ax.text(0.98, 0.02, stats_text,
           transform=ax.transAxes,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
           verticalalignment='bottom',
           horizontalalignment='right',
           fontsize=10,
           fontfamily='monospace')
    
    # Tight layout
    plt.tight_layout()
    
    # Save file
    if outfolder is not None:
        Path(outfolder).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outfolder, "scattered_light_fit_diagnostics.png")
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
        # Close the figure to prevent accumulation of open figures when called repeatedly
        plt.close(fig)
        
    return fig, ax

def plot_psf_fit_diagnostics(fit_diagnostics, output_dir="."):
    """
    Create diagnostic plots for PSF fitting results.
    
    Args:
        fit_diagnostics (dict): Dictionary containing fit data from fit_psf
        output_dir (str): Directory to save plots
    """
    from pathlib import Path
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if not fit_diagnostics['chunk_data']:
        logger.warning("No chunk data available for PSF diagnostics")
        return
        
    # Create scatter plot of r vs flux values
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Collect all data across chunks
    all_r = np.concatenate([chunk['r'] for chunk in fit_diagnostics['chunk_data']])
    all_fiber_flux = np.concatenate([chunk['fiber_flux'] for chunk in fit_diagnostics['chunk_data']])
    all_model_flux = np.concatenate([chunk['model_flux'] for chunk in fit_diagnostics['chunk_data']])
    all_waves = np.concatenate([[chunk['wavelength']] * len(chunk['r']) for chunk in fit_diagnostics['chunk_data']])
    
    # Plot 1: r vs fiber flux and model flux
    ax = axes[0, 0]
    scatter1 = ax.scatter(all_r, all_fiber_flux, c='blue', s=20, alpha=0.6, label='Fiber Flux')
    scatter2 = ax.scatter(all_r, all_model_flux, c='red', s=20, alpha=0.6, marker='x', label='Model Flux')
    ax.set_xlabel('Distance from Source (arcsec)')
    ax.set_ylabel('Flux')
    ax.set_title('PSF Fit: Observed vs Model Flux vs Radius')
    ax.legend()
    
    # Plot 2: r vs fiber flux and model flux (log scale)
    ax = axes[0, 1] 
    ax.scatter(all_r, all_fiber_flux, c='blue', s=20, alpha=0.6, label='Fiber Flux')
    ax.scatter(all_r, all_model_flux, c='red', s=20, alpha=0.6, marker='x', label='Model Flux')
    ax.set_xlabel('Distance from Source (arcsec)')
    ax.set_ylabel('Flux')
    ax.set_title('PSF Fit: Observed vs Model Flux vs Radius (Log Scale)')
    ax.set_yscale('log')
    ax.legend()
    
    # Extract source positions and FWHM from chunks
    wavelengths = fit_diagnostics['wavelengths']
    source_x = np.array([chunk['fitted_x'] for chunk in fit_diagnostics['chunk_data']])
    source_y = np.array([chunk['fitted_y'] for chunk in fit_diagnostics['chunk_data']])
    source_fwhm = np.array([chunk['fitted_fwhm'] for chunk in fit_diagnostics['chunk_data']])
    
    # Plot 3: Wavelength evolution of source position
    ax = axes[1, 0]
    ax.plot(wavelengths, source_x, 'ko-', label='X position', markersize=4)
    ax.plot(wavelengths, source_y, 'go-', label='Y position', markersize=4)
    ax.set_xlabel('Wavelength (Å)')
    ax.set_ylabel('Position (arcsec)')
    ax.set_title('Source Position vs Wavelength')
    ax.legend()
    
    # Plot 4: FWHM vs wavelength
    ax = axes[1, 1]
    ax.plot(wavelengths, source_fwhm, 'go-', markersize=4)
    ax.set_xlabel('Wavelength (Å)')
    ax.set_ylabel('FWHM (arcsec)')
    ax.set_title('Seeing vs Wavelength')
    
    plt.tight_layout()
    plt.savefig(output_path / 'psf_fit_diagnostics.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create per-chunk detailed plots
    fig, axes = plt.subplots(2, min(5, len(fit_diagnostics['chunk_data'])), figsize=(15, 8))
    if len(fit_diagnostics['chunk_data']) == 1:
        axes = axes.reshape(2, 1)
    
    for i, chunk in enumerate(fit_diagnostics['chunk_data'][:5]):  # Show first 5 chunks
        if len(fit_diagnostics['chunk_data']) > 1:
            ax1, ax2 = axes[0, i], axes[1, i]
        else:
            ax1, ax2 = axes[0, 0], axes[1, 0]
            
        # Scatter plot in x-y plane
        scatter = ax1.scatter(chunk['x_abs'], chunk['y_abs'], c=chunk['fiber_flux'], 
                            s=30, cmap='viridis')
        ax1.plot(chunk['fitted_x'], chunk['fitted_y'], 'r*', markersize=15)
        ax1.set_xlabel('X (arcsec)')
        ax1.set_ylabel('Y (arcsec)')
        ax1.set_title(f'Chunk {i+1} λ={chunk["wavelength"]:.0f}Å')
        ax1.set_aspect('equal')
        
        # Radial profile with both fiber flux and model flux
        ax2.scatter(chunk['r'], chunk['fiber_flux'], alpha=0.7, label='Fiber Flux', color='blue')
        ax2.scatter(chunk['r'], chunk['model_flux'], alpha=0.7, label='Model Flux', color='red', marker='x')
        ax2.set_xlabel('Distance (arcsec)')
        ax2.set_ylabel('Flux')
        ax2.set_title(f'Radial Profile')
        ax2.legend()
        
    plt.tight_layout()
    plt.savefig(output_path / 'psf_fit_chunks.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PSF fit diagnostics saved to {output_path}")


def plot_dar_models_comparison(wave, dx_psf, dy_psf, dx_hdr, dy_hdr, output_dir=".", align_to_ref=True, ref_wave=5000.0):
    """Plot comparison between DAR models derived from PSF-fit coefficients and from header metadata.

    To examine models consistently, this function can normalize (shift) both
    measured and modeled DAR curves so that they share the same zero-point at a
    common reference wavelength (by default 5000 Å, matching the modeled DAR's
    typical ref_wave). This removes any arbitrary constant offset between the
    two dx/dy curves due to different normalization choices.

    Args:
        wave (ndarray): Wavelength grid (Angstroms).
        dx_psf (ndarray): DAR dx from PSF-derived model (arcsec).
        dy_psf (ndarray): DAR dy from PSF-derived model (arcsec).
        dx_hdr (ndarray): DAR dx from header-derived model (arcsec).
        dy_hdr (ndarray): DAR dy from header-derived model (arcsec).
        output_dir (str, optional): Directory to save plot. Defaults to ".".
        align_to_ref (bool, optional): If True, align both models to share the
            same zero-point at ref_wave (or the closest available wavelength if
            ref_wave lies outside the grid). Defaults to True.
        ref_wave (float, optional): Reference wavelength (Å) used for alignment.
            Defaults to 5000.0.
    """
    from pathlib import Path

    wave = np.asarray(wave, dtype=float)
    dx_psf = np.asarray(dx_psf, dtype=float)
    dy_psf = np.asarray(dy_psf, dtype=float)
    dx_hdr = np.asarray(dx_hdr, dtype=float)
    dy_hdr = np.asarray(dy_hdr, dtype=float)

    # Determine anchor wavelength for alignment
    anchor_wave = None
    if align_to_ref and np.isfinite(ref_wave):
        wmin, wmax = np.nanmin(wave), np.nanmax(wave)
        if (ref_wave >= wmin) and (ref_wave <= wmax):
            anchor_wave = float(ref_wave)
        elif (5000.0 >= wmin) and (5000.0 <= wmax):
            anchor_wave = 5000.0
        else:
            # Fallback: pick the wavelength where the modeled dx is closest to zero
            try:
                idx = int(np.nanargmin(np.abs(dx_hdr)))
                anchor_wave = float(wave[idx])
            except Exception:
                # Last resort: center of wavelength range
                anchor_wave = float(0.5 * (wmin + wmax))

    # If aligning, subtract value at anchor wavelength from each curve
    def _at_anchor(y, x, x0):
        # Linear interpolation for value at anchor
        if x0 <= x[0]:
            return float(y[0])
        if x0 >= x[-1]:
            return float(y[-1])
        i = np.searchsorted(x, x0)  # index where x[i-1] < x0 <= x[i]
        i = max(1, min(i, len(x) - 1))
        x1, x2 = x[i-1], x[i]
        y1, y2 = y[i-1], y[i]
        t = (x0 - x1) / (x2 - x1) if np.isfinite(x2 - x1) and (x2 != x1) else 0.0
        return float(y1 + t * (y2 - y1))

    if anchor_wave is not None:
        dx_psf0 = _at_anchor(dx_psf, wave, anchor_wave)
        dy_psf0 = _at_anchor(dy_psf, wave, anchor_wave)
        dx_hdr0 = _at_anchor(dx_hdr, wave, anchor_wave)
        dy_hdr0 = _at_anchor(dy_hdr, wave, anchor_wave)

        dx_psf_al = dx_psf - dx_psf0
        dy_psf_al = dy_psf - dy_psf0
        dx_hdr_al = dx_hdr - dx_hdr0
        dy_hdr_al = dy_hdr - dy_hdr0
    else:
        dx_psf_al, dy_psf_al = dx_psf, dy_psf
        dx_hdr_al, dy_hdr_al = dx_hdr, dy_hdr

    outdir = Path(output_dir) / "dar_diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    # dx comparison
    ax = axes[0, 0]
    ax.plot(wave, dx_psf_al, label="dx (PSF)", color="tab:blue")
    ax.plot(wave, dx_hdr_al, label="dx (Header)", color="tab:orange", linestyle="--")
    ax.set_ylabel("dx (arcsec)")
    title = "DAR dx vs wavelength"
    if anchor_wave is not None:
        ax.axvline(anchor_wave, color="gray", linestyle=":", alpha=0.6)
        title += f" (aligned at {anchor_wave:.0f} Å)"
    ax.set_title(title)
    ax.legend()

    # dy comparison
    ax = axes[0, 1]
    ax.plot(wave, dy_psf_al, label="dy (PSF)", color="tab:green")
    ax.plot(wave, dy_hdr_al, label="dy (Header)", color="tab:red", linestyle="--")
    ax.set_ylabel("dy (arcsec)")
    title = "DAR dy vs wavelength"
    if anchor_wave is not None:
        ax.axvline(anchor_wave, color="gray", linestyle=":", alpha=0.6)
        title += f" (aligned at {anchor_wave:.0f} Å)"
    ax.set_title(title)
    ax.legend()

    # residuals (after alignment if applied)
    ax = axes[1, 0]
    ax.plot(wave, dx_psf_al - dx_hdr_al, color="tab:purple")
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("dx residual (PSF - Header) [arcsec]")

    ax = axes[1, 1]
    ax.plot(wave, dy_psf_al - dy_hdr_al, color="tab:brown")
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("dy residual (PSF - Header) [arcsec]")

    plt.tight_layout()
    plt.savefig(outdir / "dar_models_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_dither_summary(dither_info, output_dir="."):
    """
    Create summary plots for dither pattern analysis.
    
    Args:
        dither_info (dict): Dictionary containing dither information
        output_dir (str): Output directory for plots
    """

    if dither_info is None:
        logger.info("No dither information to plot")
        return
        
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Dither pattern
    ax = axes[0]
    offsets = np.array(dither_info['dither_offsets'])
    ax.scatter(offsets[:, 0], offsets[:, 1], s=100, c=range(len(offsets)), 
               cmap='viridis', edgecolors='black', linewidths=2)
    
    # Add dither numbers
    for i, (dx, dy) in enumerate(offsets):
        ax.annotate(f'{i+1}', (dx, dy), xytext=(5, 5), textcoords='offset points',
                   fontweight='bold', color='white')
    
    ax.set_xlabel('X Offset (arcsec)')
    ax.set_ylabel('Y Offset (arcsec)')
    ax.set_title('Detected Dither Pattern')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Plot 2: Offset magnitudes
    ax = axes[1]
    distances = np.sqrt(offsets[:, 0]**2 + offsets[:, 1]**2)
    ax.bar(range(1, len(distances)+1), distances, color='skyblue', edgecolor='navy')
    ax.set_xlabel('Dither Number')
    ax.set_ylabel('Offset from First Dither (arcsec)')
    ax.set_title('Dither Offset Magnitudes')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path / 'dither_summary.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Dither summary saved to {output_path}")
    
    # Print summary table
    logger.info("\nDither Summary:")
    logger.info("Dither  X_offset  Y_offset  Distance")
    logger.info("------  --------  --------  --------")
    for i, ((dx, dy), dist) in enumerate(zip(offsets, distances)):
        logger.info(f"  {i+1:2d}    {dx:6.2f}    {dy:6.2f}    {dist:6.2f}")


def plot_fiber_flux_distribution(fiber_x, fiber_y, reduced_spectra, fiber_radius,
                                 output_dir="."):
    """
    Create a scatter plot showing the distribution of light across fibers.

    Args:
        fiber_x (ndarray): X positions of fibers
        fiber_y (ndarray): Y positions of fibers
        reduced_spectra (ndarray): Reduced spectra with shape (Nfibers, Nlambda)
        fiber_radius (float): Radius of fiber in arcsec
        output_dir (str, optional): Path to save the plot. If None, shows plot.
    """
    # Collapse flux over wavelength
    collapsed_flux = np.nanmedian(reduced_spectra, axis=1)

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Convert fiber radius to marker size (this is approximate and may need adjustment)
    fiber_diameter = 2 * fiber_radius

    # Get the axis limits to calculate appropriate scaling
    x_range = np.ptp(fiber_x)  # peak-to-peak (max - min)
    y_range = np.ptp(fiber_y)
    plot_range = max(x_range, y_range)

    # Scale marker size relative to plot range
    # This creates markers that are proportional to the actual fiber size
    marker_size = (fiber_diameter / plot_range * 500) ** 2# 500 is a scaling factor

    # Create scatter plot with flux as color and size
    valid_mask = np.isfinite(collapsed_flux)
    all_indices = np.arange(len(collapsed_flux))
    interval_func = PercentileInterval(99.9)
    stretch_func = LogStretch()



    if np.any(valid_mask):
        valid_flux = collapsed_flux[valid_mask]
        norm = ImageNormalize(valid_flux, interval=interval_func, stretch=stretch_func)

        scatter = ax.scatter(fiber_x[valid_mask], fiber_y[valid_mask],
                           c=collapsed_flux[valid_mask], norm=norm,
                           s=marker_size, cmap='viridis', alpha=0.7)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Median Flux', rotation=270, labelpad=20)

        # Get the valid indices
        valid_indices = all_indices[valid_mask]

        # Mark the brightest fiber
        brightest_idx = np.nanargmax(collapsed_flux[valid_mask])
        brightest_fiber_idx = valid_indices[brightest_idx]
        ax.plot(fiber_x[brightest_fiber_idx], fiber_y[brightest_fiber_idx],
                'r*', markersize=15, label=f'Brightest fiber (#{brightest_fiber_idx})')

    # Plot all fiber positions (including low/zero flux ones)
    ax.scatter(fiber_x[~valid_mask], fiber_y[~valid_mask],
               c='lightgray', s=marker_size, alpha=0.3, label='No/invalid flux')


    # Formatting
    ax.set_xlabel('Fiber X Position (arcsec)')
    ax.set_ylabel('Fiber Y Position (arcsec)')
    ax.set_title('Distribution of Light Across Fibers\n(Median flux over wavelength)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Add some statistics
    if np.any(valid_mask):
        n_with_flux = np.sum(valid_mask)
        max_flux = np.nanmax(collapsed_flux[valid_mask])
        mean_flux = np.nanmean(collapsed_flux[valid_mask])

        stats_text = f'Fibers with flux: {n_with_flux}/{len(fiber_x)}\n'
        stats_text += f'Max flux: {max_flux:.2e}\n'
        stats_text += f'Mean flux: {mean_flux:.2e}'

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    # Save file
    if output_dir is not None:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(output_dir, "flux_distribution_fibers.png")
        plt.savefig(outfile, dpi=200, bbox_inches='tight')
