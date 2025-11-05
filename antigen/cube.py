import numpy as np
import logging

from scipy.interpolate import griddata, Rbf
from scipy.spatial import cKDTree

from antigen import fiber

logger = logging.getLogger('antigen.cube')

def fibers_to_image(fiber_x, fiber_y, fiber_flux, fiber_area, bounds=[-10., 10., -10., 10],
                    pixel_size=1.0, method="linear", rbf_func="multiquadric", k=5,
                    sigma=1.0):
    """
    Create a synthetic image from non-uniform fiber locations and flux values.

    Args:
        fiber_x (ndarray): 1D array of fiber x-locations.
        fiber_y (ndarray): 1D array of fiber y-locations.
        fiber_flux (ndarray): 1D array of fiber flux values.
        fiber_area (float): Area of fiber in square arcseconds.
        bounds (list): Boundaries of image.
        pixel_size (float): Pixel size in arcseconds for the output image.
        method (str): interpolation method:
            - "nearest" : nearest-neighbor (griddata)
            - "linear"  : linear interpolation (griddata)
            - "cubic"   : cubic interpolation (griddata)
            - "rbf"     : radial basis function interpolation
            - "gdw"     : inverse-distance weighting
        rbf_func (str): radial basis function (for method="rbf").
        k (int): number of neighbors for GDW.
        sigma (float): standard deviation for GDW.

    Returns:
        img (ndarray): synthetic image on uniform grid
        X, Y (ndarray): meshgrid coordinates
    """
    # Clean NaNs in flux and align x/y/flux
    finite_values = np.isfinite(fiber_flux)
    x, y, flux = (fiber_x[finite_values], fiber_y[finite_values], fiber_flux[finite_values])

    # Define grid
    grid_size_x = int((bounds[1] - bounds[0]) / pixel_size) + 1
    grid_size_y = int((bounds[3] - bounds[2]) / pixel_size) + 1

    xi = np.linspace(bounds[0], bounds[0] + (grid_size_x - 1) * pixel_size, grid_size_x)
    yi = np.linspace(bounds[2], bounds[2] + (grid_size_y - 1) * pixel_size, grid_size_y)
    X, Y = np.meshgrid(xi, yi)

    # If no valid input points, return zeros image
    if x.size == 0:
        img = np.zeros_like(X, dtype=float)
        # Correct for the area difference between the fiber size and the new pixel area
        flux_factor = pixel_size ** 2 / fiber_area
        img = img * flux_factor
        return img, X, Y

    flux_factor = pixel_size ** 2 / fiber_area

    if method in ["linear", "nearest", "cubic"]:
        # griddata can fail for too few points or outside convex hull; fall back to nearest
        try:
            img = griddata((x, y), flux, (X, Y), method=method)
        except Exception:
            img = griddata((x, y), flux, (X, Y), method="nearest")

    elif method == "rbf":
        # RBF interpolation
        rbf = Rbf(x, y, flux, function=rbf_func)
        img = rbf(X, Y)

    elif method == "gdw":
        # Gaussian Distance Weighting (IDW-style)
        npts = x.size
        k_eff = min(k, npts)
        if k_eff == 0:
            img = np.zeros_like(X, dtype=float)
        else:
            tree = cKDTree(np.c_[x, y])
            dist, idx = tree.query(np.c_[X.ravel(), Y.ravel()], k=k_eff)

            # Ensure 2D shapes for k_eff == 1
            if k_eff == 1:
                dist = dist[:, None]
                idx = idx[:, None]

            # Gaussian weights based on seeing sigma
            weights = np.exp(-0.5 * (dist / sigma) ** 2)

            # Normalize rows safely
            row_sums = weights.sum(axis=1, keepdims=True)
            with np.errstate(invalid='ignore', divide='ignore'):
                weights = np.where(row_sums > 0, weights / row_sums, 0.0)

            # Weighted sum of flux
            img = np.sum(weights * flux[idx], axis=1).reshape(X.shape)

    else:
        raise ValueError(f"Unknown method {method}")

    # Correct for the area difference between the fiber size and the new pixel area
    img = img * flux_factor
    return img, X, Y

def make_cube(wavelength, fiber_spectra, fiber_x, fiber_y, fiber_area,
              modeling=None, pixel_size=1.0, method="linear",
              rbf_func="multiquadric", k=5, sigma=1.0):
    """Construct a 3D datacube from fiber spectra and positions.

    This function generates a spectral cube by mapping fiber spectra onto
    a 2D grid at each wavelength. Differential atmospheric refraction (DAR)
    can be applied to adjust fiber positions as a function of wavelength.

    Args:
        wavelength (ndarray): 1D array of wavelengths.
        fiber_spectra (ndarray): 2D array of fiber flux.
        fiber_x (ndarray): 1D array of fiber x-locations.
        fiber_y (ndarray): 1D array of fiber y-locations.
        fiber_area (float): Area of fiber in square arcseconds.
        modeling (dict, optional): Dictionary containing PSF and DAR modeling results with keys:
            - dar_model: DAR model for position correction
            - sources: Detected source catalog  
            - X: Grid of X coordinates
            - Y: Grid of Y coordinates
            If None, no DAR correction is applied.
        pixel_size (float): Pixel size in arcseconds (Default is 1.0).
        method (str, optional): interpolation method
        rbf_func (str, optional): Radial basis function type if `method="rbf"`. Default is "multiquadric".
        k (int, optional): Number of nearest neighbors used if `method="idw"`. Default is 5.
        sigma (float, optional): standard deviation for GDW.

    Returns:
        cube (ndarray):
            3D datacube with shape (N_lambda, Ny, Nx), where `Ny` and `Nx`
            are determined by the spatial extent of the fibers and `pixel_size`.
        x_grid (ndarray): 2D array of x-locations.
        y_grid (ndarray): 2D array of y-locations.
    """
    # Precompute output grid extent based on fiber positions
    x_min, x_max = np.min(fiber_x), np.max(fiber_x)
    y_min, y_max = np.min(fiber_y), np.max(fiber_y)

    # Estimate output dimensions
    nx = int((x_max - x_min) / pixel_size) + 1
    ny = int((y_max - y_min) / pixel_size) + 1

    # Determine spectral dimensions and guard against mismatches
    n_wave = int(len(wavelength))
    n_spec = int(fiber_spectra.shape[1])
    n_lambda = min(n_wave, n_spec)
    if n_wave != n_spec:
        logger.warning(
            "Wavelength grid length (%d) != spectra columns (%d). Proceeding with %d planes. "
            "This often indicates an incorrect binning setting; check the --binned flag in antigen_make_cubes.py.",
            n_wave, n_spec, n_lambda
        )

    # Allocate cube
    cube = np.zeros((n_lambda, ny, nx), dtype=float)
    bounds = fiber.get_fiber_bounds(fiber_x, fiber_y)

    # Loop over wavelength bins
    for i, lam in enumerate(wavelength[:n_lambda]):
        flux = fiber_spectra[:, i]

        if modeling is not None and 'dar_model' in modeling:
            # Apply DAR model directly to each fiber position at this wavelength
            x_shift, y_shift = modeling['dar_model'](lam, fiber_x, fiber_y)
        else:
            x_shift, y_shift = fiber_x, fiber_y

        image, x_grid, y_grid = fibers_to_image(
            x_shift, y_shift, flux, fiber_area, bounds=bounds,
            pixel_size=pixel_size, method=method,
            rbf_func=rbf_func, k=k, sigma=sigma,
        )
        cube[i] = image

    return cube, x_grid, y_grid