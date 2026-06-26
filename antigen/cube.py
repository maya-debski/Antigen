import numpy as np
import logging

from scipy.interpolate import griddata, Rbf
from scipy.spatial import cKDTree

from antigen import fiber

logger = logging.getLogger('antigen.cube')

def fibers_to_image(fiber_x, fiber_y, fiber_flux, fiber_area, bounds=[-10., 10., -10., 10],
                    pixel_size=1.0, method="linear", rbf_func="multiquadric", k=5,
                    sigma=1.0, fiber_error=None, return_error=False):
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
        (optionally) err_img (ndarray): synthetic error image if return_error=True
    """
    # Coerce to 1D arrays
    fx = np.asarray(fiber_x).ravel()
    fy = np.asarray(fiber_y).ravel()
    ff = np.asarray(fiber_flux).ravel()

    # Align lengths to common minimum to avoid broadcasting issues
    lengths = [fx.size, fy.size, ff.size]
    if return_error and fiber_error is not None:
        fe = np.asarray(fiber_error).ravel()
        lengths.append(fe.size)
    else:
        fe = None
    n_common = int(min(lengths))
    if not (fx.size == fy.size == ff.size and (fe is None or fe.size == n_common)):
        logger.warning(
            "Length mismatch in fibers_to_image: x=%d, y=%d, flux=%d%s; truncating to %d",
            fx.size, fy.size, ff.size,
            " , err=" + str(fe.size) if fe is not None else "",
            n_common,
        )
    fx = fx[:n_common]; fy = fy[:n_common]; ff = ff[:n_common]
    if fe is not None:
        fe = fe[:n_common]

    # Mask non-finite entries across x, y, flux
    m = np.isfinite(fx) & np.isfinite(fy) & np.isfinite(ff)
    bad_x = int(np.count_nonzero(~np.isfinite(fx)))
    bad_y = int(np.count_nonzero(~np.isfinite(fy)))
    bad_f = int(np.count_nonzero(~np.isfinite(ff)))
    n_all = fx.size
    n_good = int(np.count_nonzero(m))
    n_bad_any = n_all - n_good
    if (bad_x + bad_y + bad_f) > 0:
        frac_bad = (n_bad_any / n_all) if n_all > 0 else 0.0
        msg = ("Non-finite entries detected in fibers_to_image: bad_x=%d, bad_y=%d, bad_flux=%d "
               "(%.1f%% of points affected); using finite subset.")
        if frac_bad >= 0.5:
            logger.warning(msg, bad_x, bad_y, bad_f, 100.0 * frac_bad)
        else:
            logger.debug(msg, bad_x, bad_y, bad_f, 100.0 * frac_bad)
    fx = fx[m]; fy = fy[m]; ff = ff[m]
    if return_error:
        if fiber_error is None:
            raise ValueError("return_error=True but fiber_error is None")
        fiber_var = np.square(fe[m]) if fe is not None else None

    # Validate pixel_size
    if not np.isfinite(pixel_size) or pixel_size <= 0:
        logger.warning("Invalid pixel_size=%s; defaulting to 1.0", str(pixel_size))
        pixel_size = 1.0

    # Validate/repair bounds if needed
    b = np.asarray(bounds, dtype=float).ravel()
    if b.size != 4 or not np.all(np.isfinite(b)) or (b[1] <= b[0]) or (b[3] <= b[2]):
        if fx.size > 0:
            xmin = float(np.nanmin(fx)); xmax = float(np.nanmax(fx))
            ymin = float(np.nanmin(fy)); ymax = float(np.nanmax(fy))
            if not (np.isfinite(xmin) and np.isfinite(xmax) and xmin < xmax):
                xmin, xmax = -5.0, 5.0
            if not (np.isfinite(ymin) and np.isfinite(ymax) and ymin < ymax):
                ymin, ymax = -5.0, 5.0
            b = np.array([xmin, xmax, ymin, ymax], dtype=float)
            logger.warning("Invalid bounds provided; using data-driven bounds [%.2f, %.2f, %.2f, %.2f]", xmin, xmax, ymin, ymax)
        else:
            b = np.array([-5.0, 5.0, -5.0, 5.0], dtype=float)
            logger.warning("Invalid bounds and no data; using default [-5,5,-5,5]")

    # Define grid sizes safely
    span_x = float(b[1] - b[0])
    span_y = float(b[3] - b[2])
    nx = max(1, int(np.floor(span_x / pixel_size)) + 1)
    ny = max(1, int(np.floor(span_y / pixel_size)) + 1)

    xi = np.linspace(b[0], b[0] + (nx - 1) * pixel_size, nx)
    yi = np.linspace(b[2], b[2] + (ny - 1) * pixel_size, ny)
    X, Y = np.meshgrid(xi, yi)

    # If no valid input points, return zeros image (and error image)
    if fx.size == 0:
        img = np.zeros_like(X, dtype=float)
        if return_error:
            err_img = np.zeros_like(X, dtype=float)
            flux_factor = pixel_size ** 2 / fiber_area
            return img * flux_factor, X, Y, err_img * flux_factor
        flux_factor = pixel_size ** 2 / fiber_area
        return img * flux_factor, X, Y

    flux_factor = pixel_size ** 2 / fiber_area

    err_img = None
    if method in ["linear", "nearest", "cubic"]:
        # griddata can fail for too few points or outside convex hull; fall back to nearest
        try:
            img = griddata((fx, fy), ff, (X, Y), method=method)
        except Exception:
            img = griddata((fx, fy), ff, (X, Y), method="nearest")
        if return_error and fiber_var is not None:
            if method == "nearest":
                # Nearest-neighbor: propagate error directly
                try:
                    nearest_idx = griddata((fx, fy), np.arange(fx.size), (X, Y), method="nearest").astype(int)
                    var_img = fiber_var[nearest_idx]
                except Exception:
                    # Fallback: simple nearest on variance values
                    var_img = griddata((fx, fy), fiber_var, (X, Y), method="nearest")
            else:
                # Approximate: interpolate variance field and take sqrt
                try:
                    var_img = griddata((fx, fy), fiber_var, (X, Y), method=method)
                except Exception:
                    var_img = griddata((fx, fy), fiber_var, (X, Y), method="nearest")
            err_img = np.sqrt(np.clip(var_img, a_min=0.0, a_max=None))

    elif method == "rbf":
        # RBF interpolation
        rbf = Rbf(fx, fy, ff, function=rbf_func)
        img = rbf(X, Y)
        if return_error and fiber_var is not None:
            rbf_var = Rbf(fx, fy, fiber_var, function=rbf_func)
            var_img = rbf_var(X, Y)
            err_img = np.sqrt(np.clip(var_img, a_min=0.0, a_max=None))

    elif method == "gdw":
        # Gaussian Distance Weighting (IDW-style)
        npts = fx.size
        k_eff = min(k, npts)
        if k_eff == 0:
            img = np.zeros_like(X, dtype=float)
            if return_error:
                err_img = np.zeros_like(X, dtype=float)
        else:
            tree = cKDTree(np.c_[fx, fy])
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
            img = np.sum(weights * ff[idx], axis=1).reshape(X.shape)
            if return_error and fiber_var is not None:
                # Variance propagation: sum(w^2 * var)
                var_val = np.sum((weights ** 2) * fiber_var[idx], axis=1).reshape(X.shape)
                err_img = np.sqrt(np.clip(var_val, a_min=0.0, a_max=None))

    else:
        raise ValueError(f"Unknown method {method}")

    # Correct for the area difference between the fiber size and the new pixel area
    img = img * flux_factor
    if return_error:
        err_img = err_img * flux_factor if err_img is not None else np.zeros_like(X, dtype=float)
        return img, X, Y, err_img
    return img, X, Y

def make_cube(wavelength, fiber_spectra, fiber_error, fiber_x, fiber_y, fiber_area,
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
        errorcube (ndarray):
            3D error cube with propagated uncertainties matching cube shape.
        x_grid (ndarray): 2D array of x-locations.
        y_grid (ndarray): 2D array of y-locations.
    """
    # Precompute output grid extent based on fiber positions (NaN-safe)
    bounds = fiber.get_fiber_bounds(fiber_x, fiber_y)
    # Validate pixel_size
    if not np.isfinite(pixel_size) or pixel_size <= 0:
        logger.warning("Invalid pixel_size=%s; defaulting to 1.0", str(pixel_size))
        pixel_size = 1.0
    # Estimate output dimensions safely
    span_x = float(bounds[1] - bounds[0])
    span_y = float(bounds[3] - bounds[2])
    nx = max(1, int(np.floor(span_x / pixel_size)) + 1)
    ny = max(1, int(np.floor(span_y / pixel_size)) + 1)

    # Diagnostics: report calculated spatial bounds, pixel size, and grid size
    try:
        m_finite = np.isfinite(fiber_x) & np.isfinite(fiber_y)
        n_finite = int(np.count_nonzero(m_finite))
    except Exception:
        n_finite = int(len(np.asarray(fiber_x).ravel()))
    logger.info(
        "Cube grid setup: bounds=[%.3f, %.3f, %.3f, %.3f], pixel_size=%.3f, method=%s, dims=(ny=%d, nx=%d), finite_fibers=%d",
        float(bounds[0]), float(bounds[1]), float(bounds[2]), float(bounds[3]),
        float(pixel_size), str(method), int(ny), int(nx), n_finite,
    )

    # Determine spectral dimensions and guard against mismatches
    n_wave = int(len(wavelength))
    n_spec = int(fiber_spectra.shape[1])
    # error spectra expected to match fiber_spectra shape
    if fiber_error is None:
        raise ValueError("fiber_error must be provided to compute error cube")
    if fiber_error.shape != fiber_spectra.shape:
        logger.warning("fiber_error shape %s != fiber_spectra shape %s; proceeding with min common length along spectral axis", fiber_error.shape, fiber_spectra.shape)
    n_spec_err = int(fiber_error.shape[1])
    n_lambda = min(n_wave, n_spec, n_spec_err)
    if not (n_wave == n_spec == n_spec_err):
        logger.warning(
            "Wavelength grid length (%d) != spectra columns (%d) and/or error columns (%d). Proceeding with %d planes.",
            n_wave, n_spec, n_spec_err, n_lambda
        )

    # Allocate cubes
    cube = np.zeros((n_lambda, ny, nx), dtype=float)
    errorcube = np.zeros((n_lambda, ny, nx), dtype=float)
    bounds = fiber.get_fiber_bounds(fiber_x, fiber_y)

    # Loop over wavelength bins
    for i, lam in enumerate(wavelength[:n_lambda]):
        flux = fiber_spectra[:, i]
        ferr = fiber_error[:, i]

        if modeling is not None and 'dar_model' in modeling:
            # Apply DAR model directly to each fiber position at this wavelength
            x_shift, y_shift = modeling['dar_model'](lam, fiber_x, fiber_y)
        else:
            x_shift, y_shift = fiber_x, fiber_y

        result = fibers_to_image(
            x_shift, y_shift, flux, fiber_area, bounds=bounds,
            pixel_size=pixel_size, method=method,
            rbf_func=rbf_func, k=k, sigma=sigma,
            fiber_error=ferr, return_error=True,
        )
        image, x_grid, y_grid, err_image = result
        cube[i] = image
        errorcube[i] = err_image

    return cube, errorcube, x_grid, y_grid