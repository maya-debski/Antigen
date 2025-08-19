import numpy as np

from scipy.interpolate import griddata, Rbf
from scipy.spatial import cKDTree

def fibers_to_image(x, y, flux, grid_size=(100, 100), method="linear", rbf_func="multiquadric", k=5):
    """
    Create a synthetic image from non-uniform fiber locations and flux values.

    Args:
        x (ndarray): 1D array of fiber x-locations.
        y (ndarray): 1D array of fiber y-locations.
        flux (ndarray): 1D array of fiber flux values.
        grid_size (tuple): (nx, ny) pixels for the output image.
        method (str): interpolation method:
            - "nearest" : nearest-neighbor (griddata)
            - "linear"  : linear interpolation (griddata)
            - "cubic"   : cubic interpolation (griddata)
            - "rbf"     : radial basis function interpolation
            - "idw"     : inverse-distance weighting
        rbf_func (str): radial basis function (for method="rbf").
        k (int): number of neighbors for IDW.

    Returns:
        img (ndarray): synthetic image on uniform grid
        X, Y (ndarray): meshgrid coordinates
    """
    # Define grid
    xi = np.linspace(np.min(x), np.max(x), grid_size[0])
    yi = np.linspace(np.min(y), np.max(y), grid_size[1])
    X, Y = np.meshgrid(xi, yi)

    if method in ["linear", "nearest", "cubic"]:
        img = griddata((x, y), flux, (X, Y), method=method)

    elif method == "rbf":
        rbf = Rbf(x, y, flux, function=rbf_func)
        img = rbf(X, Y)

    elif method == "idw":
        # Inverse Distance Weighting
        tree = cKDTree(np.c_[x, y])
        dist, idx = tree.query(np.c_[X.ravel(), Y.ravel()], k=k)
        # Avoid division by zero
        dist = np.where(dist == 0, 1e-12, dist)
        weights = 1.0 / dist
        weights /= weights.sum(axis=1)[:, None]
        img = np.sum(weights * flux[idx], axis=1).reshape(X.shape)

    else:
        raise ValueError(f"Unknown method {method}")

    return img, X, Y