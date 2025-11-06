import logging

from astropy.stats import sigma_clipped_stats
import numpy as np
from scipy.ndimage import percentile_filter
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures

logger = logging.getLogger('antigen.trace')

def _get_peak_positions(column_profile, peak_x_positions):
    """
    Refine the vertical (y) positions of fiber peaks using a 3-point quadratic interpolation.

    Args:
        column_profile (np.ndarray): 1D array representing a collapsed spatial profile (e.g., mean over chunk).
        peak_x_positions (np.ndarray): Approximate x-positions of detected peaks in the profile.

    Returns:
        refined_peaks (np.ndarray): Subpixel-accurate y-positions of fiber peaks.
    """
    y_range = np.arange(column_profile.shape[0])
    offsets = np.array([
        peak_x_positions - 1,
        peak_x_positions,
        peak_x_positions + 1
    ], dtype=int)

    # Remove any positions that fall out of bounds
    valid = (offsets >= 0) & (offsets < column_profile.shape[0])
    valid = valid.all(axis=0)
    offsets = offsets[:, valid]

    refined_peaks = y_range[offsets[1]] - (
        (column_profile[offsets[2]] - column_profile[offsets[0]]) /
        (2. * (column_profile[offsets[2]] - 2. * column_profile[offsets[1]] + column_profile[offsets[0]]))
    )
    return refined_peaks


def _estimate_offset(nominal, detected, offset_range=(-20, 20)):
    """
    Estimate a global vertical offset between nominal and detected fiber positions.

    This function tests a range of offsets and selects the one that
    maximizes the cross-correlation between the first derivatives of the nominal and detected positions.
    In effect, this is trying to match up where the gaps between the fibers are as a robust alignment.

    Args:
        nominal (np.ndarray): 1D array of expected fiber y-positions (may contain NaN).
        detected (np.ndarray): 1D array of detected peak y-positions.
        offset_range (tuple): Tuple of (min_offset, max_offset) to search: Defaults to (-20, 20).

    Returns:
        best_offset (float): The offset (in pixels) that gives the best alignment between nominal and detected positions.
    """
    # Filter out NaN values from nominal positions
    nominal = np.asarray(nominal)
    valid_nominal_mask = np.isfinite(nominal)
    
    if not np.any(valid_nominal_mask):
        logger.warning("No valid nominal positions for offset estimation, returning 0")
        return 0.0
    
    nominal_clean = nominal[valid_nominal_mask]
    detected = np.asarray(detected)
    
    if len(detected) == 0:
        logger.warning("No detected peaks for offset estimation, returning 0")
        return 0.0
    
    # Sort the cleaned arrays
    nominal_clean = np.sort(nominal_clean)
    detected = np.sort(detected)

    # Check if we have enough points for difference calculation
    if len(nominal_clean) < 2 or len(detected) < 2:
        logger.warning("Not enough points for robust offset estimation, using simple median difference")
        # Fallback: use median difference between all nominal and detected pairs
        if len(nominal_clean) > 0 and len(detected) > 0:
            # Find closest detected peak for each nominal position
            offsets = []
            for nom_pos in nominal_clean:
                closest_detected = detected[np.argmin(np.abs(detected - nom_pos))]
                offsets.append(closest_detected - nom_pos)
            return np.median(offsets)
        else:
            return 0.0

    # Interpolate spacing (first differences) onto a common grid
    nominal_diff = np.diff(nominal_clean)
    detected_diff = np.diff(detected)

    # Use y-positions of centers between fibers as the x-axis for diffs
    nominal_x = 0.5 * (nominal_clean[1:] + nominal_clean[:-1])
    detected_x = 0.5 * (detected[1:] + detected[:-1])

    # Define a common y-grid for interpolation
    x_min = int(min(nominal_x[0], detected_x[0])) - 10
    x_max = int(max(nominal_x[-1], detected_x[-1])) + 10
    x_grid = np.arange(x_min, x_max + 1)

    # Interpolate the diffs onto the common grid
    nominal_interp = np.interp(x_grid, nominal_x, nominal_diff, left=0, right=0)
    detected_interp = np.interp(x_grid, detected_x, detected_diff, left=0, right=0)

    # Cross-correlate over range of integer offsets
    scores = []
    test_range = np.linspace(offset_range[0], offset_range[1], 100)
    for offset in test_range:
        shifted_detected = np.interp(x_grid + offset, x_grid, detected_interp)
        score = np.sum(nominal_interp * shifted_detected)
        scores.append(score)

    best_offset = test_range[np.argmax(scores)]
    return best_offset


def _match_peaks_to_nominal(nominal, detected, offset, max_sep=3.0):
    """
    Match detected peaks to nominal positions using a global offset and nearest-neighbor search.

    Each nominal position is shifted by the given offset, and then the nearest
    unmatched detected peak within max_sep is assigned to that nominal fiber.

    Args:
        nominal (np.ndarray): 1D array of expected fiber y-positions.
        detected (np.ndarray): 1D array of detected peak y-positions.
        offset (float): Global vertical offset to apply to nominal positions before matching.
        max_sep (float): Maximum allowed distance (in pixels) to consider a peak a valid match.

    Returns:
        matches (np.ndarray): Array of matched peak positions (same shape as nominal), with np.nan for unmatched fibers.
    """
    corrected_nominal = nominal + offset
    matches = np.full_like(nominal, fill_value=np.nan, dtype=np.float64)
    used = set()

    #TODO: max_sep is hardcoded
    for i, nom_y in enumerate(corrected_nominal):
        diffs = np.abs(detected - nom_y)
        candidates = [j for j in np.argsort(diffs) if diffs[j] <= max_sep and j not in used]
        if candidates:
            j = candidates[0]
            matches[i] = detected[j]
            used.add(j)
    return matches

def _estimate_peak_threshold(flat_column, peak_fraction_threshold=0.1):
    """
    Estimate the peak threshold and remove large-scale structure from a data column.

    This function takes a flat field column of data and employs percentile filtering to remove
    large-scale variations, isolating the residual structure. It then estimates the noise level
    using a robust noise estimation method based on the Mean Absolute Deviation (MAD).

    Args:
        flat_column: A numpy array or sequence representing the flat field column data.
        peak_fraction_threshold: A float defining the fraction threshold for detecting peaks
            (default is 0.01).

    Returns:
        min_peak_level: Estimated noise level in the data residuals.
        background: The computed background of the data column obtained via percentile
          filtering.
    """

    
    # Remove large-scale structure with percentile filter
    background = percentile_filter(flat_column, 5, size=51)
    residuals = flat_column - background
    
    # Robust noise estimate using MAD
    min_peak_level = np.nanmax(residuals) * peak_fraction_threshold

    return min_peak_level, background

def _evaluate_trace_chunk(flat_column, x_center, chunk_index, nominal_positions, good_fiber_indices, exclude_fiber_mask,
                          total_fibers, max_sep=3.0, peak_fraction_threshold=0.1):
    """
    Evaluates a trace chunk and determines peak positions for fiber assignments.

    This function analyzes a section of a flat column in an image, estimates noise levels, and identifies
    peaks corresponding to fiber traces. The function refines the detected peaks, matches them against
    nominal positions, and computes fiber positions for both valid and excluded fibers.

    Args:
        flat_column: 1D array representing the intensity values of a column in the flattened image.
        x_center: float specifying the x-coordinate of the column being processed.
        chunk_index: int identifying the index of the current chunk being analyzed.
        nominal_positions: 1D array representing the nominal (expected) positions of fibers.
        good_fiber_indices: 1D array of indices for fibers considered "good" for peak detection.
        exclude_fiber_mask: 1D boolean array where True indicates fibers to exclude (e.g., bad fibers).
        total_fibers: int specifying the total number of fibers in the trace.
        max_sep: float, optional; the maximum allowed separation between matched peaks and nominal positions.
        peak_fraction_threshold: float, optional; the fraction of the peak used as a threshold for peak detection.

    Returns:
        ndarray: 1D array containing the fiber positions. Fibers without valid matches will have NaN values.
    """
    logger.debug(f"Processing chunk {chunk_index+1} at x={x_center:.1f}")

    # Estimate noise level for this chunk
    min_peak_height, background = _estimate_peak_threshold(flat_column, peak_fraction_threshold=peak_fraction_threshold)

    # Subtract the background
    flat_column -= background

    # Find initial peaks using derivative method
    diff = np.diff(flat_column)
    peak_indices = np.where((diff[:-1] > 0) & (diff[1:] < 0))[0]
    peak_indices = peak_indices[peak_indices > 2]  # Avoid edge artifacts

    if len(peak_indices) == 0:
        logger.warning(f"No initial peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    peak_values = flat_column[peak_indices + 1]
    
    # Use signal-to-noise threshold instead of relative peak threshold
    strong_peaks_mask = peak_values > min_peak_height
    strong_peaks = peak_indices[strong_peaks_mask] + 1

    logger.debug(f"Chunk {chunk_index+1}: Found {len(strong_peaks)} peaks above min_height={min_peak_height:.3f}")

    if len(strong_peaks) == 0:
        logger.warning(f"No strong peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    # Refine peak positions
    refined_peaks = _get_peak_positions(flat_column, strong_peaks)
    nominal_good = nominal_positions[good_fiber_indices]

    # Estimate offset and match peaks
    offset = _estimate_offset(nominal_good, refined_peaks)
    logger.debug(f"Best fiber position offset for chunk {chunk_index+1}: {offset:.2f}")
    matched_peaks = _match_peaks_to_nominal(nominal_good, refined_peaks, offset, max_sep=max_sep)

    # Build fiber trace array
    fiber_trace = np.full(total_fibers, np.nan)
    fiber_trace[good_fiber_indices] = matched_peaks  # This preserves NaN for unmatched fibers

    # Count successful matches for diagnostics
    successful_matches = np.sum(np.isfinite(matched_peaks))
    logger.debug(f"Chunk {chunk_index+1}: Matched {successful_matches}/{len(good_fiber_indices)} good fibers")

    # Extrapolate excluded fiber positions
    for excluded_idx in np.where(exclude_fiber_mask)[0]:
        nearest_good = np.argmin(np.abs(excluded_idx - good_fiber_indices))
        ref_idx = good_fiber_indices[nearest_good]
        if np.isfinite(fiber_trace[ref_idx]):  # Check for valid (not NaN) reference
            fiber_trace[excluded_idx] = (
                fiber_trace[ref_idx] +
                nominal_positions[excluded_idx] -
                nominal_positions[ref_idx]
            )

    return fiber_trace


def get_trace(flat_image, nominal_trace_positions, exclude_fiber_mask, num_chunks=80, 
              peak_fraction_threshold=0.1, trace_poly_order=5):
    """
    Extracts fiber traces from a flat image by performing robust polynomial fits to fiber positions
    in segmented image chunks. This function processes the image to identify fiber traces,
    allows exclusion of specified fibers, and fits robust polynomials to model fiber positions.

    Args:
        flat_image (np.ndarray): The image array from which fiber traces are extracted.
        nominal_trace_positions (np.ndarray): Estimated initial positions for each fiber.
        exclude_fiber_mask (np.ndarray): Boolean mask indicating fibers to exclude from processing.
        num_chunks (int, optional): Number of chunks to divide the image columns for processing.
            Defaults to 80.
        peak_fraction_threshold (float, optional): Fractional threshold for peak identification in
            trace estimation. Defaults to 0.1.
        trace_poly_order (int, optional): Polynomial order for smoothing fiber traces. Defaults to 5.

    Returns:
        trace_array (np.ndarray): Smoothed fiber traces for all fibers over the full image.
        good_fiber_mask (np.ndarray): Mask indicating which fibers were successfully traced.
        raw_trace_matrix (np.ndarray): Raw fiber traces extracted from image chunks.
        x_chunk_centers (np.ndarray): Center positions of each image chunk.

    Raises:
        ValueError: If robust fitting fails and fallback fitting also fails.
        ImportError: If required libraries for robust fitting are unavailable.
        np.linalg.LinAlgError: If the fallback polynomial fitting encounters singular matrix issues.
    """
    total_fibers = len(nominal_trace_positions)
    good_fiber_indices = np.where(~exclude_fiber_mask)[0].astype(int)

    logger.info("Beginning fiber trace extraction")
    image = flat_image

    chunk_indices = np.array_split(np.arange(image.shape[1]), num_chunks)
    x_chunk_centers = np.array([np.mean(chunk) for chunk in chunk_indices])
    chunked_image = np.array_split(image, num_chunks, axis=1)

    # Discard edges to avoid poor flat data
    chunked_image = chunked_image[1:-1]
    x_chunk_centers = x_chunk_centers[1:-1]

    flattened_chunks = [np.median(chunk, axis=1) for chunk in chunked_image]
    raw_trace_matrix = np.full((total_fibers, len(flattened_chunks)), np.nan)  # Initialize with NaN

    center_idx = len(flattened_chunks) // 2
    center_column = flattened_chunks[center_idx]
    x_center = x_chunk_centers[center_idx]

    # Use nominal positions as-is
    center_trace = _evaluate_trace_chunk(
        center_column, x_center, center_idx, nominal_trace_positions,
        good_fiber_indices, exclude_fiber_mask, total_fibers,
        peak_fraction_threshold=peak_fraction_threshold
    )
    raw_trace_matrix[:, center_idx] = center_trace
    initial_offset = np.nanmedian(center_trace - nominal_trace_positions)
    logger.info(f"Best fiber position offset for central chunk: {initial_offset:.2f}")

    # Estimate global offset once using center

    for direction in [+1, -1]:
        for step in range(1, len(flattened_chunks)):
            k = center_idx + direction * step
            if not (0 <= k < len(flattened_chunks)):
                continue

            flat_column = flattened_chunks[k]
            x_center = x_chunk_centers[k]

            nominal_estimate = raw_trace_matrix[:, k - direction]

            trace = _evaluate_trace_chunk(
                    flat_column, x_center, k, nominal_estimate, good_fiber_indices,
                    exclude_fiber_mask, total_fibers,
                    peak_fraction_threshold=peak_fraction_threshold
                    )
            raw_trace_matrix[:, k] = trace

    logger.info(f"Smoothing fiber traces with {trace_poly_order}th-order polynomial")
    trace_array = np.zeros((total_fibers, image.shape[1]))
    x_full = np.arange(image.shape[1])
    x_max = np.max(x_full)

    good_fits = 0
    total_outliers = 0
    total_points = 0

    for fiber_index in range(total_fibers):
        valid = np.isfinite(raw_trace_matrix[fiber_index])
        if np.sum(valid) < 2 * trace_poly_order + 1:
            logger.debug(f"Fiber {fiber_index} has insufficient valid trace values "
                         f"({np.sum(valid)} points for order-{trace_poly_order} polynomial).")
            continue
        x_valid = x_chunk_centers[valid]
        y_valid = raw_trace_matrix[fiber_index, valid]

        # Robust polynomial fitting using sklearn
        poly_coeffs, inlier_mask = _fit_robust_polynomial_sklearn(
            x_valid / x_max, y_valid,
            poly_order=trace_poly_order
        )

        trace_array[fiber_index] = np.polyval(poly_coeffs, x_full / x_max)
        good_fits += 1

        # Track outlier statistics
        total_points += len(y_valid)
        total_outliers += np.sum(~inlier_mask)



    outlier_fraction = total_outliers / max(total_points, 1)
    logger.info(f"Trace extraction complete: {good_fits}/{total_fibers} successful fits, "
                f"{outlier_fraction:.1%} points identified as outliers")

    good_fiber_mask = ~exclude_fiber_mask
    return trace_array, good_fiber_mask, raw_trace_matrix, x_chunk_centers


def _fit_robust_polynomial_sklearn(x_data, y_data, poly_order=5):
    """
    Fit a polynomial to data using sklearn's robust regression methods.

    Args:
        x_data (np.ndarray): Independent variable (x-coordinates)
        y_data (np.ndarray): Dependent variable (y-coordinates)
        poly_order (int): Order of polynomial to fit

    Returns:
        coeffs (np.ndarray): Polynomial coefficients (highest order first, numpy convention)
        inlier_mask (np.ndarray): Boolean mask of inliers used in final fit

    Raises:
        ValueError: If insufficient points for the polynomial order
    """
    min_points = max(2 * poly_order + 1, 5)
    if len(x_data) < min_points:
        raise ValueError(f"Need at least {min_points} points for order-{poly_order} polynomial, got {len(x_data)}")

    regressor = RANSACRegressor()
    poly_features = PolynomialFeatures(poly_order)
    pipeline = Pipeline([
        ('poly', poly_features),
        ('robust', regressor)
    ])

    # Fit the model
    X_reshaped = x_data.reshape(-1, 1)
    pipeline.fit(X_reshaped, y_data)

    robust_regressor = pipeline.named_steps['robust']

    # Extract coefficients and convert to numpy polyfit convention (highest power first)
    fitted_model = robust_regressor.estimator_
    sklearn_coeffs = fitted_model.coef_[1:] # Remove bias column
    intercept = fitted_model.intercept_

    # sklearn gives coefficients in ascending order of powers, numpy wants descending
    poly_coeffs = np.concatenate([[intercept], sklearn_coeffs])[::-1]
    inlier_mask = robust_regressor.inlier_mask_

    return poly_coeffs, inlier_mask
