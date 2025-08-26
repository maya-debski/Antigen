import numpy as np
import logging

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
        nominal (np.ndarray): 1D array of expected fiber y-positions.
        detected (np.ndarray): 1D array of detected peak y-positions.
        offset_range (tuple): Tuple of (min_offset, max_offset) to search: Defaults to (-20, 20).

    Returns:
        best_offset (float): The offset (in pixels) that gives the best alignment between nominal and detected positions.
    """
    nominal = np.sort(np.asarray(nominal))
    detected = np.asarray(detected)

    # Interpolate spacing (first differences) onto a common grid
    nominal_diff = np.diff(nominal)
    detected_diff = np.diff(detected)

    # Use y-positions of centers between fibers as the x-axis for diffs
    nominal_x = 0.5 * (nominal[1:] + nominal[:-1])
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

    for i, nom_y in enumerate(corrected_nominal):
        diffs = np.abs(detected - nom_y)
        candidates = [j for j in np.argsort(diffs) if diffs[j] <= max_sep and j not in used]
        if candidates:
            j = candidates[0]
            matches[i] = detected[j]
            used.add(j)
    return matches

def _evaluate_trace_chunk(flat_column, x_center, chunk_index, nominal_positions, good_fiber_indices, exclude_fiber_mask,
                          total_fibers, max_sep=3.0):
    """
    Process a single image chunk: detect peaks, estimate offset, match to nominal,
    and extrapolate excluded fiber positions.

    Args:
        flat_column (np.ndarray): 1D median profile of the current chunk.
        x_center (float): X position (center) of this chunk.
        chunk_index (int): Index of the chunk.
        nominal_positions (np.ndarray): Expected positions for each fiber.
        good_fiber_indices (np.ndarray): Indices of fibers not excluded.
        exclude_fiber_mask (np.ndarray): Boolean mask for excluded fibers.
        total_fibers (int): Total number of fibers.
        max_sep (float): Maximum matching distance.

    Returns:
        fiber_trace: 1D array of matched/refined fiber positions (length = total_fibers).
    """
    logger.debug(f"Processing chunk {chunk_index+1} at x={x_center:.1f}")

    diff = np.diff(flat_column)
    peak_indices = np.where((diff[:-1] > 0) & (diff[1:] < 0))[0]
    peak_indices = peak_indices[peak_indices > 2]

    if len(peak_indices) == 0:
        logger.warning(f"No initial peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    peak_values = flat_column[peak_indices + 1]
    strong_peaks = peak_indices[peak_values > 0.3 * np.median(peak_values)] + 1

    if len(strong_peaks) == 0:
        logger.warning(f"No strong peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    refined_peaks = _get_peak_positions(flat_column, strong_peaks)
    nominal_good = nominal_positions[good_fiber_indices]

    offset = _estimate_offset(nominal_good, refined_peaks)
    matched_peaks = _match_peaks_to_nominal(nominal_good, refined_peaks, offset, max_sep=max_sep)

    fiber_trace = np.zeros(total_fibers)
    fiber_trace[good_fiber_indices] = np.where(np.isnan(matched_peaks), 0.0, matched_peaks)

    for excluded_idx in np.where(exclude_fiber_mask)[0]:
        nearest_good = np.argmin(np.abs(excluded_idx - good_fiber_indices))
        ref_idx = good_fiber_indices[nearest_good]
        if fiber_trace[ref_idx] != 0.0:
            fiber_trace[excluded_idx] = (
                fiber_trace[ref_idx] +
                nominal_positions[excluded_idx] -
                nominal_positions[ref_idx]
            )

    return fiber_trace


def get_trace(flat_image, nominal_trace_positions, exclude_fiber_mask, num_chunks=80):
    """
    Extract smoothed fiber traces from a twilight flat field image.

    Args:
        flat_image (np.ndarray): 2D array representing the flat field image.
        nominal_trace_positions (np.ndarray): 1D array of expected y-positions for each fiber.
        exclude_fiber_mask (np.ndarray): 1D boolean array marking fibers to exclude from tracing.
        num_chunks (int): Number of chunks to bin image in dispersion direction

    Returns:
        trace_array (np.ndarray): 2D array (fibers × x-pixels) of smoothed fiber traces.
        good_fiber_mask (np.ndarray): 1D boolean array indicating usable fibers.
        raw_trace_matrix (np.ndarray): 2D array (fibers × chunks) of raw trace positions.
        x_chunk_centers (np.ndarray): 1D array of x-axis chunk centers used for trace extraction.
    """
    total_fibers = len(nominal_trace_positions)
    good_fiber_indices = np.where(~exclude_fiber_mask)[0]

    logger.info("Beginning fiber trace extraction")
    image = flat_image

    chunk_indices = np.array_split(np.arange(image.shape[1]), num_chunks)
    x_chunk_centers = np.array([np.mean(chunk) for chunk in chunk_indices])
    chunked_image = np.array_split(image, num_chunks, axis=1)

    # Discard edges to avoid poor flat data
    chunked_image = chunked_image[1:-1]
    x_chunk_centers = x_chunk_centers[1:-1]

    flattened_chunks = [np.median(chunk, axis=1) for chunk in chunked_image]
    raw_trace_matrix = np.zeros((total_fibers, len(flattened_chunks)))

    center_idx = len(flattened_chunks) // 2
    center_column = flattened_chunks[center_idx]
    x_center = x_chunk_centers[center_idx]

    # Use nominal positions as-is
    center_trace = _evaluate_trace_chunk(center_column, x_center, center_idx, nominal_trace_positions,
                                         good_fiber_indices,  exclude_fiber_mask, total_fibers)
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

            trace = _evaluate_trace_chunk(flat_column, x_center, k, nominal_estimate, good_fiber_indices,
                                          exclude_fiber_mask, total_fibers)
            raw_trace_matrix[:, k] = trace

    logger.info("Smoothing fiber traces with 7th-order polynomial")
    trace_array = np.zeros((total_fibers, image.shape[1]))
    x_full = np.arange(image.shape[1])

    for fiber_index in range(total_fibers):
        valid = raw_trace_matrix[fiber_index] != 0.0
        if np.any(valid):
            poly_coeffs = np.polyfit(x_chunk_centers[valid], raw_trace_matrix[fiber_index, valid], deg=7)
            trace_array[fiber_index] = np.polyval(poly_coeffs, x_full)
        else:
            logger.debug(f"Fiber {fiber_index} has no valid trace values.")
    good_fiber_mask = ~exclude_fiber_mask
    logger.info("Trace extraction complete")
    return trace_array, good_fiber_mask, raw_trace_matrix, x_chunk_centers