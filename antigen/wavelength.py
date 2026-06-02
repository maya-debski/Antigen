"""
Class-based refactor of get_wavelength_registration with clean, fail-fast logic.

This module provides a class WavelengthCalibrator that re-implements the
wavelength calibration pipeline from scratch (informed by prior behavior),
without importing helpers from antigen.wavelength. It preserves the legacy
public API via a top-level get_wavelength_registration wrapper with the same
signature and return values.

Design goals
- Deterministic flow with explicit validations; avoid try/except masking.
- Clear separation of responsibilities across small private methods.
- Only rely on shared domain utilities outside this module where appropriate
  (e.g., sky continuum identification from antigen.sky, peak refinement from
  antigen.trace), not legacy wavelength helpers.

Returned tuple
(wavelength_solution, residuals, fitted_trace_positions, arc_pixel_locations)
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from astropy.stats import mad_std, biweight_midvariance
from astropy.convolution import Gaussian1DKernel, convolve

from antigen.sky import identify_sky_pixels
from antigen.trace import _get_peak_positions
from antigen.registration import (
    match_width as reg_match_width,
    build_registration_image as reg_build_registration_image,
    estimate_subpixel_shift_1d_xcorr as reg_estimate_subpixel_shift,
    fit_poly2d_least_squares as reg_fit_poly2d,
)
import logging

logger = logging.getLogger(__name__)


Array = np.ndarray


@dataclass(slots=True)
class _PCCResult:
    # Shift surface fit: dx that maps SOURCE (rectified clip) -> REFERENCE (registration image)
    eval_src_to_ref: callable
    # Convenience inverse mapping (approximate): REFERENCE -> SOURCE
    eval_ref_to_src: callable
    # Optional dense map along all fibers/columns for ref->all mapping
    dx_map_ref_to_all: Array | None
    # Origin row of the rectified clip within detector Y coordinates
    y_start: int


class WavelengthCalibrator:
    def __init__(
        self,
        spectra: Array,
        trace_positions: Array,
        valid_fibers: Array,
        arc_pixel_guesses: Array,
        arc_wavelengths: Array,
        *,
        peak_threshold: float = 5.0,
        reference_fiber_index: int = 130,
        good_arc_residual_limit: float = 0.2,
        binned: bool = False,
        out_folder: str | Path | None = None,
        detector_dimensions: dict | None = None,
        enable_plots: bool = True,
        # Optional overrides (also supported via detector_dimensions['wavelength'] keys)
        detection_window: int | None = None,
        poly_order: int | None = None,
        pcc_tile_h: int | None = None,
        pcc_tile_w: int | None = None,
        pcc_step_y: int | None = None,
        pcc_step_x: int | None = None,
        pcc_poly_order: int | None = None,
        pcc_max_corr_lag: int | None = None,
    ) -> None:
        # Basic shape validations (fail fast)
        if spectra.ndim != 2:
            raise ValueError("spectra must be 2D [num_fibers, num_pixels]")
        if trace_positions.ndim != 2:
            raise ValueError("trace_positions must be 2D [num_fibers, num_pixels]")
        if spectra.shape != trace_positions.shape:
            raise ValueError("spectra and trace_positions must have the same shape")
        if arc_pixel_guesses.ndim != 1:
            raise ValueError("arc_pixel_guesses must be 1D [num_lines]")
        if arc_wavelengths.ndim != 1:
            raise ValueError("arc_wavelengths must be 1D [num_lines]")
        if arc_pixel_guesses.shape[0] != arc_wavelengths.shape[0]:
            raise ValueError("arc_pixel_guesses and arc_wavelengths must have same length")
        if valid_fibers.shape[0] != spectra.shape[0]:
            raise ValueError("valid_fibers length must equal number of fibers")

        self.spectra = np.asarray(spectra, dtype=float)
        self.trace_positions = np.asarray(trace_positions, dtype=float)
        self.valid_fibers = np.asarray(valid_fibers, dtype=bool)
        self.arc_pixel_guesses = np.asarray(arc_pixel_guesses, dtype=float)
        self.arc_wavelengths = np.asarray(arc_wavelengths, dtype=float)
        self.peak_threshold = float(peak_threshold)
        self.reference_fiber_index = int(reference_fiber_index)
        self.good_arc_residual_limit = float(good_arc_residual_limit)
        self.binned = bool(binned)
        self.out_folder = Path(out_folder) if out_folder is not None else None
        self.detector_dimensions = detector_dimensions or {}
        self.enable_plots = bool(enable_plots)
        # Thresholds and model orders (configurable)
        self.sn_thresh = float(peak_threshold)
        self.poly_order = 4
        # Detection window (pixels) for matching guesses to peaks
        self.detection_window = 7
        # PCC tiling/fit parameters
        self.pcc_tile_h = 256
        self.pcc_tile_w = 256
        self.pcc_step_y = 128
        self.pcc_step_x = 128
        self.pcc_poly_order = 3
        # Default PCC max correlation lag: solid default per spectrograph is 24 px
        self.pcc_max_corr_lag = 24
        # Apply optional overrides from config first, then constructor params override
        self._apply_overrides_from_config()
        if detection_window is not None:
            self.detection_window = int(detection_window)
        if poly_order is not None:
            self.poly_order = int(poly_order)
        if pcc_tile_h is not None:
            self.pcc_tile_h = int(pcc_tile_h)
        if pcc_tile_w is not None:
            self.pcc_tile_w = int(pcc_tile_w)
        if pcc_step_y is not None:
            self.pcc_step_y = int(pcc_step_y)
        if pcc_step_x is not None:
            self.pcc_step_x = int(pcc_step_x)
        if pcc_poly_order is not None:
            self.pcc_poly_order = int(pcc_poly_order)
        if pcc_max_corr_lag is not None:
            self.pcc_max_corr_lag = int(pcc_max_corr_lag)
        # Derived
        self.num_fibers, self.num_pixels = self.spectra.shape
        if not (0 <= self.reference_fiber_index < self.num_fibers):
            raise ValueError("reference_fiber_index out of range")

    def _apply_overrides_from_config(self) -> None:
        """Apply optional parameter overrides from detector_dimensions['wavelength'].

        Supported keys:
        - detection_window, poly_order, pcc_tile_h, pcc_tile_w, pcc_step_y, pcc_step_x,
          pcc_poly_order, pcc_max_corr_lag
        """
        cfg = self.detector_dimensions or {}
        wl = cfg.get('wavelength') if isinstance(cfg, dict) else None
        if not isinstance(wl, dict):
            return
        self.detection_window = int(wl.get('detection_window', self.detection_window))
        self.poly_order = int(wl.get('poly_order', self.poly_order))
        self.pcc_tile_h = int(wl.get('pcc_tile_h', self.pcc_tile_h))
        self.pcc_tile_w = int(wl.get('pcc_tile_w', self.pcc_tile_w))
        self.pcc_step_y = int(wl.get('pcc_step_y', self.pcc_step_y))
        self.pcc_step_x = int(wl.get('pcc_step_x', self.pcc_step_x))
        self.pcc_poly_order = int(wl.get('pcc_poly_order', self.pcc_poly_order))
        self.pcc_max_corr_lag = wl.get('pcc_max_corr_lag', self.pcc_max_corr_lag)
        if self.pcc_max_corr_lag is not None:
            self.pcc_max_corr_lag = int(self.pcc_max_corr_lag)

    # ------------------------- Public API -------------------------
    def run(self) -> Tuple[Array, Array, Array, Array]:
        """Execute the wavelength calibration pipeline and return outputs.

        Returns (wavelength_solution, residuals, fitted_trace_positions, arc_pixel_locations).
        """
        # 1) Adjust for binning (no-op here but validates)
        self._adjust_for_binning()

        # 2) Select good fibers and build normalized array used for rectification
        good_mask, norm = self._select_good_fibers()

        # 2b) Fast rectification on good fibers to build source image for PCC
        rect_clip, y_start = self._fast_rectify_per_column(norm[good_mask], self.trace_positions[good_mask])
        # 3) Build reference registration image from raw spectra[ref]
        ref_profile = self._build_ref_profile(self.spectra[self.reference_fiber_index])
        # Match width to rectified clip
        if rect_clip.shape[1] != ref_profile.shape[0]:
            ref_profile = reg_match_width(ref_profile, rect_clip.shape[1])
        reg_img = reg_build_registration_image(ref_profile, rect_clip.shape[0])

        # 3: PCC diagnostics and shift surface
        pcc_res = self._run_pcc_diagnostics(rect_clip, reg_img, y_start, good_mask)

        # 3a: Detect reference arc lines on the ref profile (1D)
        ref_detected = self._detect_ref_arclines(ref_profile, self.arc_pixel_guesses, window=self.detection_window)

        # 3b: Build base ref positions and propagate to all fibers using dx map
        arc_pixel_locations = self._build_arc_positions_from_dx_map(ref_detected, pcc_res, good_mask)
        # Keep a copy of pre-refinement arc positions for diagnostics
        pre_refined_arc_pixel_locations = np.array(arc_pixel_locations, dtype=float, copy=True)

        # 3c: Refine arc positions per good fiber by re-detecting peaks near current guesses
        arc_pixel_locations = self._refine_arc_positions_on_good_fibers(arc_pixel_locations, good_mask, window=self.detection_window)

        # Optional diagnostics: quick plot (two x-limits) — not using try/except
        if self.enable_plots and self.out_folder is not None:
            self._plot_norm_spectra_with_arc_positions(norm, good_mask, arc_pixel_locations)

        # 4) Fit across fibers (trace-space), then final wavelength solution
        fitted_trace_positions, residuals = self._fit_arc_line_positions(arc_pixel_locations, pre_refined_arc_pixels=pre_refined_arc_pixel_locations)

        wavelength_solution = self._fit_wavelength_polynomials(
            fitted_trace_positions, arc_pixel_locations
        )


        return wavelength_solution, residuals, fitted_trace_positions, arc_pixel_locations

    # ------------------------- Step helpers -------------------------
    def _adjust_for_binning(self) -> None:
        # Adjust arc line pixel guesses for 2x binned data
        if self.binned not in (True, False):
            raise ValueError("binned must be boolean")
        if self.binned:
            # Halve the arc pixel guesses so they are in binned-pixel units
            self.arc_pixel_guesses = self.arc_pixel_guesses / 2.0

    def _select_good_fibers(self) -> Tuple[Array, Array]:
        """Identify a robust set of fibers for registration/fitting and return normalized spectra.

        Applies a light continuum subtraction per fiber before normalization to improve
        registration (PCC) robustness, then scales by the per-fiber 98th percentile of
        the absolute continuum-subtracted profile. Selects fibers that are finite across
        pixels and flagged valid. Raises if too few remain for stable PCC.
        """
        # Continuum subtraction per fiber (fail-safe)
        cont_sub = np.empty_like(self.spectra, dtype=float)
        n_fail = 0
        for f in range(self.num_fibers):
            spec = self.spectra[f]
            try:
                mask, cont = identify_sky_pixels(spec)
                prof = spec - cont
            except Exception:
                # Fallback: no continuum subtraction for this fiber
                prof = spec.astype(float, copy=True)
                n_fail += 1
            cont_sub[f] = prof
        if n_fail > 0:
            logger.debug("Continuum subtraction failed on %d/%d fibers in _select_good_fibers; using raw spectra for those.", n_fail, self.num_fibers)

        # Robust per-fiber scale using 98th percentile of absolute values
        perc98 = np.nanpercentile(np.abs(cont_sub), 98, axis=1)
        scale = np.where((perc98 > 0) & np.isfinite(perc98), perc98, 1.0)
        norm = (cont_sub.T / scale).T

        # Basic validity: all finite across pixels and marked valid
        finite_rows = np.all(np.isfinite(norm), axis=1)
        good_mask = self.valid_fibers & finite_rows
        if np.count_nonzero(good_mask) < max(5, self.num_fibers // 20):
            # Require at least some fibers for robust PCC
            raise ValueError("Too few good fibers after normalization to proceed")
        return good_mask, norm

    def _fast_rectify_per_column(self, norm_good: Array, y_locs_good: Array) -> Tuple[Array, int]:
        """Fast per-column interpolation from scattered Y locations (per fiber) to an integer Y grid.

        Mirrors the legacy wavelength._fast_rectify_per_column behavior for image quality and
        non-finite handling:
        - Determines output grid size from detector_dimensions if provided, otherwise from
          global trace_positions/spectra.
        - Per column: averages duplicate Y samples, then interpolates onto the full integer Y grid
          with NaNs outside the convex hull of valid samples.
        - Crops to rows fully covered (no NaNs across columns) to avoid NaN borders and returns
          the starting Y offset (y_start) of the cropped clip in detector coordinates.
        """
        # Determine output grid size (nrows, ncols)
        if not self.detector_dimensions or ('X' not in self.detector_dimensions or 'Y' not in self.detector_dimensions):
            ncols = int(self.trace_positions.shape[1])
            if np.isfinite(self.trace_positions).any():
                nrows = int(np.nanmax(self.trace_positions) + 1)
            else:
                nrows = int(self.spectra.shape[0])
        else:
            ncols = int(self.detector_dimensions['X'])
            nrows = int(self.detector_dimensions['Y'])

        y_grid = np.arange(nrows, dtype=float)
        rectified = np.full((nrows, ncols), np.nan, dtype=float)

        # Per-column interpolation from scattered fiber Y to integer grid
        max_x = min(ncols, y_locs_good.shape[1])
        for x in range(max_x):
            y = y_locs_good[:, x]
            v = norm_good[:, x]
            m = np.isfinite(y) & np.isfinite(v)
            if np.count_nonzero(m) == 0:
                continue
            ys = y[m].astype(float)
            vs = v[m].astype(float)
            if ys.size == 1:
                yi = int(np.clip(int(np.round(ys[0])), 0, nrows - 1))
                rectified[yi, x] = float(vs[0])
                continue
            order = np.argsort(ys)
            ys = ys[order]
            vs = vs[order]
            u, inv = np.unique(ys, return_inverse=True)
            sums = np.bincount(inv, weights=vs)
            counts = np.bincount(inv)
            vs_avg = sums / np.maximum(counts, 1)
            rectified[:, x] = np.interp(y_grid, u, vs_avg, left=np.nan, right=np.nan)

        # Clip to rows fully covered (avoid NaN borders)
        finite_counts = np.sum(np.isfinite(rectified), axis=1)
        valid_rows = np.where(finite_counts == ncols)[0]
        if valid_rows.size >= 2:
            y_start = int(valid_rows[0])
            y_end = int(valid_rows[-1] + 1)
            rect_clip_local = rectified[y_start:y_end, :]
        else:
            y_start = 0
            rect_clip_local = rectified
        return rect_clip_local, y_start

    def _build_ref_profile(self, spectrum: Array) -> Array:
        # Continuum subtraction using sky pixel identification
        mask, cont = identify_sky_pixels(spectrum)
        prof = spectrum - cont
        # Robust normalization
        p98 = np.nanpercentile(np.abs(prof), 98)
        if not np.isfinite(p98) or p98 <= 0:
            raise ValueError("Invalid reference profile scale")
        return prof / p98


    def _run_pcc_diagnostics(self, rect_img: Array, reg_img: Array, y_start: int, good_mask: Array) -> _PCCResult:
        if rect_img.shape != reg_img.shape:
            raise ValueError("rectified and registration images must have identical shapes")
        H, W = rect_img.shape
        # Grid-sample tiles and estimate dx via subpixel phase cross-correlation
        tile_h, tile_w, step_y, step_x = self.pcc_tile_h, self.pcc_tile_w, self.pcc_step_y, self.pcc_step_x
        logger.debug("Running PCC diagnostics on %d x %d rectified image", H, W)

        def _sample_tiles(max_lag_val: int):
            coords_y, coords_x, dx_meas = [], [], []
            used_tiles = 0
            skipped_tiles = 0
            boundary_hits = 0
            for y0 in range(0, H - tile_h + 1, step_y):
                for x0 in range(0, W - tile_w + 1, step_x):
                    src = rect_img[y0:y0 + tile_h, x0:x0 + tile_w]
                    ref = reg_img[y0:y0 + tile_h, x0:x0 + tile_w]
                    # Optional diagnostic: when the chosen lag is near the boundary, plot c vs lag
                    diag_cb = None
                    if self.enable_plots and self.out_folder is not None:
                        def diag_cb(lags, c, info, y0=y0, x0=x0, tile_h=tile_h, tile_w=tile_w, max_lag_val=max_lag_val):
                            try:
                                import matplotlib
                                matplotlib.use("Agg", force=True)
                                import matplotlib.pyplot as plt
                                from pathlib import Path
                                outdir = Path(self.out_folder) / "pcc_xcorr"
                                outdir.mkdir(parents=True, exist_ok=True)
                                fig, ax = plt.subplots(figsize=(5, 3))
                                ax.plot(lags, c, marker='o', ms=3, lw=0.8)
                                k = int(info.get('k', 0)); lag0 = info.get('lag0', 0.0); mlag = info.get('max_lag', None)
                                if 0 <= k < len(c):
                                    ax.plot([lags[k]], [c[k]], 'ro', ms=5, label=f"k={k}, lag={lag0:.1f}")
                                if mlag is not None:
                                    ax.axvline(+int(mlag), color='orange', ls='--', lw=0.8)
                                    ax.axvline(-int(mlag), color='orange', ls='--', lw=0.8, label='±max_lag')
                                ax.set_xlabel('Lag (px)')
                                ax.set_ylabel('Cross-correlation')
                                ax.set_title(f"PCC xcorr at tile center=({y0 + tile_h/2:.0f},{x0 + tile_w/2:.0f})")
                                ax.grid(alpha=0.2)
                                ax.legend(loc='best', fontsize=8)
                                fname = outdir / f"xcorr_y{y0}_x{x0}_lag{int(max_lag_val)}.png"
                                fig.tight_layout(); fig.savefig(fname, dpi=140); plt.close(fig)
                                # Also save the raw curve for inspection
                                import numpy as _np  # local alias to avoid top-level deps
                                _np.savetxt(outdir / f"xcorr_y{y0}_x{x0}_lag{int(max_lag_val)}.csv", _np.column_stack([lags, c]), delimiter=",", header="lag,c", comments="")
                                logger.debug("Saved PCC xcorr diagnostic: %s", fname)
                            except Exception as e:
                                logger.debug("Failed to save PCC xcorr diagnostic at y0=%d x0=%d: %s", y0, x0, e)
                    dy, dx = reg_estimate_subpixel_shift(ref, src, max_lag=int(max_lag_val), diag_callback=diag_cb)  # shift to map src -> ref
                    logger.debug("PCC tile at (y0=%d, x0=%d) center=(%.1f, %.1f): dy=%.3f, dx=%.3f", y0, x0, y0 + tile_h / 2.0, x0 + tile_w / 2.0, dy, dx)
                    if abs(dy) > tile_h * 0.5:
                        skipped_tiles += 1
                        logger.debug("Skipping tile at (y0=%d, x0=%d): |dy|=%.3f exceeds %.1f (tile_h/2)", y0, x0, abs(dy), tile_h * 0.5)
                        continue
                    coords_y.append(y0 + tile_h / 2.0)
                    coords_x.append(x0 + tile_w / 2.0)
                    dx_meas.append(dx)
                    used_tiles += 1
                    # Flag boundary hits when dx is at the correlation search edge
                    if int(max_lag_val) > 0 and abs(dx) >= (int(max_lag_val) - 0.75):
                        boundary_hits += 1
            return coords_y, coords_x, dx_meas, used_tiles, skipped_tiles, boundary_hits

        # Initial lag: explicit from config if provided; otherwise solid default 24 px
        initial_lag = int(self.pcc_max_corr_lag) if self.pcc_max_corr_lag is not None else 24
        coords_y, coords_x, dx_meas, used_tiles, skipped_tiles, boundary_hits = _sample_tiles(initial_lag)
        total_tiles = used_tiles + skipped_tiles
        boundary_frac = (boundary_hits / max(used_tiles, 1)) if used_tiles else 0.0
        logger.debug("PCC sampling summary (lag=%d): used %d, skipped %d, boundary_hits=%d (%.1f%%)",
                     initial_lag, used_tiles, skipped_tiles, boundary_hits, 100.0 * boundary_frac)

        # If many tiles are saturating the lag boundary, expand the lag once and re-sample
        final_lag = initial_lag
        if used_tiles >= 6 and boundary_frac > 0.2:
            # Expand but cap to avoid pathological values
            cap = max(24, W // 2)
            proposed = max(initial_lag + 8, initial_lag * 2)
            new_lag = int(min(cap, proposed))
            if new_lag > final_lag:
                logger.info("PCC: expanding max_lag from %d to %d due to %.0f%% boundary hits", initial_lag, new_lag, 100.0 * boundary_frac)
                coords_y, coords_x, dx_meas, used_tiles, skipped_tiles, boundary_hits = _sample_tiles(new_lag)
                final_lag = new_lag
                boundary_frac = (boundary_hits / max(used_tiles, 1)) if used_tiles else 0.0
                logger.debug("PCC sampling summary (lag=%d): used %d, skipped %d, boundary_hits=%d (%.1f%%)",
                             final_lag, used_tiles, skipped_tiles, boundary_hits, 100.0 * boundary_frac)

        if len(dx_meas) < 6:
            raise ValueError("Insufficient PCC samples to fit shift surface")
        x = np.asarray(coords_x)
        y = np.asarray(coords_y)
        z = np.asarray(dx_meas)
        eval_src_to_ref = reg_fit_poly2d(x, y, z, order=self.pcc_poly_order)

        def eval_ref_to_src(X, Y):
            return -eval_src_to_ref(X, Y)

        # Build dx map from reference to all fibers along their trace y(x)
        # We compute E(y,x) for every fiber (using its y(x)) and subtract the ref's shift
        # so that dx_ref_to_fiber(x) = E_ref(x) - E_fiber(x) maps ref x to fiber x.
        y_all = self.trace_positions
        # Use only valid finite rows; others will be NaN in the map
        E_all = np.full_like(y_all, np.nan, dtype=float)
        X = np.tile(np.arange(W), (self.num_fibers, 1))
        m = np.isfinite(y_all)
        E_all[m] = eval_src_to_ref(X[m], y_all[m] - y_start)
        E_ref = E_all[self.reference_fiber_index]
        dx_map_ref_to_all = (E_ref[None, :] - E_all)
        return _PCCResult(eval_src_to_ref=eval_src_to_ref, eval_ref_to_src=eval_ref_to_src, dx_map_ref_to_all=dx_map_ref_to_all, y_start=y_start)


    def _detect_ref_arclines(self, ref_profile: Array, guesses: Array, window: int = 7) -> Array:
        # Smooth with a modest Gaussian to suppress noise
        kernel = Gaussian1DKernel(1.0)
        sm = convolve(ref_profile, kernel)

        noise = np.sqrt(biweight_midvariance(sm, ignore_nan=True))
        # Simple peak detection via derivative sign change
        d = sm[1:] - sm[:-1]
        loc = np.where((d[:-1] > 0) & (d[1:] < 0))[0] + 1
        if loc.size == 0:
            return np.full_like(guesses, np.nan, dtype=float)
        # Sub-pixel refinement
        signal = sm[loc]
        loc = _get_peak_positions(sm, loc)

        # Match nearest to each guess; enforce proximity window
        out = np.full_like(guesses, np.nan, dtype=float)
        for j, g in enumerate(np.asarray(guesses, dtype=float)):
            if not np.isfinite(g):
                continue
            k = int(np.argmin(np.abs(loc - g)))
            if signal[k] / noise < self.sn_thresh:
                logger.debug("Rejecting ref arc line %d: signal/noise (%.2f < %.2f)",
                            j, signal[k] / noise, self.sn_thresh)
                continue
            if np.abs(loc[k] - g) <= max(2.0, window):
                out[j] = loc[k]
        return out

    def _build_arc_positions_from_dx_map(self, base_ref: Array, pcc_res: _PCCResult, good_mask: Array) -> Array:
        """Propagate reference-fiber arc detections to all fibers using the PCC dx map.

        Only lines detected in the reference are propagated; non-good fibers are masked to NaN
        to exclude them from later fitting stages. Results are clipped to detector bounds.
        """
        n_lines = base_ref.shape[0]
        out = np.full((self.num_fibers, n_lines), np.nan, dtype=float)
        dx_map = pcc_res.dx_map_ref_to_all
        if dx_map is None:
            # No propagation possible; only detected lines in ref are left as NaN (cannot tile safely)
            return out
        # Only propagate lines detected in the reference
        detected = np.isfinite(base_ref)
        cols = np.rint(np.clip(base_ref[detected], 0, self.num_pixels - 1)).astype(int)
        for j_idx, j in enumerate(np.where(detected)[0]):
            xj = cols[j_idx]
            out[:, j] = base_ref[j] + dx_map[:, xj]
        # Mask non-good fibers so they never enter fits
        out[~good_mask, :] = np.nan
        # Clip to bounds
        np.clip(out, 0, self.num_pixels - 1, out=out)
        return out

    def _refine_arc_positions_on_good_fibers(self, arc_pos: Array, good_mask: Array, window: int = 7) -> Array:
        if arc_pos.shape[0] != self.num_fibers:
            raise ValueError("arc_pos has incompatible number of fibers")
        updated = np.array(arc_pos, dtype=float, copy=True)
        num_lines = arc_pos.shape[1]
        # Process only good fibers
        idx_good = np.where(good_mask & self.valid_fibers)[0]
        for f in idx_good:
            spectrum = self.spectra[f]
            try:
                profile = self._build_ref_profile(spectrum)
            except Exception:
                # If profile building fails, skip this fiber without altering values
                continue
            guesses = updated[f]
            # If all guesses are NaN, skip
            if not np.any(np.isfinite(guesses)):
                continue
            detections = self._detect_ref_arclines(profile, guesses, window=window)
            # Replace only where detection is finite; keep NaNs where not found
            m = np.isfinite(detections)
            if np.any(m):
                updated[f, m] = detections[m]
        # Clip bounds to detector width
        np.clip(updated, 0, self.num_pixels - 1, out=updated)
        return updated

    def _plot_norm_spectra_with_arc_positions(self, norm: Array, good_mask: Array, arc_pos: Array) -> None:
        """Plot small windows (~50A) around each arc line for the reference fiber.

        Replaces the previous multi-fiber overview. For each line in the arc list, we:
        - Determine the expected pixel on the reference fiber (from arc_pos[ref] or fall back to guesses).
        - Estimate dispersion (A/pixel) from available (pixel, wavelength) pairs to size the window (~50A total).
        - Plot the normalized spectrum segment around that pixel and draw a vertical line at the expected pixel.
        """
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        # Choose the reference fiber
        ref = int(self.reference_fiber_index)
        if not (0 <= ref < norm.shape[0]):
            raise ValueError("reference_fiber_index out of range for plotting")
        n_pix = norm.shape[1]

        # Build expected pixel locations for ref fiber
        expected_px = None
        if arc_pos is not None and arc_pos.ndim == 2 and arc_pos.shape[0] > ref:
            expected_px = np.array(arc_pos[ref], dtype=float)
        if expected_px is None or expected_px.shape[0] != self.arc_wavelengths.shape[0]:
            # Fallback to guesses if arc_pos missing/mismatched
            expected_px = np.array(self.arc_pixel_guesses, dtype=float)

        # Estimate dispersion (A/pixel) using linear fit wl ~ a + b * pixel
        wl = np.array(self.arc_wavelengths, dtype=float)
        m = np.isfinite(wl) & np.isfinite(expected_px)
        slope_A_per_pix = np.nan
        if np.count_nonzero(m) >= 2:
            try:
                b, a = np.polyfit(expected_px[m], wl[m], 1)  # returns [b, a] for y = b*x + a
                slope_A_per_pix = float(b)
            except Exception:
                slope_A_per_pix = np.nan
        # Convert ~50A total window to pixels
        total_A = 50.0
        if np.isfinite(slope_A_per_pix) and abs(slope_A_per_pix) > 1e-6:
            half_width_pix = 0.5 * total_A / abs(slope_A_per_pix)
            half_width_pix = float(np.clip(half_width_pix, 8.0, 200.0))
        else:
            # Safe fallback half-width in pixels if dispersion unknown
            half_width_pix = 20.0

        # Prepare plotting grid
        n_lines = int(wl.size)
        # Only include lines with finite expected pixel inside detector bounds
        valid_lines = [j for j in range(n_lines) if np.isfinite(expected_px[j]) and 0 <= expected_px[j] < n_pix]
        if not valid_lines:
            # Nothing to plot
            return
        cols = 4
        rows = int(np.ceil(len(valid_lines) / cols))
        fig_w = 3.0 * cols
        fig_h = 2.2 * rows
        fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h), squeeze=False, sharex=False, sharey=False)

        x = np.arange(n_pix)
        for k, j in enumerate(valid_lines):
            r = k // cols
            c = k % cols
            ax = axes[r, c]
            cx = float(expected_px[j])
            lo = int(max(0, np.floor(cx - half_width_pix)))
            hi = int(min(n_pix - 1, np.ceil(cx + half_width_pix)))
            if hi <= lo:
                hi = min(n_pix - 1, lo + 1)
            seg_x = x[lo:hi + 1]
            seg_y = norm[ref, lo:hi + 1]
            ax.plot(seg_x, seg_y, lw=0.7)
            ax.axvline(cx, color='r', alpha=0.6, linestyle='--', label='expected px')
            lbl_wl = float(wl[j]) if np.isfinite(wl[j]) else np.nan
            if np.isfinite(lbl_wl):
                ax.set_title(f"λ={lbl_wl:.1f}Å @ px={cx:.2f}")
            else:
                ax.set_title(f"line {j} @ px={cx:.2f}")
            ax.grid(alpha=0.2)
        # Hide any unused axes
        for idx in range(len(valid_lines), rows * cols):
            r = idx // cols
            c = idx % cols
            axes[r, c].axis('off')

        fig.tight_layout()
        out = Path(self.out_folder) / "ref_fiber_arc_windows.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)

    def _fit_arc_line_positions(self, arc_pixels: Array, pre_refined_arc_pixels: Array | None = None) -> Tuple[Array, Array]:
        num_fibers, num_lines = arc_pixels.shape
        fitted_positions = np.full_like(arc_pixels, np.nan, dtype=float)
        residuals = np.full(num_lines, np.nan, dtype=float)

        # Precompute per-fiber trace-position coordinate at every detector column
        trace_pos_map = self.trace_positions  # shape: (num_fibers, num_pixels)

        # Collect per-line data for optional diagnostics (fit scatter and residuals)
        diag_data = []  # list of dicts per fitted line

        for line_idx in range(num_lines):
            measured_pixels_for_line = arc_pixels[:, line_idx]
            # Valid if measurement exists and fiber is flagged as valid globally
            valid_mask = np.isfinite(measured_pixels_for_line) * self.valid_fibers
            if valid_mask.sum() > self.poly_order + 1:
                # Use the average measured pixel to select a representative column
                approx_column = int(np.nanmean(measured_pixels_for_line))
                approx_column = int(np.clip(approx_column, 0, trace_pos_map.shape[1] - 1))
                # Trace-position coordinate for all fibers at that column
                trace_coords_all_fibers = trace_pos_map[:, approx_column]
                # Trace-position coordinates for fibers where this line was measured (post-refine)
                measured_cols_int = np.array(measured_pixels_for_line[valid_mask], dtype=int)
                measured_cols_int = np.clip(measured_cols_int, 0, trace_pos_map.shape[1] - 1)
                trace_coords_valid_fibers = trace_pos_map[valid_mask, measured_cols_int]

                # Optional: pre-refinement measurements for overlay
                y_valid_pre = None
                x_meas_pre = None
                if pre_refined_arc_pixels is not None and pre_refined_arc_pixels.shape == arc_pixels.shape:
                    pre_vec = pre_refined_arc_pixels[:, line_idx]
                    pre_mask = np.isfinite(pre_vec) * self.valid_fibers
                    if np.count_nonzero(pre_mask) >= 2:
                        pre_cols_int = np.array(pre_vec[pre_mask], dtype=int)
                        pre_cols_int = np.clip(pre_cols_int, 0, trace_pos_map.shape[1] - 1)
                        y_valid_pre = trace_pos_map[pre_mask, pre_cols_int].astype(float)
                        x_meas_pre = pre_vec[pre_mask].astype(float)

                # Normalize coordinates by detector width for numerical stability
                detector_width = float(self.detector_dimensions.get('X', self.num_pixels))
                x_norm = measured_pixels_for_line[valid_mask] / detector_width
                y_norm = trace_coords_valid_fibers / detector_width

                # Fit polynomial: pixel_x ~ P(trace_coord)
                poly_coef = np.polyfit(y_norm, x_norm, self.poly_order)

                # Evaluate for all fibers at the representative column, then un-normalize
                fitted_positions[:, line_idx] = (
                    np.polyval(poly_coef, trace_coords_all_fibers / detector_width) * detector_width
                )

                # Robust residual scatter for the fibers used in the fit
                fit_residuals = measured_pixels_for_line[valid_mask] - fitted_positions[valid_mask, line_idx]
                residuals[line_idx] = mad_std(fit_residuals)

                # Save diagnostic info for this line
                try:
                    wl_val = float(self.arc_wavelengths[line_idx]) if line_idx < self.arc_wavelengths.size and np.isfinite(self.arc_wavelengths[line_idx]) else np.nan
                except Exception:
                    wl_val = np.nan
                # Dense curve for plotting
                y_min = np.nanmin(trace_coords_all_fibers)
                y_max = np.nanmax(trace_coords_all_fibers)
                if not np.isfinite(y_min) or not np.isfinite(y_max) or y_max <= y_min:
                    y_min, y_max = 0.0, float(self.num_fibers - 1)
                y_dense = np.linspace(y_min, y_max, 200)
                x_fit_dense = np.polyval(poly_coef, y_dense / detector_width) * detector_width
                diag_entry = {
                    'line_idx': line_idx,
                    'wavelength': wl_val,
                    'y_valid': trace_coords_valid_fibers.astype(float),
                    'x_meas': measured_pixels_for_line[valid_mask].astype(float),
                    'y_dense': y_dense.astype(float),
                    'x_fit_dense': x_fit_dense.astype(float),
                    'mad_resid': float(residuals[line_idx]),
                    'poly_coef': np.array(poly_coef, dtype=float),
                    'detector_width': float(detector_width),
                }
                if y_valid_pre is not None and x_meas_pre is not None:
                    diag_entry['y_valid_pre'] = y_valid_pre
                    diag_entry['x_meas_pre'] = x_meas_pre
                diag_data.append(diag_entry)

        # Optional: produce a multi-page PDF with one page per line showing fit and residuals
        if self.enable_plots and self.out_folder is not None and len(diag_data) > 0:
            try:
                import matplotlib
                matplotlib.use("Agg", force=True)
                import matplotlib.pyplot as plt
                from matplotlib.backends.backend_pdf import PdfPages

                pdf_path = Path(self.out_folder) / "xcoord_fit_per_line.pdf"
                with PdfPages(pdf_path) as pdf:
                    for item in diag_data:
                        fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
                        ax0, ax1 = axes
                        # Left: measured vs trace coord with fitted curve
                        ax0.plot(item['y_dense'], item['x_fit_dense'], 'r-', lw=1.2, label='poly fit (post)')
                        ax0.plot(item['y_valid'], item['x_meas'], 'o', ms=3, alpha=0.9, label='measured (post)')
                        # Overlay pre-refine points if available
                        if 'y_valid_pre' in item and 'x_meas_pre' in item:
                            ax0.plot(item['y_valid_pre'], item['x_meas_pre'], 'x', ms=4, alpha=0.8, color='tab:green', label='measured (pre)')
                        ax0.set_xlabel('Trace coord (Y)')
                        ax0.set_ylabel('Pixel x')
                        title_wl = f"λ={item['wavelength']:.2f}Å" if np.isfinite(item['wavelength']) else "λ=unknown"
                        ax0.set_title(f"Line {item['line_idx']} | {title_wl}")
                        ax0.grid(alpha=0.2)
                        ax0.legend(loc='best', fontsize=8)

                        # Right: residuals vs trace coord using the post-fit polynomial
                        det_w = float(item.get('detector_width', self.detector_dimensions.get('X', self.num_pixels)))
                        # Model at post-refine valid coords
                        model_at_valid = np.polyval(item['poly_coef'], item['y_valid'] / det_w) * det_w
                        resid = item['x_meas'] - model_at_valid
                        ax1.axhline(0.0, color='k', lw=0.8)
                        limit = float(self.good_arc_residual_limit)
                        ax1.axhline(+limit, color='orange', lw=0.8, ls='--')
                        ax1.axhline(-limit, color='orange', lw=0.8, ls='--', label=f"±limit ({limit:.2f}px)")
                        ax1.plot(item['y_valid'], resid, 'o', ms=3, alpha=0.9, label='post')
                        # Pre-refine residuals relative to same model
                        if 'y_valid_pre' in item and 'x_meas_pre' in item:
                            model_at_pre = np.polyval(item['poly_coef'], item['y_valid_pre'] / det_w) * det_w
                            resid_pre = item['x_meas_pre'] - model_at_pre
                            ax1.plot(item['y_valid_pre'], resid_pre, 'x', ms=4, alpha=0.8, color='tab:green', label='pre')
                        ax1.set_xlabel('Trace coord (Y)')
                        ax1.set_ylabel('Residual (px)')
                        ax1.set_title(f"MAD(post)={item['mad_resid']:.3f} px")
                        ax1.grid(alpha=0.2)
                        ax1.legend(loc='best', fontsize=8)

                        fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)
            except Exception as e:
                logger.error("Failed to write trace-fit diagnostics PDF: %s", e)

        return fitted_positions, residuals

    def _fit_wavelength_polynomials(self, fitted_trace_positions: Array, arc_pixels: Array) -> Array:
        # Determine good lines globally based on residuals between fitted positions and measurements
        res = np.abs(fitted_trace_positions - arc_pixels)
        line_res = mad_std(res, axis=0, ignore_nan=True)

        # Compute residual limit (using requested variable name good_arc_residual_limt) and derive mask
        good_arc_residual_limt = float(np.nanmedian(line_res) * 3.0)
        good_arc_residual_limit = good_arc_residual_limt  # keep both names for clarity/back-compat in this scope
        finite_mask = np.isfinite(line_res)
        good_lines = finite_mask & (line_res < good_arc_residual_limit)

        # Log summary statistics of line residuals and the chosen limit
        if np.any(finite_mask):
            lr = line_res[finite_mask]
            p5, p50, p95 = np.nanpercentile(lr, [5, 50, 95])
            logger.info(
                "Arc-line residuals: n_finite=%d, n_good=%d/%d, limit=%.4f px | min=%.4f, p50=%.4f, p95=%.4f, max=%.4f",
                int(np.count_nonzero(finite_mask)), int(np.count_nonzero(good_lines)), int(line_res.size),
                good_arc_residual_limit, float(np.nanmin(lr)), float(p50), float(p95), float(np.nanmax(lr))
            )
        else:
            logger.error("Arc-line residuals: all NaN; cannot select good lines")

        # If too few good lines, log a compact table of problematic lines and (optionally) save CSV
        n_good = int(np.count_nonzero(good_lines))
        if n_good < 3:
            n_total = int(line_res.size)
            over_by = line_res - good_arc_residual_limit
            # Build table rows for logging (sort by residual descending, show up to 20)
            indices = np.arange(n_total)
            order = np.argsort(np.where(np.isfinite(line_res), line_res, -np.inf))[::-1]
            rows_to_show = []
            for idx in order:
                if len(rows_to_show) >= 20:
                    break
                wl_val = float(self.arc_wavelengths[idx]) if idx < self.arc_wavelengths.size and np.isfinite(self.arc_wavelengths[idx]) else np.nan
                res_val = float(line_res[idx]) if np.isfinite(line_res[idx]) else np.nan
                is_good = bool(good_lines[idx])
                rows_to_show.append((idx, wl_val, res_val, good_arc_residual_limit, res_val - good_arc_residual_limit if np.isfinite(res_val) else np.nan, is_good))
            logger.warning("Too few good arc lines for wavelength fit: %d good out of %d (limit=%.4f px). Showing up to 20 lines (idx, wl_A, res_px, limit_px, over_by_px, is_good):",
                           n_good, n_total, good_arc_residual_limit)
            for r in rows_to_show:
                logger.warning("%4d, %9.2f, %7.4f, %7.4f, %+8.4f, %s", r[0], r[1], r[2], r[3], r[4], str(r[5]))

            # Save full CSV if out_folder is available
            if self.out_folder is not None:
                try:
                    out_csv = Path(self.out_folder) / "arc_line_residuals.csv"
                    data = np.column_stack([
                        indices,
                        self.arc_wavelengths if self.arc_wavelengths.size == n_total else np.full(n_total, np.nan, dtype=float),
                        line_res,
                        np.full(n_total, good_arc_residual_limit, dtype=float),
                        over_by,
                        good_lines.astype(int),
                    ])
                    header = "index,wavelength_A,line_residual_px,limit_px,over_by_px,is_good(1/0)"
                    np.savetxt(out_csv, data, delimiter=",", header=header, comments="", fmt=["%d","%.6f","%.6f","%.6f","%.6f","%d"])
                except Exception as e:
                    logger.error("Failed to write arc_line_residuals.csv: %s", e)

        wavelength_solution = np.full((self.num_fibers, self.num_pixels), np.nan, dtype=float)
        x = np.arange(self.num_pixels, dtype=float)
        wl = self.arc_wavelengths[good_lines]
        norm = 1000.
        for f in range(self.num_fibers):
            y = fitted_trace_positions[f, good_lines]
            m = np.isfinite(y)
            if np.count_nonzero(m) < 3:
                continue
            # Degree selection: prefer 5th, lower if not enough points
            for deg in (5, 4, 3, 2, 1):
                if np.count_nonzero(m) >= deg + 1:
                    coeff = np.polyfit(y[m] / norm, wl[m] / norm, deg)
                    wavelength_solution[f] = np.polyval(coeff, x / norm) * norm
                    break
        return wavelength_solution

# ------------------------- Public wrapper -------------------------

def get_wavelength_registration(
    spectra,
    trace_positions,
    valid_fibers,
    arc_pixel_guesses,
    arc_wavelengths,
    peak_threshold=5,
    reference_fiber_index=130,
    good_arc_residual_limit=0.2,
    binned=False,
    out_folder=None,
    detector_dimensions=None,
    *,
    enable_plots=True,
):
    calibrator = WavelengthCalibrator(
        spectra=spectra,
        trace_positions=trace_positions,
        valid_fibers=valid_fibers,
        arc_pixel_guesses=arc_pixel_guesses,
        arc_wavelengths=arc_wavelengths,
        peak_threshold=peak_threshold,
        reference_fiber_index=reference_fiber_index,
        good_arc_residual_limit=good_arc_residual_limit,
        binned=binned,
        out_folder=out_folder,
        detector_dimensions=detector_dimensions,
        enable_plots=enable_plots,
    )
    return calibrator.run()

def get_rectified_wavelength(start_wavelength, end_wavelength, detector_dimensions):
    """Generate a uniformly spaced wavelength grid.

    Creates a linear wavelength grid between specified start and end wavelengths
    with the number of points equal to the X dimension of the detector.

    Args:
        start_wavelength (float): Starting wavelength value
        end_wavelength (float): Ending wavelength value
        detector_dimensions (dict): Dictionary containing detector dimensions with keys:
            - X (int): Number of points in the wavelength grid (detector X dimension)
            - Y (int): Y dimension of the detector (not used)

    Returns:
        def_wavelength (ndarray): 1D rectified wavelength array with uniform spacing
    """
    def_wavelength = np.linspace(start_wavelength, end_wavelength, detector_dimensions['X'])
    return def_wavelength