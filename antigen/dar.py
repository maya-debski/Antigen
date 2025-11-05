import numpy as np
import logging
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroplan import Observer, FixedTarget

logger = logging.getLogger('antigen.dar')

class DARModel:
    """Differential Atmospheric Refraction model.

    This class provides utilities to compute DAR vectors under site/atmospheric
    conditions and to represent a fitted polynomial model of dx, dy versus
    wavelength that can be applied to fiber coordinates.

    The class mirrors the older module-level API. Helper physics functions and
    site defaults are encapsulated here, while legacy module-level names remain
    for backward compatibility.
    """
    # ---- Site: Harlan J. Smith Telescope (McDonald Observatory, Mt. Locke) ----
    DEFAULT_LAT = 30.671745  # deg
    DEFAULT_LON = -104.021973  # deg  (255.978027 E)
    DEFAULT_ELEV_M = 2007.7  # meters

    def __init__(
        self,
        coeffs_x=None,
        coeffs_y=None,
        wave=None,
        ref_wave=5000.0,
        header=None,
        ra_deg=None,
        dec_deg=None,
        obstime=None,
        pressure_hPa=600.0,
        temperature_C=7.0,
        rel_humidity=0.20,
        altitude_deg=None,
        airmass=None,
        location=None,
        poly_order=2,
    ):
        """Initialize DARModel.

        You can initialize with precomputed polynomial coefficients (legacy
        behavior), or omit them and provide sufficient information to compute
        them from observing conditions. FITS headers are supported.

        Args:
            coeffs_x (array-like, optional): Polynomial coefficients for dx
                (np.polyval format). If provided along with ``coeffs_y``, these
                are used directly.
            coeffs_y (array-like, optional): Polynomial coefficients for dy
                (np.polyval format).
            wave (array-like, optional): Wavelength grid in Angstroms. Required
                if coefficients must be fitted and cannot be inferred from
                ``header``.
            ref_wave (float, optional): Reference wavelength in Angstroms for
                DAR computation. Defaults to 5000.0.
            header (Mapping, optional): FITS-like header providing keywords for
                wavelength solution and observing conditions. If ``wave`` is
                not given, it will be derived from standard WCS keywords when
                available.
            ra_deg, dec_deg, obstime, pressure_hPa, temperature_C, rel_humidity,
            altitude_deg, airmass, location (optional): Parameters forwarded to
                ``compute_dar_ifu`` when computing coefficients.
            poly_order (int, optional): Polynomial order for fitting dx, dy vs
                wavelength. Defaults to 2.

        Raises:
            ValueError: If neither coefficients are supplied nor sufficient
                information is available to compute them.
        """
        # If coefficients are provided, use them (backward compatible path)
        if coeffs_x is not None and coeffs_y is not None:
            self.coeffs_x = np.asarray(coeffs_x)
            self.coeffs_y = np.asarray(coeffs_y)
            return

        # Otherwise, we need to compute coeffs from provided inputs/header
        # 1) Derive wavelength grid if not supplied
        if wave is None and header is not None:
            wave = self._wave_from_header(header)
        if wave is None:
            raise ValueError("DARModel requires coeffs_x/coeffs_y or enough information (wave/header + conditions) to compute them.")

        # 2) Pull observing conditions either from arguments or header
        ra_val, dec_val, time_val = ra_deg, dec_deg, obstime
        P_val, T_val, RH_val = pressure_hPa, temperature_C, rel_humidity
        alt_val, am_val = altitude_deg, airmass
        loc_val = location if location is not None else (self.DEFAULT_LAT, self.DEFAULT_LON, self.DEFAULT_ELEV_M)

        if header is not None:
            # Only fill missing items from header
            if ra_val is None or dec_val is None:
                ra_hdr, dec_hdr = self._radec_from_header(header)
                ra_val = ra_val if ra_val is not None else ra_hdr
                dec_val = dec_val if dec_val is not None else dec_hdr
            if time_val is None:
                time_val = self._time_from_header(header)
            if am_val is None and alt_val is None:
                alt_hdr, am_hdr = self._altair_from_header(header)
                alt_val = alt_val if alt_val is not None else alt_hdr
                am_val = am_val if am_val is not None else am_hdr
            if P_val is None:
                P_val = self._first_float(header, ["PRESSURE", "PRES", "BAROM", "PRESS", "PRESSURE_HPA"]) or pressure_hPa
            if T_val is None:
                T_val = self._first_float(header, ["TEMP", "TEMPERAT", "AMBTEMP", "TAMBIENT", "OUTTEMP"]) or temperature_C
            if RH_val is None:
                RH_val = self._first_float(header, ["HUMIDITY", "RH", "RELHUM", "REHUM"])
                if RH_val is not None and RH_val > 1.0:
                    RH_val = RH_val / 100.0
                if RH_val is None:
                    RH_val = rel_humidity

        # 3) Compute dx, dy and fit polynomials
        dx_dy_qZ = self.compute_dar_ifu(
            wave,
            ref_wave=ref_wave,
            ra_deg=ra_val,
            dec_deg=dec_val,
            obstime=time_val,
            pressure_hPa=P_val,
            temperature_C=T_val,
            rel_humidity=RH_val,
            altitude_deg=alt_val,
            airmass=am_val,
            location=loc_val,
            return_parallactic=True,
        )
        dx, dy, (q_rad, Z_rad) = dx_dy_qZ

        # Log diagnostics about header-derived inputs and computed angles
        try:
            obstime_str = (
                time_val.isot if isinstance(time_val, Time) else str(time_val)
            ) if (time_val is not None) else None
        except Exception:
            obstime_str = str(time_val) if time_val is not None else None

        logger.info(
            "Building DARModel from header/inputs: RA(deg)=%s, DEC(deg)=%s, obstime=%s, "
            "airmass=%s, altitude(deg)=%s, P(hPa)=%s, T(C)=%s, RH=%s, q(deg)=%.3f, Z(deg)=%.3f",
            f"{ra_val:.6f}" if ra_val is not None else None,
            f"{dec_val:.6f}" if dec_val is not None else None,
            obstime_str,
            f"{am_val:.3f}" if am_val is not None else None,
            f"{alt_val:.3f}" if alt_val is not None else None,
            f"{P_val:.2f}" if P_val is not None else None,
            f"{T_val:.2f}" if T_val is not None else None,
            f"{RH_val:.3f}" if RH_val is not None else None,
            np.deg2rad(0.0) if False else np.rad2deg(q_rad),
            np.rad2deg(Z_rad),
        )

        if (ra_val is None) or (dec_val is None) or (time_val is None):
            logger.warning(
                "Parallactic angle q defaulted to 0 because RA/DEC/obstime was incomplete."
            )

        self.coeffs_x, self.coeffs_y = self.fit_polynomials(wave, dx, dy, order=poly_order)

    def __call__(self, wavelength, fiber_x, fiber_y):
        """Apply DAR correction to fiber positions.

        Args:
            wavelength (array-like): Wavelengths in Angstroms.
            fiber_x (array-like): Input fiber x positions (arcsec).
            fiber_y (array-like): Input fiber y positions (arcsec).

        Returns:
            tuple:
                corrected_x (ndarray): fiber_x - dx(wavelength)
                corrected_y (ndarray): fiber_y + dy(wavelength)
        Note:
            Sign conventions:
            - Physics model computes (dx, dy) in East and North (sky) directions.
            - Instrument axes use +X = West, +Y = North. To move fibers opposite
              to the DAR in X we subtract dx (East), which shifts toward East in
              instrument coordinates. In Y, North is positive in both frames, so
              we add dy.
        """
        dx = np.polyval(self.coeffs_x, wavelength)
        dy = np.polyval(self.coeffs_y, wavelength)
        return fiber_x - dx, fiber_y + dy

    # --------------------------
    # Helper physics (encapsulated)
    @staticmethod
    def _sat_vapor_pressure_hPa(T_C):
        """Saturation vapor pressure over water (Tetens) in hPa.

        Args:
            T_C (float): Temperature in Celsius.

        Returns:
            float: Saturation vapor pressure in hPa.
        """
        return 6.1078 * 10.0 ** (7.5 * T_C / (T_C + 237.3))

    @classmethod
    def _n_air_filippenko(cls, lambda_um, P_hPa, T_C, RH):
        """Refractive index of moist air at local conditions.

        Implements Filippenko (1982) dry-air refractivity scaled to local
        pressure/temperature with a humidity reduction term.

        Args:
            lambda_um (array-like): Wavelengths in microns.
            P_hPa (float): Pressure in hPa.
            T_C (float): Temperature in Celsius.
            RH (float): Relative humidity [0,1].

        Returns:
            ndarray: Refractive index n(λ).
        """
        lam = np.array(lambda_um, dtype=float)
        sigma2 = (1.0 / lam) ** 2
        R_dry_SL = 64.328 + 29498.1 / (146.0 - sigma2) + 255.4 / (41.0 - sigma2)

        P_mmHg = P_hPa * 0.75006156
        T = T_C
        scale_PT = (P_mmHg / 720.883) * (288.15 / (273.15 + T))

        e_s_hPa = cls._sat_vapor_pressure_hPa(T)
        f_mmHg = (RH * e_s_hPa) * 0.75006156
        humid_term = (0.0624 - 0.000680 * sigma2) * f_mmHg / (1.0 + 0.003661 * T)

        R_local_ppm = R_dry_SL * scale_PT - humid_term
        n = 1.0 + R_local_ppm * 1e-6
        return n

    # --------------------------
    # Public utilities
    @classmethod
    def compute_dar_ifu(
        cls,
        wave,
        *,
        ref_wave=5000.0,
        ra_deg=None,
        dec_deg=None,
        obstime=None,
        pressure_hPa=600.0,
        temperature_C=7.0,
        rel_humidity=0.20,
        altitude_deg=None,
        airmass=None,
        location=(DEFAULT_LAT, DEFAULT_LON, DEFAULT_ELEV_M),
        return_parallactic=False,
    ):
        """Compute DAR vectors (dx, dy) under given conditions.

        Args:
            wave (array-like): Wavelengths in Angstroms.
            ref_wave (float, optional): Reference wavelength in Angstroms.
            ra_deg (float, optional): RA in deg (ICRS).
            dec_deg (float, optional): Dec in deg (ICRS).
            obstime (str | astropy.time.Time, optional): Observation time.
            pressure_hPa (float, optional): Pressure in hPa.
            temperature_C (float, optional): Temperature in Celsius.
            rel_humidity (float, optional): Relative humidity [0,1].
            altitude_deg (float, optional): Telescope altitude in degrees.
            airmass (float, optional): Airmass (sec z). Overrides altitude.
            location (tuple, optional): (lat_deg, lon_deg, elev_m).
            return_parallactic (bool, optional): Include (q, Z) if True.

        Returns:
            tuple: (dx, dy) or (dx, dy, (q, Z)).
        """
        wave = np.atleast_1d(wave).astype(float)
        ref_wave = float(ref_wave)

        lat, lon, elev_m = location
        observer = Observer(longitude=lon * u.deg, latitude=lat * u.deg, elevation=elev_m * u.m)

        # Zenith distance
        if airmass is not None:
            Z_rad = np.arccos(1.0 / np.asarray(airmass, dtype=float))
        elif altitude_deg is not None:
            Z_rad = np.deg2rad(90.0 - float(altitude_deg))
        else:
            if obstime is None or ra_deg is None or dec_deg is None:
                raise ValueError("Need (airmass) or (altitude_deg) or (ra_deg, dec_deg, obstime) to determine Z.")
            t = Time(obstime) if not isinstance(obstime, Time) else obstime
            target = FixedTarget(coord=SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs"))
            altaz = observer.altaz(t, target.coord)
            Z_rad = (90.0 * u.deg - altaz.alt).to_value(u.rad)

        # Parallactic angle (East of North)
        if (ra_deg is not None) and (dec_deg is not None) and (obstime is not None):
            t = Time(obstime) if not isinstance(obstime, Time) else obstime
            target = FixedTarget(coord=SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs"))
            q = observer.parallactic_angle(t, target.coord).to_value(u.rad)
        else:
            q = 0.0

        lam_um = wave * 1e-4
        lam0_um = ref_wave * 1e-4
        n_lam = cls._n_air_filippenko(lam_um, pressure_hPa, temperature_C, rel_humidity)
        n_lam0 = cls._n_air_filippenko(lam0_um, pressure_hPa, temperature_C, rel_humidity)

        dR_arcsec = 206265.0 * (n_lam - n_lam0) * np.tan(Z_rad)

        dx = dR_arcsec * np.sin(q)
        dy = dR_arcsec * np.cos(q)

        return (dx, dy, (q, Z_rad)) if return_parallactic else (dx, dy)

    # --------------------------
    # Header parsing helpers
    @staticmethod
    def _first_float(header, keys):
        for k in keys:
            if k in header:
                try:
                    return float(header[k])
                except Exception:
                    continue
        return None

    @classmethod
    def _wave_from_header(cls, header):
        """Construct wavelength array from FITS-like WCS keywords.

        Supported keys: CRVAL1, CDELT1 or CD1_1, CRPIX1, NAXIS1.
        Assumes Angstroms if units are unspecified.
        """
        if header is None:
            return None
        crval = header.get("CRVAL1")
        cdelt = header.get("CDELT1", header.get("CD1_1"))
        crpix = header.get("CRPIX1", 1.0)
        npix = header.get("NAXIS1")
        if crval is None or cdelt is None or npix is None:
            return None
        try:
            crval = float(crval)
            cdelt = float(cdelt)
            crpix = float(crpix)
            npix = int(npix)
        except Exception:
            return None
        indices = np.arange(1, npix + 1, dtype=float)
        wave = crval + (indices - crpix) * cdelt
        return wave

    @staticmethod
    def _radec_from_header(header):
        """Extract RA/DEC from header in degrees.

        Supports:
        - Numeric degrees (float or string convertible)
        - Sexagesimal strings for RA (hms) and DEC (dms), e.g., '12:34:56.7', '+12:34:56'
        """
        if header is None:
            return None, None
        ra = header.get("RA")
        dec = header.get("DEC")
        if ra is None or dec is None:
            return None, None
        # First try simple float conversion
        try:
            ra_val = float(ra)
            dec_val = float(dec)
            return ra_val, dec_val
        except Exception:
            pass
        # Try sexagesimal parsing
        try:
            # Common case: RA in hours, DEC in degrees
            c = SkyCoord(str(ra), str(dec), unit=(u.hourangle, u.deg), frame="icrs")
            return c.ra.deg, c.dec.deg
        except Exception:
            try:
                # Fallback: both in degrees-strings
                c = SkyCoord(str(ra), str(dec), unit=(u.deg, u.deg), frame="icrs")
                return c.ra.deg, c.dec.deg
            except Exception:
                logger.warning(f"Failed to parse RA/DEC from header values RA={ra!r}, DEC={dec!r}")
                return None, None

    @staticmethod
    def _time_from_header(header):
        if header is None:
            return None
        date_obs = header.get("DATE-OBS") or header.get("DATEOBS")
        mjd_obs = header.get("MJD-OBS") or header.get("MJDOBS")
        if mjd_obs is not None:
            try:
                return Time(float(mjd_obs), format="mjd")
            except Exception:
                pass
        if date_obs is not None:
            try:
                return Time(str(date_obs))
            except Exception:
                return None
        return None

    @staticmethod
    def _altair_from_header(header):
        if header is None:
            return None, None
        alt = header.get("ALTITUDE")
        if alt is None:
            alt = header.get("ALT")
        am = header.get("AIRMASS")
        if am is None:
            am = header.get("SECZ")
        try:
            alt_val = float(alt) if alt is not None else None
        except Exception:
            alt_val = None
        try:
            am_val = float(am) if am is not None else None
        except Exception:
            am_val = None
        return alt_val, am_val

    @staticmethod
    def fit_polynomials(wave, dx, dy, order=2):
        """Fit polynomials to dx, dy as functions of wavelength.

        Args:
            wave (array-like): Wavelengths in Angstroms.
            dx (array-like): East component of DAR (arcsec).
            dy (array-like): North component of DAR (arcsec).
            order (int, optional): Polynomial order. Defaults to 2.

        Returns:
            tuple:
                coeffs_x (ndarray): Coefficients for dx(wave).
                coeffs_y (ndarray): Coefficients for dy(wave).
        """
        wave = np.asarray(wave, dtype=float)
        coeffs_x = np.polyfit(wave, np.asarray(dx, dtype=float), order)
        coeffs_y = np.polyfit(wave, np.asarray(dy, dtype=float), order)
        return coeffs_x, coeffs_y