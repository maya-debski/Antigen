import numpy as np

class DARModel:
    """Differential Atmospheric Refraction model."""

    def __init__(self, coeffs_x, coeffs_y):
        """
        Args:
            coeffs_x (array-like): Polynomial coefficients for x, in np.polyval format.
            coeffs_y (array-like): Polynomial coefficients for y, in np.polyval format.
        """
        self.coeffs_x = coeffs_x
        self.coeffs_y = coeffs_y

    def __call__(self, wavelength, fiber_x, fiber_y):
        """Apply DAR correction to fiber positions."""
        dx = np.polyval(self.coeffs_x, wavelength)
        dy = np.polyval(self.coeffs_y, wavelength)

        return fiber_x - dx, fiber_y - dy