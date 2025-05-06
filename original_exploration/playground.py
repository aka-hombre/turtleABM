from scipy.interpolate import RegularGridInterpolator
from typing import Tuple, Optional
from utilities import getdvx
import numpy as np

class CurrentField:
    """
    Interpolates ocean current vectors (u, v) from a rectilinear NetCDF grid using scipy.
    
    This class enables spatially continuous sampling of current vectors using interpolation 
    over a preloaded velocity field, such as one extracted with `getdvx`.

    Attributes:
        u_interp: Interpolator for the u-component (east-west) velocity.
        v_interp: Interpolator for the v-component (north-south) velocity.
        bounds: Optional geographic bounding box of the field (lon_min, lon_max, lat_min, lat_max).
    
    Args:
        lon2d: 2D array of longitudes.
        lat2d: 2D array of latitudes.
        u_data: 2D array of u-component velocities.
        v_data: 2D array of v-component velocities.
        bounds: Optional tuple specifying the data bounds.

    Example:
        >>> lon2d, lat2d, u, v, bounds = getdvx("currents.nc")
        >>> field = CurrentField(lon2d, lat2d, u, v, bounds)
        >>> u_val, v_val = field.get_vector(-90.5, 29.8)
    """
    def __init__(
        self,
        lon2d: np.ndarray,
        lat2d: np.ndarray,
        u_data: np.ndarray,
        v_data: np.ndarray,
        bounds: Optional[Tuple[float, float, float, float]] = None
    ):
        self.bounds = bounds

        # Assume grid is rectilinear, so 1D extraction is valid
        self.lon = lon2d[0, :]
        self.lat = lat2d[:, 0]

        # Ensure shape consistency
        assert u_data.shape == (len(self.lat), len(self.lon)), "u_data shape mismatch"
        assert v_data.shape == (len(self.lat), len(self.lon)), "v_data shape mismatch"

        # Build interpolators (lat, lon) order
        self.u_interp = RegularGridInterpolator((self.lat, self.lon), u_data, bounds_error=False, fill_value=0.0)
        self.v_interp = RegularGridInterpolator((self.lat, self.lon), v_data, bounds_error=False, fill_value=0.0)

    def get_vector(self, x: float, y: float) -> Tuple[float, float]:
        """
        Sample the interpolated (u, v) current vector at a given geographic coordinate.
        
        Args:
            x: Longitude in decimal degrees.
            y: Latitude in decimal degrees.

        Returns:
            Tuple of (u, v) current components at the specified location. Returns (0.0, 0.0)
            if outside the field or if interpolation fails.
        """
        try:
            u_val = float(self.u_interp((y, x)))  # lat, lon
            v_val = float(self.v_interp((y, x)))
            return u_val, v_val
        except Exception as e:
            print(f"Interpolation failed at ({x}, {y}): {e}")
            return 0.0, 0.0

lon2d, lat2d, u_data, v_data, bounds = getdvx("hycom2016/expt_91_uv3z.nc")

# Create the CurrentField interpolator
current_field = CurrentField(lon2d, lat2d, u_data, v_data, bounds)

# Define a geographic point (longitude, latitude)
x,y = [-65,
          34]
# Get the vector
u, v = current_field.get_vector(x, y)

# Print the vector
print(f"Current vector at ({x}, {y}): u = {u:.4f}, v = {v:.4f}")