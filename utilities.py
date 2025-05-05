import xarray as xr
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from typing import Tuple, List, Optional
import numpy as np
from pyproj import Geod
import json
from processing import get_bounds

def getdvx(fullfilename,  
                  Longitude=None,
                  Latitude=None,)-> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[Tuple[float, float, float, float]]]:
    """
    Extract velocity data and positional grids from a netCDF file. 
    Note: This function works specifically for /thredds/ncss/GLBu0.08/expt_91.2/uv3 datastructure. 
    Upon submitting the request for the data, the request for the rectagualar grid was specified with coordinates. The original
    slicing logic was using the Gulf database, and working with that specific .nc data structure.: 
    

    Args:
        fullfilename: Path to the netCDF file containing ocean current data.
        Longitude: Slice object defining the longitude range (min, max). If None, uses full range.
        Latitude: Slice object defining the latitude range (min, max). If None, uses full range.
    Returns:
        lon2d: 2D array of longitudes.
        lat2d: 2D array of latitudes.
        u_data: 2D array of u-component velocities.
        v_data: 2D array of v-component velocities.
        bounds: Optional tuple of (lon_min, lon_max, lat_min, lat_max).
        Returns `None` if slicing fails due to invalid coordinates.


    Raises:
        KeyError: If required variables ('u_barotropic_velocity', 'v_barotropic_velocity') are missing.
        ValueError: If Latitude/Longitude slices are out of bounds.

    Example:
        >>> lonlat, velocities, bounds = getvelocities("currents.nc", Latitude=slice(29.5, 30.5))
        >>> print(velocities[0].shape)  # u-component shape
    """
    ds = xr.open_dataset(fullfilename)

   #  slice(0, 10) as default
    if Latitude is not None:
        lat_slice = slice(*Latitude)
    else:
        lat_slice = slice(None) 
    if Longitude is not None:
        lon_slice = slice(*Longitude)
    else:
        lon_slice = slice(None)
    try:
        u = ds['water_u']
        v = ds['water_v']
    except KeyError as e:
        print(f"Error: Variable not found in dataset - {str(e)}")
    except ValueError as e:
        # Handle case where the coordinate values are out of bounds
        if "Latitude" in str(e) or "Latitude" in str(e):
            print("Error: Latitude out of bounds")
        elif "Longitude" in str(e):
            print("Error: Longitude out of bounds")
        else:
            print(f"Error: {str(e)}")

    lat = u['lat'].values
    lon = u['lon'].values

    lon2d, lat2d = np.meshgrid(lon, lat)

    # Kate mansfield paper: "Turtles largely remained ino ceanicwaters off the Coninental Shelf (greater than 200 mdepth)"
    u_data = u.isel(depth=200, time=0).values
    v_data = v.isel(depth=200, time=0).values

    bounds = [lon.min(), lon.max(), lat.min(), lat.max()] if Latitude and Longitude else None
    return lon2d, lat2d, u_data, v_data, bounds

def plotdvx(lon: np.ndarray, 
            lat: np.ndarray, 
            du: np.ndarray, 
            dv: np.ndarray, 
            bounds: Optional[List[float]] = None,
            scale: Optional[float] = 5,
            skip: Optional[int] = None) -> None:
    """
    Plots 2D velocity vectors on a geographic map using longitude and latitude grids.

    Args:
        lon: 2D array of longitudes (meshgrid format).
        lat: 2D array of latitudes (meshgrid format).
        du: 2D array of u-component velocities (eastward).
        dv: 2D array of v-component velocities (northward).
        bounds: Optional list/tuple of [lon_min, lon_max, lat_min, lat_max] for plot extent.
               If None, uses the min/max of input lon/lat arrays.
        scale: Optional float values of scale for quiver function.
        skip: Optional integer specifying stride for plotting vectors (plot every nth point).
              If None or 1, plots all vectors. Must be ≥ 1.

    Returns:
        None: Displays the plot directly using matplotlib.

    Raises:
        ValueError: If input arrays have incompatible shapes, bounds are invalid, or skip < 1.
        TypeError: If bounds is provided but not in the correct format.
        RuntimeError: If plotting fails due to Cartopy/matplotlib issues.

    Example:
        >>> # Plot all vectors
        >>> plotdvx(lon2d, lat2d, u, v, bounds)
        >>> # Plot every 5th vector
        >>> plotdvx(lon2d, lat2d, u, v, bounds, skip=5)
    """
    try:
        # Validate inputs
        if lon.shape != lat.shape or lon.shape != du.shape or lon.shape != dv.shape:
            raise ValueError("All input arrays must have the same shape")
            
        if bounds is not None:
            if len(bounds) != 4:
                raise ValueError("bounds must contain exactly 4 elements [lon_min, lon_max, lat_min, lat_max]")
            if bounds[0] >= bounds[1] or bounds[2] >= bounds[3]:
                raise ValueError("Invalid bounds: min values must be less than max values")

        if skip is not None:
            if not isinstance(skip, int) or skip < 1:
                raise ValueError("skip must be an integer ≥ 1")

        # Create figure
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        try:
            ax.coastlines()
            ax.set_extent(bounds if bounds is not None else 
                         [lon.min(), lon.max(), lat.min(), lat.max()], 
                         crs=ccrs.PlateCarree())

            # Apply skip if specified
            if skip and skip > 1:
                lon_sub = lon[::skip, ::skip]
                lat_sub = lat[::skip, ::skip]
                du_sub = du[::skip, ::skip]
                dv_sub = dv[::skip, ::skip]
            else:
                lon_sub, lat_sub, du_sub, dv_sub = lon, lat, du, dv

            # Plot vectors
            quiver = ax.quiver(lon_sub, lat_sub, du_sub, dv_sub, 
                             scale=scale, transform=ccrs.PlateCarree())

            title = "Barotropic Velocity Vectors"
            if bounds is not None:
                title += f"\nLon.=({bounds[0]:.2f}, {bounds[1]:.2f}) Lat.=({bounds[2]:.2f}, {bounds[3]:.2f})"
            if skip and skip > 1:
                title += f" (showing every {skip}th vector)"
            plt.title(title)
            
            plt.show()
            
        except Exception as e:
            plt.close(fig)
            raise RuntimeError(f"Plotting failed: {str(e)}") from e
            
    except ValueError as e:
        print(f"ValueError: {str(e)}")
        raise
    except TypeError as e:
        print(f"TypeError: {str(e)}")
        raise
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        raise

def geoarea(bounds:np.ndarray)-> float:
    """
    Takes the area of a ploygon on the earth (m^2). 
    Note: this is approximation.

    Args:
        bounds:  An array containg (lon_min, lon_max, lat_min, lat_max)
    Returns:
        area: float of the area of the polygon created by the bounds
    Example:
        >>> lonlat, velocities, bounds = getvelocities("currents.nc", Latitude=slice(29.5, 30.5))
        >>> print(geoarea(bounds))
    """
    lon_min, lon_max, lat_min, lat_max = bounds
    
    # Define polygon edges (closed loop)
    lons = [lon_min, lon_min, lon_max, lon_max, lon_min] 
    lats = [lat_min, lat_max, lat_max, lat_min, lat_min]
    
    # Compute area using cartopy's
    geod = Geod(ellps="WGS84")
    area, _ = geod.polygon_area_perimeter(lons, lats)
    
    return abs(area)



'''
# Example usage

with open("geojsondata/gspoly.geojson", "r") as f:
    geojson_data = json.load(f)


for feature in geojson_data["features"]:
    if feature.get("id") == 2:
        coords = feature["geometry"]["coordinates"]
        feature["geometry"]["coordinates"] = coords[::-1]  # Reverse the outer list
        break

# Save to file
with open("modified.geojson", "w") as f:
    json.dump(geojson_data, f, indent=2)


b = get_bounds('geojsondata/model_area.geojson')
lon, lat, du, dv, bounds = getdvx("hycom2016/expt_91_uv3z.nc", Longitude=(b[0], b[2]), Latitude=(b[1], b[3]))

plotdvx(lon, lat, du, dv, bounds, 10, 15)

import xarray as xr
from scipy.interpolate import RegularGridInterpolator
import numpy as np

class CurrentField:
    def __init__(self, nc_file, u_var="u", v_var="v", lat_var="latitude", lon_var="longitude", time_index=0):
        # Load data from .nc file
        self.ds = xr.open_dataset(nc_file)

        # Extract variables
        self.lats = self.ds[lat_var].values
        self.lons = self.ds[lon_var].values

        # Extract U and V components for one time step (if time dimension exists)
        u_data = self.ds[u_var]
        v_data = self.ds[v_var]

        if "time" in u_data.dims:
            u_data = u_data.isel(time=time_index)
            v_data = v_data.isel(time=time_index)

        # Handle shape: (lat, lon)
        self.interp_u = RegularGridInterpolator((self.lats, self.lons), u_data.values, bounds_error=False, fill_value=0.0)
        self.interp_v = RegularGridInterpolator((self.lats, self.lons), v_data.values, bounds_error=False, fill_value=0.0)

    def get_vector(self, lon, lat):
        """
        Sample the current vector (u, v) at a given lon, lat.
        """
        point = np.array([[lat, lon]])  # interpolators use (lat, lon) order
        u = self.interp_u(point)[0]
        v = self.interp_v(point)[0]
        return u, v
'''