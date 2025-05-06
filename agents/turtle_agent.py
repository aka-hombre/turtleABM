import mesa
import mesa_geo as mg

from shapely.geometry import Point
import random
from mesa.visualization import SolaraViz
from pyproj import Geod

class Turtle(mg.GeoAgent):
    """
    Agent to model turtle behavior. v1 incorperates rudimentry movement.
    """
    def __init__(self, model,geometry, deg_min, deg_max, crs):
        super().__init__( model, geometry, crs)
        self.deg_min = deg_min 
        self.deg_max = deg_max 
    def __repr__(self):
        return f"Person {self.unique_id}"
    
    

    def move(self):
        """
        This funciton moves a turtle over one hour of time. 
        Assumptions: 
            There is no prefrence for angle or direction that a turtle will move in. 
        """
        geod = Geod(ellps="WGS84")
        distance_km = random.uniform(1.091, 1.239) #This is picking uniformly at random from the data in turtle_data.ipynb
        distance_m = distance_km * 1000  #pyproj uses meters

        lon, lat = self.geometry.x, self.geometry.y

        # Random direction in degrees
        angle = random.uniform(self.deg_min, self.deg_max)

        # Compute destination using pyproj.Geod
        dest_lon, dest_lat, _ = geod.fwd(lon, lat, angle, distance_m)
        new_position = Point(dest_lon, dest_lat)

        # Check if within polygon, but not the boundary
        if not self.model.boundary_polygon.boundary.intersects(new_position) and self.model.boundary_polygon.contains(new_position): 
            self.geometry = new_position
        else:
            #print("First attempt failed — trying to bounce outward.")
            #At this point a turtle is close to a boundary. This will move the turtle away from the boundary
            boundary = self.model.boundary_polygon.boundary
            nearest_boundary_point = boundary.interpolate(boundary.project(self.geometry))
            #print(f"Nearest boundary point: {nearest_boundary_point}")

            #This ensures turtle will move in reverse on the earth's surface
            # Compute geodesic azimuth from boundary to current point
            b_lon, b_lat = nearest_boundary_point.x, nearest_boundary_point.y
            _, azimuth_to_agent, _ = geod.inv(b_lon, b_lat, lon, lat)

            # Bounce back: reverse the azimuth
            outward_azimuth = (azimuth_to_agent + 180) % 360

            # Compute new destination in outward direction
            bounce_lon, bounce_lat, _ = geod.fwd(lon, lat, outward_azimuth, distance_m)
            new_position = Point(bounce_lon, bounce_lat)
            #print(f"Second attempt destination: {new_position}")

            if not self.model.boundary_polygon.boundary.intersects(new_position) and self.model.boundary_polygon.contains(new_position):
                self.geometry = new_position
                #print("Second attempt succeeded.")
            else:
                print("Failed movement.")

    def step(self):
        self.move()