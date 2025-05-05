from shapely.geometry import Point, LineString
import mesa
from mesa.visualization import SolaraViz
import mesa_geo as mg
from mesa_geo.visualization import make_geospace_component
from mesa_geo import GeoSpace
import random, math
import geopandas as gpd

#for accurately moving agents
from geopy.distance import distance
from geopy import Point as GeoPyPoint
from pyproj import Geod

class Turtle(mg.GeoAgent):
    """
    Agent to model turtle behavior
    """
    def __init__(self, model, geometry, crs):
        super().__init__( model, geometry, crs)
    def __repr__(self):
        return f"Person {self.unique_id}"
    
    

    def move(self):
        geod = Geod(ellps="WGS84")
        distance_km = random.uniform(1.091, 1.239)
        distance_m = distance_km * 1000  # pyproj uses meters

        print(f"\nCurrent location: {self.geometry.x}, {self.geometry.y}")
        lon, lat = self.geometry.x, self.geometry.y

        # Random direction in degrees
        angle = random.uniform(0, 360)

        # Compute destination using pyproj.Geod
        dest_lon, dest_lat, _ = geod.fwd(lon, lat, angle, distance_m)
        new_position = Point(dest_lon, dest_lat)

        # Check if within polygon
        if not self.model.boundary_polygon.boundary.intersects(new_position) and self.model.boundary_polygon.contains(new_position):
            self.geometry = new_position
        else:
            print("First attempt failed — trying to bounce outward.")
            boundary = self.model.boundary_polygon.boundary
            nearest_boundary_point = boundary.interpolate(boundary.project(self.geometry))
            print(f"Nearest boundary point: {nearest_boundary_point}")

            # Compute geodesic azimuth from boundary to current point
            b_lon, b_lat = nearest_boundary_point.x, nearest_boundary_point.y
            _, azimuth_to_agent, _ = geod.inv(b_lon, b_lat, lon, lat)

            # Bounce back: reverse the azimuth
            outward_azimuth = (azimuth_to_agent + 180) % 360

            # Compute new destination in outward direction
            bounce_lon, bounce_lat, _ = geod.fwd(lon, lat, outward_azimuth, distance_m)
            new_position = Point(bounce_lon, bounce_lat)
            print(f"Second attempt destination: {new_position}")

            if not self.model.boundary_polygon.boundary.intersects(new_position) and self.model.boundary_polygon.contains(new_position):
                self.geometry = new_position
                print("Second attempt succeeded.")
            else:
                print("Second attempt failed.")

    def step(self):
        self.move()



class MovingModel(mesa.Model):
    def __init__(self, num_agents=15):
        super().__init__()

        self.gdf = gpd.read_file("geojsondata/gspoly.geojson")
        self.boundary_polygon = self.gdf.union_all()

        self.space = GeoSpace(crs=self.gdf.crs, warn_crs_conversion=False)
        self.num_agents = num_agents
        self.start_lat = 26
        self.start_lon = -79
        self.start_noise = 1

        for _ in range(num_agents):
            while True:
                x = self.start_lon+random.uniform(-self.start_noise, self.start_noise)
                y = self.start_lat+random.uniform(-self.start_noise, self.start_noise)
                point = Point(x, y)
                if self.boundary_polygon.contains(point):
                    agent = Turtle(self, point, self.gdf.crs)
                    self.space.add_agents(agent)
                    break
    def step(self):
        self.agents.do("step")



model = MovingModel()

def Turtle_draw(agent):
    return {"color": "Green", "radius": 10} if isinstance(agent, Turtle) else {}

model_params = {
    "num_agents": {
        "type": "SliderInt",
        "value": 1,
        "label": "Number of Agents",
        "min": 1,
        "max": 100,
        "step": 1,
    }
}

page = SolaraViz(
    model,
    name="Turtle ABM Model",
    model_params=model_params,
    components=[
        make_geospace_component(
            Turtle_draw,
            zoom=4,
            scroll_wheel_zoom=False
        ),
    ],
)
