import geopandas as gpd
import mesa

from shapely.geometry import Point
from mesa_geo import GeoSpace
from agents.turtle_agent import Turtle
from mesa.datacollection import DataCollector

class MovingModel(mesa.Model):
    """
    Model for turtle movement, bounded by the gulf stream. Agents are assumed to be dropped off at the same location (-79, 29).

    """
    def __init__(self, num_agents=1, deg_min=0, deg_max=360, seed=None):
        super().__init__(seed=seed)

        self.gdf = gpd.read_file("geojsondata/gspoly_w_startingpoly.geojson").to_crs("EPSG:4326")
        self.boundary_polygon = self.gdf.union_all()

        self.space = GeoSpace(crs=self.gdf.crs, warn_crs_conversion=False)
        self.num_agents = num_agents
        self.start_lat = 26
        self.start_lon = -79
        #self.start_noise = 0.001

        for _ in range(num_agents):
            while True:
                x = self.start_lon 
                y = self.start_lat
                point = Point(x, y)
                if self.boundary_polygon.contains(point):
                    agent = Turtle(self, point, deg_min, deg_max, self.gdf.crs)
                    self.space.add_agents(agent)
                    break
        
        self.datacollector = DataCollector(
            agent_reporters={
                "lon": lambda a: a.geometry.x,
                "lat": lambda a: a.geometry.y
            }
        )

    def step(self):
        self.agents.do("step")
        self.datacollector.collect(self)

