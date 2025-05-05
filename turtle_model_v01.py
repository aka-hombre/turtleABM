from shapely.geometry import Point
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

class Turtle(mg.GeoAgent):
    """
    Agent to model turtle behavior
    """
    def __init__(self, model, geometry, crs):
        super().__init__( model, geometry, crs)
    def __repr__(self):
        return f"Person {self.unique_id}"
    
    def move(self):
        print(f"{self.geometry}")
        angle = random.uniform(0, 2 * math.pi)
        distance_km = random.uniform(1.091, 1.239) 

        origin = GeoPyPoint(self.geometry.y, self.geometry.x)  # lat, lon

        destination = distance(kilometers=distance_km).destination(origin, angle)

        new_position = Point(destination.longitude, destination.latitude)

        if self.model.boundary_polygon.contains(new_position):
            self.geometry = new_position

    def step(self):
        #print(f"{self.geometry.x}, {self.geometry.y}")
        self.move()
        print(f"{self.geometry.x}, {self.geometry.y}")



class MovingModel(mesa.Model):
    def __init__(self, num_agents=15):
        super().__init__()

        self.gdf = gpd.read_file("geojsondata/model_area.geojson")
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
