from shapely.geometry import Point
import mesa
from mesa.visualization import SolaraViz
import mesa_geo as mg
from mesa_geo.visualization import make_geospace_component
from mesa_geo import GeoSpace
import random
import geopandas as gpd

class Turtle(mg.GeoAgent):
    """
    Agent to model turtle behavior
    """
    def __init__(self, model, geometry, move, crs):
        super().__init__( model, geometry, crs)
        self.move = move
    def __repr__(self):
        return f"Person {self.unique_id}"
    
    def step(self):
        dx = random.uniform(-self.move, self.move)
        dy = random.uniform(-self.move, self.move)
        new_position = Point(self.geometry.x + dx, self.geometry.y + dy)
        if self.model.boundary_polygon.contains(new_position):
            self.geometry = new_position



class MovingModel(mesa.Model):
    def __init__(self, num_agents=15, move=5):
        super().__init__()

        self.gdf = gpd.read_file("geojsondata/model_area.geojson")
        self.boundary_polygon = self.gdf.union_all()

        self.space = GeoSpace(crs=self.gdf.crs, warn_crs_conversion=False)
        self.num_agents = num_agents

        for _ in range(num_agents):
            while True:
                minx, miny, maxx, maxy = self.boundary_polygon.bounds
                x = random.uniform(minx, maxx)
                y = random.uniform(miny, maxy)
                point = Point(x, y)
                if self.boundary_polygon.contains(point):
                    agent = Turtle(self, point, move, self.gdf.crs)
                    #agent.pos = point
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
        "value": 10,
        "label": "Number of Agents",
        "min": 1,
        "max": 100,
        "step": 1,
    },
    "move":{
        "type": "SliderInt",
        "value": 11,
        "label": "move",
        "min": 1,
        "max": 100000,
        "step": 100,
    }
}

page = SolaraViz(
    model,
    name="GeoSim",
    model_params=model_params,
    components=[
        make_geospace_component(
            Turtle_draw,
            zoom=4,
            scroll_wheel_zoom=False
        ),
    ],
)
