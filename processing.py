import geopandas as gpd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from typing import Union
from shapely.geometry import shape, mapping, Polygon, LineString, box
from shapely.ops import unary_union, split



def cutoutleft(left, right):
    """
    Takes a left region (polygon), a right region (polygon), cuts them and returs the new polygon with the

    """
    left_geo = left.geometry.union_all()
    right_geo = right.geometry.union_all()
    overlap_left = right_geo.intersection(left_geo)
    newshape = right_geo.difference(overlap_left)
    return gpd.GeoDataFrame(geometry=[newshape], crs=right.crs)

def cutoutright(left, right):
    """
    Takes a left region (polygon), a right region (polygon or multilinestring),
    cuts them and returns the new polygon with the right parts removed.
    
    Args:
        left: GeoDataFrame with polygon geometry
        right: GeoDataFrame with polygon or MultiLineString geometry
    
    Returns:
        GeoDataFrame with the modified left geometry
    """
    left_geo = left.geometry.union_all()
    
    # Handle both Polygon and MultiLineString cases
    right_geo = right.geometry.union_all()
    
    # Initialize with the original left geometry
    result = left_geo
    
    # If right is a MultiLineString, we need to handle each line separately
    if right_geo.geom_type == 'MultiLineString':
        for line in right_geo.geoms:
            line_as_poly = line.buffer(1e-8)
            result = result.difference(line_as_poly)
    else:
        # Original polygon difference operation
        overlap = left_geo.intersection(right_geo)
        result = result.difference(overlap)
    
    return gpd.GeoDataFrame(geometry=[result], crs=left.crs)


def get_bounds(filename):
    area = gpd.read_file(filename)
    area = area.union_all()
    return area.bounds

print('Bounds:'+str(get_bounds('geojsondata/model_area.geojson')))


'''
#right = gpd.read_file('geojsondata/gspoly.geojson')

#newmiddle = cutoutleft(left, middle)
#newmiddle = cutoutright(newmiddle, right)
#newmiddle.to_file("shore_to_gulfstrem.geojson", driver="GeoJSON")

fig, ax = plt.subplots(figsize=(10, 8))
#left.plot(ax=ax, color='lightblue', edgecolor='darkblue', linewidth=0.5)

newshore.plot(ax=ax, color='salmon', edgecolor='magenta', linewidth=0.5)
#lgs.plot(ax=ax, color='lightblue', edgecolor='darkblue', linewidth=0.5)
#right.plot(ax=ax, color='orange', edgecolor='red', linewidth=0.5)
#sargasso.plot(ax=ax, color='salmon', edgecolor='magenta', linewidth=0.5)
plt.title('Sea and FL')
plt.axis('off')
plt.show()
coords = [
    (-76.57702860978428, 38.514803910605934),
    (-67.4271326757245, 48.01014307612323)
]

slice_line = LineString(coords)
slice_gdf = gpd.GeoDataFrame(geometry=[slice_line], crs="EPSG:4326")

gs = gpd.read_file('geojsondata/gulfstream.geojson').to_crs('EPSG:32617')

shore = gpd.read_file('geojsondata/shore_to_gulfstrem.geojson')

newshore = horizontal_slice(shore, gs)
#lgs = gpd.read_file('geojsondata/leftofgs.geojson').to_crs('EPSG:32617')
#sargasso = gpd.read_file('geojsondata/sargasso.geojson').to_crs('EPSG:32617')

#left = gpd.read_file('geojsondata/helper/custom.geojson')
#middle= gpd.read_file('ge
'''