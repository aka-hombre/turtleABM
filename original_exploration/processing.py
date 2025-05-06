import geopandas as gpd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from typing import Union
from shapely.geometry import shape, mapping, Polygon, LineString, box
from shapely.ops import unary_union, split
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap




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
current_location = (-78.95833000903956,32.034504148477325)
g = gpd.read_file('geojsondata/gspoly.geojson')
line = gpd.GeoSeries(LineString([current_location,(-78.96191293028112, 32.03808706971887)] ))
nline = gpd.GeoSeries(LineString([current_location,(-78.95779329244226, 32.045614147619176)] ))
'''

eastof_gs = gpd.read_file('geojsondata/shore_to_gulfstrem.geojson').to_crs("EPSG:4326")
gulfstream = gpd.read_file('geojsondata/gspoly.geojson').to_crs("EPSG:4326")

# Create a plot
fig, ax = plt.subplots(figsize=(16.8, 10))
cmap = LinearSegmentedColormap.from_list("ocean_gradient", ["#ADD8E6", "#00008B"])

# Plot each layer with specific colors for East of Gulf Stream and Gulf Stream
eastof_gs_plot = eastof_gs.plot(ax=ax, color=cmap(0.3), edgecolor='darkblue', linewidth=0.5)  
gulfstream_plot = gulfstream.plot(ax=ax, color=cmap(0.8), edgecolor='navy', linewidth=0.5)

# Create custom legend handles
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='East of Gulf Stream to Shore', 
           markerfacecolor='lightblue', markersize=10, linewidth=0),
    Line2D([0], [0], marker='o', color='w', label='Polygon of Gulf Stream', 
           markerfacecolor='darkslategray', markersize=10, linewidth=0)
]

# Add the legend manually
ax.legend(handles=legend_handles, loc='upper right')

# Add labels and title
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_title("Model Area with CRS: EPSG:4326")

plt.savefig('model_area.png', format='png', transparent=True, dpi=600)
plt.show()


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
#s
plt.title('Sea and FL')
plt.axis('off')
plt.show()
cargasso.plot(ax=ax, color='salmon', edgecolor='magenta', linewidth=0.5)oords = [
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