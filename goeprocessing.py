import geopandas as gpd
import matplotlib.pyplot as plt
from typing import Union
from shapely.geometry import shape, mapping, Polygon, LineString, box
from shapely.ops import unary_union, split

def multiline_json_to_geojson(input_json_path: str, output_geojson_path: str) -> None:
    """
    Convert a JSON file with multiline strings to a Polygon GeoJSON file.
    
    Args:
        input_json_path: Path to the input JSON file containing multiline strings
        output_geojson_path: Path where the output GeoJSON will be saved
    """
    # Load the input JSON file
    with open(input_json_path, 'r') as f:
        data = json.load(f)
    
    # Check if geopandas is available
    try:
        import geopandas as gpd
        USE_GEOPANDAS = True
    except ImportError:
        USE_GEOPANDAS = False
    
    # Process the multiline strings into polygons
    features = []
    
    if isinstance(data, dict):
        data = [data]  # Convert single feature to list for uniform processing
    
    for feature in data:
        if 'geometry' not in feature or 'coordinates' not in feature['geometry']:
            continue
            
        coords = feature['geometry']['coordinates']
        
        # Handle both single polygon and multipolygon cases
        polygons = []
        if isinstance(coords[0][0][0], (list, tuple)):  # Multipolygon case
            for polygon_coords in coords:
                polygons.append(Polygon(polygon_coords[0], [hole[0] for hole in polygon_coords[1:]]))
        else:  # Single polygon case
            polygons.append(Polygon(coords[0], [hole[0] for hole in coords[1:]]))
        
        # Create geometry - MultiPolygon if multiple, Polygon if single
        geometry = MultiPolygon(polygons) if len(polygons) > 1 else polygons[0]
        
        # Create feature properties (excluding geometry fields)
        properties = {k: v for k, v in feature.items() if k not in ['geometry', 'coordinates']}
        
        features.append({
            'type': 'Feature',
            'geometry': geometry.__geo_interface__,
            'properties': properties
        })
    
    # Create GeoJSON FeatureCollection
    geojson = {
        'type': 'FeatureCollection',
        'features': features
    }
    
    # Save using geopandas if available (better formatting)
    if USE_GEOPANDAS:
        gdf = gpd.GeoDataFrame.from_features(geojson['features'])
        gdf.to_file(output_geojson_path, driver='GeoJSON')
    else:
        with open(output_geojson_path, 'w') as f:
            json.dump(geojson, f)


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

gs = gpd.read_file('geojsondata/gulfstream.geojson').to_crs('EPSG:32617')
#lgs = gpd.read_file('geojsondata/leftofgs.geojson').to_crs('EPSG:32617')
#sargasso = gpd.read_file('geojsondata/sargasso.geojson').to_crs('EPSG:32617')

left = gpd.read_file('geojsondata/helper/custom.geojson')
middle= gpd.read_file('geojsondata/helper/middle.geojson')
right = gpd.read_file('geojsondata/gspoly.geojson')

newmiddle = cutoutleft(left, middle)
#newmiddle = cutoutright(newmiddle, right)

fig, ax = plt.subplots(figsize=(10, 8))
#left.plot(ax=ax, color='lightblue', edgecolor='darkblue', linewidth=0.5)

newmiddle.plot(ax=ax, color='salmon', edgecolor='magenta', linewidth=0.5)
#lgs.plot(ax=ax, color='lightblue', edgecolor='darkblue', linewidth=0.5)
right.plot(ax=ax, color='orange', edgecolor='red', linewidth=0.5)
#sargasso.plot(ax=ax, color='salmon', edgecolor='magenta', linewidth=0.5)
plt.title('Sea and FL')
plt.axis('off')  
plt.show()