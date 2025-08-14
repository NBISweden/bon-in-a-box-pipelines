import pandas as pd
from geopy.geocoders import Nominatim
from shapely.geometry import Point, Polygon, MultiPolygon, shape, box
from pyproj import CRS, Transformer
import numpy as np
import requests
import json
import pycountry
import geopandas # Ensure geopandas is imported for spatial operations
import os


#Import inputs

input = biab_inputs()
coordinates = input['coordinates']
input_crs = input['crs_input']
output_crs = input['crs_output']
grid_space = input['grid_space'] # Spacing in km for grid points
generate_grid = input['generate_grid'] #Boolean to control grid generation
table_local = input['table_local'] #Boolean to generate a table for local geotiff cropping


# Functions

#A function to transform coordinates
def transform_crs(longitudes, latitudes, input_crs_str, output_crs_str):
    """
    Transforms coordinates from an input CRS to an output CRS.

    Args:
        longitudes (np.array or list): Array or list of longitudes.
        latitudes (np.array or list): Array or list of latitudes.
        input_crs_str (str): String representation of the input CRS (e.g., "EPSG:4326").
        output_crs_str (str): String representation of the output CRS (e.g., "EPSG:3857").

    Returns:
        tuple: (transformed_longitudes, transformed_latitudes)
    """
    try:
        transformer = Transformer.from_crs(CRS(input_crs_str), CRS(output_crs_str), always_xy=True)
        transformed_longitudes, transformed_latitudes = transformer.transform(longitudes, latitudes)
        return transformed_longitudes, transformed_latitudes
    except Exception as e:
        biab_error_stop(f"Error transforming CRS from {input_crs_str} to {output_crs_str}: {e}")
        return None, None

#A function to convert the country codes from alpha-2 to alpha-3 using pycountries. Alpha-3 are used in geoboundaries
def convert_alpha2_to_alpha3(alpha2_code):
    try:
        country = pycountry.countries.get(alpha_2=alpha2_code)
        if country:
            return country.alpha_3
        else:
            print(f"Could not find alpha-3 for alpha-2 code: {alpha2_code}")
            return None
    except KeyError:
        print(f"Invalid alpha-2 code provided: {alpha2_2_code}")
        return None

#A function to get the country 2 letters code based on coordinates
def get_country_info_from_coords(latitude, longitude):
    """
    Gets the country alpha-2 code and bounding box for a given latitude and longitude. First does reverse geocoding
    to get the coutry and then forward geocoding to get the bounding box. The grid of points will be later generated for
    the bounding box and then cropped by the country shape.
    Returns (country_alpha2_code, bounding_box) where bounding_box is [south, north, west, east].
    Input latitude and longitude for this function MUST be in WGS84 (EPSG:4326).
    """
    geolocator = Nominatim(user_agent="country_grid_generator_app", timeout=10) # Increased timeout
    
    country_code_alpha2 = None
    country_name = None
    bounding_box = None

    try:
        # Step 1: Reverse geocode to get country name and alpha-2 code
        location_reverse = geolocator.reverse((latitude, longitude), exactly_one=True, language='en')
        if location_reverse and location_reverse.raw and 'address' in location_reverse.raw:
            country_code_alpha2 = location_reverse.raw['address'].get('country_code', '').upper()
            country_name = location_reverse.raw['address'].get('country')
            # print(f"Identified country (from reverse geocode): {country_name} (alpha-2: {country_code_alpha2})")
        else:
            print(f"Could not find country name or alpha-2 code for {latitude}, {longitude}")
            return None, None

        if country_name:
            # Step 2: Forward geocode the country name to get the country's overall bounding box
            # Use a longer timeout for this request as it might be for a large country
            location_country = geolocator.geocode(country_name, exactly_one=True, language='en', timeout=30) # Increased timeout
            if location_country and location_country.raw and 'boundingbox' in location_country.raw:
                bounding_box = location_country.raw['boundingbox'] # [south, north, west, east]
                # print(f"Identified country bounding box (from forward geocode): {bounding_box}")
            else:
                print(f"Could not find bounding box for country '{country_name}' via forward geocoding.")
        else:
            print("No country name found to perform forward geocoding for bounding box.")

        return country_code_alpha2, bounding_box

    except Exception as e:
        print(f"Error during geocoding: {e}")
        return None, None

#Function to get the GeoJSON polygon for a country using its alpha-3 code from geoboundaries.org
def get_country_geojson_polygon(alpha3_code):
    """
    Fetches the GeoJSON boundary for a given ISO 3166-1 alpha-3 country code
    from geoboundaries.org and returns it as a Shapely Polygon or MultiPolygon.
    """
    try:
        url = f"https://www.geoboundaries.org/api/current/gbOpen/{alpha3_code}/ADM0"
        response = requests.get(url, timeout=30)
        response.raise_for_status() # Raise an HTTPError for bad responses (4xx or 5xx)
        data = response.json()

        if "gjDownloadURL" not in data:
            print(f"Error: No GeoJSON download URL found for {alpha3_code}.")
            return None

        geojson_url = data["gjDownloadURL"]
        geojson_response = requests.get(geojson_url, timeout=60)
        geojson_response.raise_for_status()
        geojson_data = geojson_response.json()

        # Convert GeoJSON to Shapely geometry
        country_shape = shape(geojson_data['features'][0]['geometry'])
        print(f"Successfully fetched GeoJSON for {alpha3_code}.")
        return country_shape
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GeoJSON for {alpha3_code} from GeoBoundaries: {e}")
        return None
    except (json.JSONDecodeError, KeyError, IndexError) as e:
        print(f"Error parsing GeoJSON data for {alpha3_code}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while getting GeoJSON for {alpha3_code}: {e}")
        return None

#A function to determine the UTM region of a point. We will use this to build the grid of points later.
def get_utm_crs_for_point(latitude, longitude):
    """
    Determines the appropriate UTM CRS for a given latitude and longitude.
    Note: Input latitude and longitude for this function MUST be in WGS84 (EPSG:4326).
    """
    # UTM zones are 6 degrees wide.
    utm_zone = int(np.floor((longitude + 180) / 6) + 1)
    # Northern hemisphere is positive, Southern is negative
    southern_hemisphere = latitude < 0
    
    # Create the CRS string for the UTM zone
    utm_crs = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    
    # Adjust for southern hemisphere if necessary
    if southern_hemisphere:
        utm_crs = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    
    print(f"Determined UTM CRS: {utm_crs.to_string()}")
    return utm_crs

#A function to generate a grid of points within a bounding box in WGS84. We adjust manually for latitude to get approximately the same space between points. 
def generate_grid_points_in_wgs84_bounding_box(bounding_box_coords, spacing_km=100):
    """
    Generates a grid of points within a WGS84 bounding box,
    with an approximate spacing in kilometers, adjusting for latitude.

    I decided to use WGS84 because is the most common and global CRS, and the one in which we get the bounding box per country.
    The problem is that the grid distances change with latitude, but we approximate to how many degrees are the grid spacing
    through latitudes. At the equator, 1 degree latitude is approximately 111 km, and it shrinks towards the poles.

    We start creating the points at the south west corner of the bounding box. We create until the east edge 
    for a specific latitude, and then we move to the next latitude until we reach the north edge of the bounding box.

    bounding_box_coords: [south, north, west, east] in degrees.
    """
    south, north, west, east = map(float, bounding_box_coords)

    grid_points = []

    # Approximate degrees per 100km at the equator
    # 1 degree latitude is approx 111 km
    lat_step_deg = spacing_km / 111.0 
    
    current_lat = south #Start at the south edge of the bounding box
    while current_lat <= north: #Create points until the north edge
        # Calculate longitude step at current latitude
        # 1 degree longitude is approx 111 km * cos(latitude)
        # Avoid division by zero or very small numbers near poles
        cos_lat = np.cos(np.radians(current_lat)) #Numpy expects radians
        if abs(cos_lat) < 1e-6: # If very close to pole, use a small fixed step or skip
            lon_step_deg = lat_step_deg # Fallback if cos_lat is too small
        else:
            lon_step_deg = spacing_km / (111.0 * cos_lat)

        current_lon = west
        while current_lon <= east:
            point = Point(current_lon, current_lat) # Shapely Point expects (x, y) = (lon, lat)
            grid_points.append(point)
            current_lon += lon_step_deg
        current_lat += lat_step_deg
    
    print(f"Generated {len(grid_points)} approximate grid points within the bounding box (WGS84).")
    return grid_points

# Create a main function putting everything together

def create_country_grid(user_coordinates_df, grid_spacing_km=100, input_crs="EPSG:4326", output_crs="EPSG:4326", generate_grid=True):
    """
    Main function to create grids of points for all unique countries identified in the table of coordinates.
    It expands the most recent 'start_date' and 'end_date' to the grid of points. 
    The grid points are cropped by the actual country shape fetched online.

    """
    # Check that the first three columns have the proper names
    if list(user_coordinates_df.columns[:3]) != ['id', 'latitude', 'longitude']:
        biab_error_stop("Error: The first three columns of the input DataFrame must be named 'id', 'latitude', 'longitude'.")

    # Make a copy to avoid modifying the original DataFrame directly
    processed_user_coords_df = user_coordinates_df.copy()

    # Check for 'start_date' and 'end_date' columns and find the most recent dates
    has_date_cols = False
    most_recent_start_date = None
    most_recent_end_date = None

    if 'start_date' in processed_user_coords_df.columns and 'end_date' in processed_user_coords_df.columns:
        has_date_cols = True
        try:
            # Convert to datetime objects for proper comparison
            processed_user_coords_df['start_date_dt'] = pd.to_datetime(processed_user_coords_df['start_date'])
            processed_user_coords_df['end_date_dt'] = pd.to_datetime(processed_user_coords_df['end_date'])
            
            most_recent_start_date = processed_user_coords_df['start_date_dt'].max().strftime('%Y-%m-%d')
            most_recent_end_date = processed_user_coords_df['end_date_dt'].max().strftime('%Y-%m-%d')
            print(f"Found 'start_date' and 'end_date' columns. Most recent start_date: {most_recent_start_date}, end_date: {most_recent_end_date}")
        except Exception as e:
            print(f"Warning: Could not parse 'start_date' or 'end_date' to datetime. Dates will not be propagated. Error: {e}")
            has_date_cols = False
        finally:
            # Delete temporary datetime columns
            if 'start_date_dt' in processed_user_coords_df.columns:
                processed_user_coords_df = processed_user_coords_df.drop(columns=['start_date_dt'])
            if 'end_date_dt' in processed_user_coords_df.columns:
                processed_user_coords_df = processed_user_coords_df.drop(columns=['end_date_dt'])


    # Transform the input coordinates to WGS84 for geocoding. No matter the input CRS, we need WGS84 for geocoding.
    longitudes_wgs84, latitudes_wgs84 = transform_crs(
        processed_user_coords_df['longitude'].to_numpy(),
        processed_user_coords_df['latitude'].to_numpy(),
        input_crs,
        "EPSG:4326" # Geopy (Nominatim) requires WGS84
    )
    if longitudes_wgs84 is None: # Check for transformation errors
        return None

    # Transform the input coordinates directly to the desired output_crs for new columns
    longitude_trans_output_crs, latitude_trans_output_crs = transform_crs(
        processed_user_coords_df['longitude'].to_numpy(),
        processed_user_coords_df['latitude'].to_numpy(),
        input_crs,
        output_crs
    )
    if longitude_trans_output_crs is None: # Check for transformation errors
        return None

    # Add the transformed user coordinates to the initial dataframe
    processed_user_coords_df[f"latitude_{output_crs.replace(':', '_')}"] = latitude_trans_output_crs
    processed_user_coords_df[f"longitude_{output_crs.replace(':', '_')}"] = longitude_trans_output_crs

    #If we do not need to generate the greed, we can return the coordinates and finish here
    if not generate_grid:
        print("Grid generation skipped as requested. Returning only transformed user coordinates.")
        return processed_user_coords_df

    # Identify unique countries from all user coordinates
    #I will iterate through the coordinates table (after transformed to WGS84) and fill a dictionary with the country alpha-2 codes and their bounding box coordinates. 
    unique_countries = {} # Stores {alpha2_code: bounding_box_coords}
    for i in range(len(latitudes_wgs84)):
        lat_wgs84 = latitudes_wgs84[i]
        lon_wgs84 = longitudes_wgs84[i]
        
        country_alpha2_code, bounding_box_coords = get_country_info_from_coords(lat_wgs84, lon_wgs84)
        
        if country_alpha2_code and bounding_box_coords:
            if country_alpha2_code not in unique_countries:
                unique_countries[country_alpha2_code] = bounding_box_coords
                print(f"Identified unique country: {country_alpha2_code} with bounding box: {bounding_box_coords}")
        else:
            print(f"Could not identify country or bounding box for coordinate: {lat_wgs84}, {lon_wgs84}. Skipping.")

    if not unique_countries:
        print("No unique countries identified from the input coordinates. Cannot generate grids.")
        return processed_user_coords_df # Return user data even if no grids generated

    all_generated_grid_dfs = []

    # Generate grid for each unique country
    for country_alpha2_code, bounding_box_coords in unique_countries.items():
        print(f"\n--- Generating and cropping grid for country (alpha-2): {country_alpha2_code} ---")
        
        # Convert alpha-2 to alpha-3 for GeoBoundaries API
        country_alpha3_code = convert_alpha2_to_alpha3(country_alpha2_code)
        if not country_alpha3_code:
            print(f"Could not get alpha-3 code for {country_alpha2_code}. Skipping grid generation for this country.")
            continue

        # Fetch country GeoJSON
        country_shape_wgs84 = get_country_geojson_polygon(country_alpha3_code)
        if country_shape_wgs84 is None:
            print(f"Could not fetch GeoJSON for {country_alpha3_code}. Skipping grid generation for this country.")
            continue

        # Generate grid points within the WGS84 bounding box (initial dense grid)
        initial_grid_points_wgs84 = generate_grid_points_in_wgs84_bounding_box(bounding_box_coords, grid_spacing_km)
        
        if not initial_grid_points_wgs84:
            print(f"No initial grid points generated within bounding box for {country_alpha2_code}. Skipping this country.")
            continue

        # Convert initial grid points to a GeoDataFrame
        initial_grid_gdf = geopandas.GeoDataFrame(
            {'geometry': initial_grid_points_wgs84},
            crs="EPSG:4326"
        )
        
        # Create a GeoDataFrame for the country shape
        country_shape_gdf = geopandas.GeoDataFrame(
            {'geometry': [country_shape_wgs84]},
            crs="EPSG:4326"
        )

        # Perform spatial join to filter grid points that are within the country shape
        # Use 'within' predicate for point-in-polygon test
        cropped_grid_gdf = geopandas.sjoin(initial_grid_gdf, country_shape_gdf, predicate='within', how='inner')
        
        # Extract the filtered Shapely Points
        cropped_grid_points_wgs84 = cropped_grid_gdf['geometry'].tolist()

        print(f"Cropped grid to {len(cropped_grid_points_wgs84)} points within the actual country shape for {country_alpha2_code}.")

        if not cropped_grid_points_wgs84:
            print(f"No grid points remained after cropping for {country_alpha2_code}. Skipping this country.")
            continue

        # Convert list of cropped Shapely Points to DataFrame
        grid_df_wgs84_cropped = pd.DataFrame([{'latitude': p.y, 'longitude': p.x} for p in cropped_grid_points_wgs84])

        # Transform the final cropped grid points to the specified output_crs
        final_grid_longitudes, final_grid_latitudes = transform_crs(
            grid_df_wgs84_cropped['longitude'].to_numpy(),
            grid_df_wgs84_cropped['latitude'].to_numpy(),
            "EPSG:4326", # Cropped grid points are currently in WGS84
            output_crs
        )
        if final_grid_longitudes is None:
            continue # Skip this country if transformation fails

        country_grid_data = {
            'id': [f"grid_{country_alpha2_code}_{i}" for i in range(len(final_grid_latitudes))], # Unique IDs per country grid
            'latitude': np.nan, # Original latitude for grid points is not meaningful
            'longitude': np.nan, # Original longitude for grid points is not meaningful
            f"latitude_{output_crs.replace(':', '_')}": final_grid_latitudes,
            f"longitude_{output_crs.replace(':', '_')}": final_grid_longitudes,
        }
        
        # Add start_date and end_date to grid points if present in user data
        if has_date_cols:
            country_grid_data['start_date'] = most_recent_start_date
            country_grid_data['end_date'] = most_recent_end_date

        country_grid_df = pd.DataFrame(country_grid_data)
        
        print(f"Generated {len(country_grid_df)} cropped grid points for {country_alpha2_code} in {output_crs}.")
        all_generated_grid_dfs.append(country_grid_df)

    # Combine processed user coordinates with all generated grid points
    if all_generated_grid_dfs:
        # Concatenate all country-specific grids first
        combined_grids_df = pd.concat(all_generated_grid_dfs, ignore_index=True)

        # Ensure all columns from processed_user_coords_df are present in combined_grids_df
        # This handles propagation of any extra columns from user_coordinates_df
        missing_cols_in_grids = set(processed_user_coords_df.columns) - set(combined_grids_df.columns)
        for col in missing_cols_in_grids:
            combined_grids_df[col] = np.nan
        
        # Reorder columns to match processed_user_coords_df for clean concatenation
        combined_grids_df = combined_grids_df[processed_user_coords_df.columns]

        combined_df = pd.concat([processed_user_coords_df, combined_grids_df], ignore_index=True)
    else:
        combined_df = processed_user_coords_df # Only user data if no grids were generated

    print(f"\nTotal combined points in final DataFrame: {len(combined_df)}")
    return combined_df

# RUN FUNCTION

coords_df = pd.read_csv(coordinates)

#Table for google earth engine
print("transforming coordinates to output CRS")
final_coords_df = create_country_grid(coords_df, grid_spacing_km=grid_space,
                                input_crs=input_crs, output_crs=output_crs,
                                generate_grid=generate_grid) 


output_crs_table=os.path.join(output_folder, "transformed_coord_grid.csv")
final_coords_df.to_csv(output_crs_table, index=False)
biab_output("coordinates_trans_gee", output_crs_table)

#Table for local geotiff

if table_local:
    final_coords_df_local = create_country_grid(coords_df, grid_spacing_km=grid_space,
                                    input_crs=input_crs, output_crs='EPSG:4326',
                                    generate_grid=generate_grid)
    output_crs_table_local=os.path.join(output_folder, "transformed_coord_grid_local.csv")
    final_coords_df_local.to_csv(output_crs_table_local, index=False)
    biab_output("coordinates_trans_local", output_crs_table_local) 

#Export the output CRS for the pipeline
biab_output("crs_output", output_crs)

