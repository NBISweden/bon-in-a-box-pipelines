import pandas as pd
from geopy.geocoders import Nominatim
from shapely.geometry import Point, Polygon, MultiPolygon, shape, box
from pyproj import CRS, Transformer
import numpy as np
import requests
import json
import pycountry
import geopandas
import os
import sys

# Get all inputs from biab
input = biab_inputs()
coordinates = input['coordinates']
input_crs = input['crs_input']
output_crs = input['crs_output']
grid_space = input['grid_space']
country_code = input['country_code']
start_date = input['start_date']
end_date = input['end_date']
create_grid = input['create_grid']

print(f"Coordinates: {coordinates}")
# Helper Functions

# A function to transform coordinates
def transform_crs(longitudes, latitudes, input_crs_str, output_crs_str):
    """
    Transforms coordinates from an input CRS to an output CRS. We need this because no matter the input 
    CRS, we need to transform to WGS84 for reverse geocoding

    """
    transformer = Transformer.from_crs(CRS(input_crs_str), CRS(output_crs_str), always_xy=True)
    transformed_longitudes, transformed_latitudes = transformer.transform(longitudes, latitudes)
    return transformed_longitudes, transformed_latitudes

# A function to convert the country codes from alpha-2 to alpha-3 using pycountries. Alpha-3 are used in geoboundaries
def convert_alpha2_to_alpha3(alpha2_code):
    """
    Converts an alpha-2 country code to an alpha-3 country code using the pycountry library.
    the alpha-2 is good to get the bounding box, but we need to cut the grid within the bounding box
    with the country shape. This country shape we get it using the alpha3 code and the function get_country_geojson_polygon
    """
    try:
        country = pycountry.countries.get(alpha_2=alpha2_code)
        if country:
            return country.alpha_3
        else:
            print(f"Could not find alpha-3 for alpha-2 code: {alpha2_code}")
            return None
    except KeyError:
        print(f"Invalid alpha-2 code provided: {alpha2_code}")
        return None

# A function to get the country 2 letters code based on coordinates
def get_country_info_from_coords(latitude, longitude):
    """
    A function for reverse geocoding
    Gets the country alpha-2 code and bounding box for a given latitude and longitude.
    Returns (country_alpha2_code, bounding_box) where bounding_box is [south, north, west, east].
    Input latitude and longitude for this function MUST be in WGS84 (EPSG:4326). 
    """
    geolocator = Nominatim(user_agent="country_grid_generator_app", timeout=60) 
    
    country_code_alpha2 = None
    country_name = None
    bounding_box = None

    try:
        # Step 1: Reverse geocode to get country name and alpha-2 code
        location_reverse = geolocator.reverse((latitude, longitude), exactly_one=True, language='en')
        if location_reverse and location_reverse.raw and 'address' in location_reverse.raw:
            country_code_alpha2 = location_reverse.raw['address'].get('country_code', '').upper()
            country_name = location_reverse.raw['address'].get('country')
        else:
            print(f"Could not find country name or alpha-2 code for {latitude}, {longitude}")
            return None, None

        if country_name:
            # Step 2: Get the country's overall bounding box
            location_country = geolocator.geocode(country_name, exactly_one=True, language='en', timeout=60)
            if location_country and location_country.raw and 'boundingbox' in location_country.raw:
                bounding_box = location_country.raw['boundingbox'] # [south, north, west, east]
            else:
                print(f"Could not find bounding box for country '{country_name}' via forward geocoding.")
        else:
            print("No country name found to perform forward geocoding for bounding box.")

        return country_code_alpha2, bounding_box

    except Exception as e:
        print(f"Error during geocoding: {e}")
        return None, None

# Function to get the GeoJSON polygon for a country using its alpha-3 code from geoboundaries.org
def get_country_geojson_polygon(alpha3_code):
    """
    Fetches the GeoJSON boundary for a given alpha-3 country code
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
    #Some error cases for this function
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GeoJSON for {alpha3_code} from GeoBoundaries: {e}")
        return None
    except (json.JSONDecodeError, KeyError, IndexError) as e:
        print(f"Error parsing GeoJSON data for {alpha3_code}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while getting GeoJSON for {alpha3_code}: {e}")
        return None

# A function to generate a grid of points within a bounding box in WGS84.
def generate_grid_points_in_wgs84_bounding_box(bounding_box_coords, spacing_km=100):
    """
    Generates a grid of points within a WGS84 bounding box.
    bounding_box_coords: [south, north, west, east] in degrees.
    Because the distances for a certain degree changes across latitudes, 
    we adjust for that manually 
    """
    south, north, west, east = map(float, bounding_box_coords)

    grid_points = []
    lat_step_deg = spacing_km / 111.0 #If we assume a distance of 100km between points ,we approximate how many degrees 100km in latitude will be
    
    current_lat = south #We start selecting the lowest latitude
    #Now we build two loops, the first one is to calculate how many degrees 100km in longitude are 
    #This is depending on the latitude where we are
    while current_lat <= north:
        cos_lat = np.cos(np.radians(current_lat))
        if abs(cos_lat) < 1e-6:
            lon_step_deg = lat_step_deg
        else:
            lon_step_deg = spacing_km / (111.0 * cos_lat)
    #Here we start creating the grid points at the specified longitude based on the latitude we are
        current_lon = west
        while current_lon <= east:#The loop will run until the current longitude exceeds the eastern boundary
            point = Point(current_lon, current_lat)
            grid_points.append(point)
            current_lon += lon_step_deg
        current_lat += lat_step_deg
    
    print(f"Generated {len(grid_points)} approximate grid points within the bounding box (WGS84).")
    return grid_points

# Main function to create the grid, now with updated logic to handle different inputs
def create_country_grid(user_coordinates_df, input_crs, output_crs, grid_spacing_km=100, user_country_code=None, user_start_date=None, user_end_date=None):
    """
    Main function to create grids of points. It can use a coordinate table and/or a country code.
    If using a table, it checks if the table is already transformed. It then expands the grid points.
    If the user does not want to create the grid, it returns the original table and the output crs
    """

    if not create_grid:
        print("`create_grid` is FALSE. No grid will be generated.")
        if user_coordinates_df is not None:
            print("Returning the user-provided table for next steps of the pipeline.")
            return user_coordinates_df
        else:
            # If no grid is created and no table is provided (very unlikely case), we trigger an error.
            biab_error_stop("Error: `create_grid` is FALSE and no coordinate table was provided.")
            return None

    print("`create_grid` is TRUE. Proceeding with grid generation process.")
    #We create empty dataframe and empty objects to later take the user provided dates
    processed_user_coords_df = pd.DataFrame()
    has_date_cols = False
    final_start_date = None
    final_end_date = None
    
    # A dictionary that will contain the unique countries to create a grid across
    unique_countries = {}
    
    # Take dates provided by user if any
    if user_start_date and user_end_date:
        final_start_date = user_start_date
        final_end_date = user_end_date
        has_date_cols = True
        print(f"Using user-provided dates: {final_start_date} to {final_end_date}.")

    # Scenario 1: User provides a coordinate table. We then copy it and process it
    if user_coordinates_df is not None:
        print("Using provided coordinate table...")
        processed_user_coords_df = user_coordinates_df.copy()

        # Check for transformed coordinates first.
        input_lat_col = f"latitude_{input_crs.replace(':', '_')}"
        input_lon_col = f"longitude_{input_crs.replace(':', '_')}"

        if input_lat_col in processed_user_coords_df.columns and input_lon_col in processed_user_coords_df.columns:
            print(f"Found transformed coordinates ({input_lat_col}, {input_lon_col}). Transforming them to EPSG:4326 for geocoding.")
            geocoding_lons_wgs84, geocoding_lats_wgs84 = transform_crs(
                processed_user_coords_df[input_lon_col].to_numpy(),
                processed_user_coords_df[input_lat_col].to_numpy(),
                input_crs,
                "EPSG:4326"
            )
        else:
            print("Transformed coordinates not found. Using original 'latitude' and 'longitude' columns for geocoding.")
            # Transform the original coordinates to WGS84 for geocoding
            geocoding_lons_wgs84, geocoding_lats_wgs84 = transform_crs(
                processed_user_coords_df['longitude'].to_numpy(),
                processed_user_coords_df['latitude'].to_numpy(),
                input_crs,
                "EPSG:4326"
            )
        
        # Now, add WGS84 columns to the user's dataframe for consistency, if they don't exist
        if 'latitude' not in processed_user_coords_df.columns:
            processed_user_coords_df['latitude'] = geocoding_lats_wgs84
        if 'longitude' not in processed_user_coords_df.columns:
            processed_user_coords_df['longitude'] = geocoding_lons_wgs84
        
        if geocoding_lats_wgs84 is None:
            return None
        
        # Identify unique countries from all user coordinates, by reverse geocoding
        print("Identifying unique countries from user coordinates...")
        for i in range(len(geocoding_lats_wgs84)):
            lat_wgs84 = geocoding_lats_wgs84[i]
            lon_wgs84 = geocoding_lons_wgs84[i]
            country_alpha2_code, bounding_box_coords = get_country_info_from_coords(lat_wgs84, lon_wgs84)
            if country_alpha2_code and bounding_box_coords: #We add it to the dictionary 
                unique_countries[country_alpha2_code] = bounding_box_coords
        
        # In the case that the user does not provide the dates, to avoid an error on the pipeline, we can just add the last dates found on the table
        if not has_date_cols and 'start_date' in processed_user_coords_df.columns and 'end_date' in processed_user_coords_df.columns:
            has_date_cols = True
            try:
                processed_user_coords_df['start_date_dt'] = pd.to_datetime(processed_user_coords_df['start_date'])
                processed_user_coords_df['end_date_dt'] = pd.to_datetime(processed_user_coords_df['end_date'])
                final_start_date = processed_user_coords_df['start_date_dt'].max().strftime('%Y-%m-%d')
                final_end_date = processed_user_coords_df['end_date_dt'].max().strftime('%Y-%m-%d')
                print(f"Found date columns in table. Using most recent dates: {final_start_date} to {final_end_date}.")
            except Exception as e:
                print(f"Warning: Could not parse dates. Error: {e}")
                has_date_cols = False
            finally:
                if 'start_date_dt' in processed_user_coords_df.columns: processed_user_coords_df = processed_user_coords_df.drop(columns=['start_date_dt'])
                if 'end_date_dt' in processed_user_coords_df.columns: processed_user_coords_df = processed_user_coords_df.drop(columns=['end_date_dt'])

        # Transform the WGS84 points to output CRS and add new columns
        print(f"Transforming user input points from WGS84 to {output_crs}...")
        final_lons, final_lats = transform_crs(
            geocoding_lons_wgs84,
            geocoding_lats_wgs84,
            "EPSG:4326",
            output_crs
        )
        if final_lons is None: return None
        
        processed_user_coords_df[f"latitude_{output_crs.replace(':', '_')}"] = final_lats#We return the column with the final lat and lon values at the desired crs
        processed_user_coords_df[f"longitude_{output_crs.replace(':', '_')}"] = final_lons

    # Scenario 2: User provides one or more country codes (can be in addition to a table or alone)
    if user_country_code:
        # Ensure country_code is a list
        if isinstance(user_country_code, str):
            user_country_code = [user_country_code]
        #Iterate through the provided country codes
        for code in user_country_code:
            print(f"Adding grid for provided country code '{code}'.")
            # Perform forward geocoding on the country code to get the bounding box
            geolocator = Nominatim(user_agent="country_grid_generator_app", timeout=60)
            location_country = geolocator.geocode(code, exactly_one=True, language='en')
            if location_country and location_country.raw and 'boundingbox' in location_country.raw:
                bounding_box_coords = location_country.raw['boundingbox']
                unique_countries[code.upper()] = bounding_box_coords#We add the bounding box to the country code element of the dictionary
                print(f"Found bounding box for {code}: {bounding_box_coords}")
            else:
                print(f"Warning: Could not find bounding box for country code '{code}'. Skipping.")
                continue

    # If no coordinates or country code are provided, stop the process
    if not unique_countries:
        biab_error_stop("Error: No coordinate table or country code provided. Cannot generate grid.")
        return None

    all_generated_grid_dfs = []

    # Generate grid for each unique country by iterating through the dictionary
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
        #This step is very important, because the grid sometimes includes points in the sea or other regions
        #We intersect the grid with the country shape
        cropped_grid_gdf = geopandas.sjoin(initial_grid_gdf, country_shape_gdf, predicate='within', how='inner')
        cropped_grid_points_wgs84 = cropped_grid_gdf['geometry'].tolist()

        print(f"Cropped grid to {len(cropped_grid_points_wgs84)} points within the actual country shape for {country_alpha2_code}.")

        if not cropped_grid_points_wgs84:
            print(f"No grid points remained after cropping for {country_alpha2_code}. Skipping this country.")
            continue

        grid_df_wgs84_cropped = pd.DataFrame([{'latitude': p.y, 'longitude': p.x} for p in cropped_grid_points_wgs84])

        # Transform the final cropped grid points to the specified output_crs
        final_grid_longitudes, final_grid_latitudes = transform_crs(
            grid_df_wgs84_cropped['longitude'].to_numpy(),
            grid_df_wgs84_cropped['latitude'].to_numpy(),
            "EPSG:4326", # Cropped grid points are currently in WGS84
            output_crs
        )
        if final_grid_longitudes is None:
            continue

        country_grid_data = {#We create the grid data that will be the dataframe later
            'id': [f"grid_{country_alpha2_code}_{i}" for i in range(len(final_grid_latitudes))],
            f"latitude_{output_crs.replace(':', '_')}": final_grid_latitudes,
            f"longitude_{output_crs.replace(':', '_')}": final_grid_longitudes,
        }
        
        # Add start_date and end_date to grid points if present
        if has_date_cols:
            country_grid_data['start_date'] = final_start_date
            country_grid_data['end_date'] = final_end_date

        country_grid_df = pd.DataFrame(country_grid_data)
        
        print(f"Generated {len(country_grid_df)} cropped grid points for {country_alpha2_code} in {output_crs}.")
        all_generated_grid_dfs.append(country_grid_df)

    # Here we append the grids to the dataframe of the user
    final_df = pd.DataFrame()
    if not processed_user_coords_df.empty:
        all_dataframes = [processed_user_coords_df] + all_generated_grid_dfs
        final_df = pd.concat(all_dataframes, ignore_index=True, sort=False)
    elif all_generated_grid_dfs:
        final_df = pd.concat(all_generated_grid_dfs, ignore_index=True, sort=False)
    
    print(f"\nTotal combined points in final DataFrame: {len(final_df)}")
    return final_df

# RUN FUNCTION

coords_df = None
if coordinates is not None:
    coords_df = pd.read_csv(coordinates)
    print("Coordinates table found. Running grid generation based on table...")
    #Check the first three columns
    if list(coords_df.columns[:3]) != ['id', 'latitude', 'longitude']:
        biab_error_stop("Error: The first three columns of the table must be named 'id', 'latitude', 'longitude'.")
    #Now I check if the input table has dates in the right format. If they dont, but have other colums with "dates" I trigger an error to reformat table
    if 'start_date' not in coords_df.columns or 'end_date' not in coords_df.columns:
        if any('date' in col.lower() for col in coords_df.columns):
            biab_error_stop("Error: Date columns must be named 'start_date' and 'end_date' for processing.")

#After checking if there an input table or not, and if it is in the right format, I run the main function

final_coords_df = create_country_grid(
    coords_df, 
    input_crs=input_crs, 
    output_crs=output_crs,
    grid_spacing_km=grid_space,
    user_country_code=country_code,
    user_start_date=start_date,
    user_end_date=end_date
)

#Export in bon in a box format
if final_coords_df is not None:
    # Save the final combined table
    output_table_path = os.path.join(output_folder, "transformed_coord_grid.csv")
    final_coords_df.to_csv(output_table_path, index=False)
    
    # Export the final table and the output CRS
    biab_output("transformed_grid", output_table_path)
    biab_output("crs_output", output_crs)
    print("Final grid and coordinates table successfully created.")
