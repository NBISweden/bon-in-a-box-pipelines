import sys, json;
import os;
import pandas as pd
from pyproj import CRS, Transformer
import numpy as np


#import inputs
input = biab_inputs()
coordinates = input['coordinates']
input_crs = input['crs_input']
output_crs = input['crs_output']



# Define a function to transform coordinates
def transform_crs(longitudes, latitudes, input_crs_str, output_crs_str):
    """
    Transforms coordinates from an input CRS to an output CRS.

    """
    transformer = Transformer.from_crs(CRS(input_crs_str), CRS(output_crs_str), always_xy=True)
    transformed_longitudes, transformed_latitudes = transformer.transform(longitudes, latitudes)
    return transformed_longitudes, transformed_latitudes


# Run the transformation
print("Loading coordinates and transforming to the output CRS...")
coords_df = pd.read_csv(coordinates)

# Check that the first three columns have the proper names
if list(coords_df.columns[:3]) != ['id', 'latitude', 'longitude']:
    biab_error_stop("Error: The first three columns of the table must be named 'id', 'latitude', 'longitude'.")



#2. TRANSFORM TO THE SPECIFIED OUTPUT CRS

# Transform the coordinates to the specified output_crs
transformed_longitudes, transformed_latitudes = transform_crs(
    coords_df['longitude'].to_numpy(),
    coords_df['latitude'].to_numpy(),
    input_crs,
    output_crs
)

# Add the transformed coordinates to the DataFrame with a specific column name
coords_df[f"latitude_{output_crs.replace(':', '_')}"] = transformed_latitudes
coords_df[f"longitude_{output_crs.replace(':', '_')}"] = transformed_longitudes

# Save the transformed table for GEE
output_table_path = os.path.join(output_folder, "transformed_coordinates.csv")
coords_df.to_csv(output_table_path, index=False)

# Export the transformed table for the next script in the pipeline
biab_output("transformed_coordinates", output_table_path)
biab_output("crs_output", output_crs) # Export the main output CRS for the pipeline

print("Successfully transformed coordinates and saved the new table.")
