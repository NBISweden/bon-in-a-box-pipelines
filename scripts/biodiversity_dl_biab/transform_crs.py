import sys, json;
import os;
import pandas as pd;
from pyproj import Transformer
from pathlib import Path




#import inputs
input = biab_inputs()
longitude = input['longitude']
latitude = input['latitude']
coordinates = input['coordinates']
input_crs = input['crs_input']
output_crs = input['crs_output']

#export output crs
biab_output("crs_output", output_crs)

if coordinates is not None and (longitude is not None or latitude is not None):
    # Mixed input: both coordinate table and long or lat are provided
    #raise ValueError("Error: Mixed input provided. Please provide either longitude and latitude or coordinate table.")
    biab_error_stop("Error: Mixed input provided. Please provide either longitude and latitude or coordinate table.")

if (longitude is None or latitude is None) and coordinates is None:
    # No coordinates table and not complete lat and long
    #raise ValueError("Error: No coordinates provided. Please provide either longitude and latitude or coordinates.")
    biab_error_stop("Error: No coordinates provided. Please provide either longitude and latitude or coordinates.")

if(input_crs is None or output_crs is None):
    # No input or output CRS provided
    biab_error_stop("Error: No input or output CRS provided. Please provide both input and output CRS.")


#Function to transform coordinates

def transform_crs(longitude, latitude, input_crs, output_crs):
    transformer = Transformer.from_crs(input_crs, output_crs, always_xy=True)
    # Transform arrays of longitude and latitude
    longitude_trans, latitude_trans = transformer.transform(longitude, latitude)
    return(longitude_trans, latitude_trans)


###RUN FUNCTIONS###

#For single point

if coordinates is None:
    print("converting single point coordinates...")
    #Transform the coordinates
    longitude_trans, latitude_trans = transform_crs(longitude, latitude, input_crs, output_crs)
    #create output
    biab_output("longitude_trans", longitude_trans)
    biab_output("latitude_trans", latitude_trans)
    biab_output("coordinates_trans", None)


##For coordinates table
if coordinates is not None:
    print("transforming coordinates in coordinate table...")
    #read the coordinates table
    coordinates_df = pd.read_csv(coordinates)
    #Check if the four first columns are the right ones
    colnames = list(coordinates_df.columns[0:3])
    if colnames != ['id','latitude', 'longitude']:
        biab_error_stop("Error: The first three columns must be named 'id', 'latitude', 'longitude'.")
    #Transform the coordinates
    longitude= coordinates_df['longitude'].to_numpy()
    latitude = coordinates_df['latitude'].to_numpy()
    longitude_trans, latitude_trans = transform_crs(longitude, latitude, input_crs, output_crs)
    coordinates_df["latitude_"+output_crs] = latitude_trans
    coordinates_df["longitude_"+output_crs] = longitude_trans
    #export the transformed coordinates
    #output_crs=os.path.join(os.getcwd(),"output", "biodiversity_dl_biab","transform_crs", "transformed_coordinates.csv")
    #output_crs=os.path.join("/","output", "biodiversity_dl_biab","transform_crs", "transformed_coordinates.csv")
    output_crs=os.path.join(output_folder, "transformed_coordinates.csv")
    coordinates_df.to_csv(output_crs, index=False)
    biab_output("coordinates_trans", output_crs)
    biab_output("longitude_trans", None)#Adding this lines to avoid error on the pipeline
    biab_output("latitude_trans", None)#Adding this lines to avoid error on the pipeline

