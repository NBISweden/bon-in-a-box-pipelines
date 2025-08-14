import openeo;
import sys, json;
import os;
import pandas as pd;


#check that openeo it is properly installed
print("dependencies loaded succesfully")
print("openEO version: ", openeo.__version__)


#import inputs
input = biab_inputs()
#select inputs
longitude = input['longitude']
latitude = input['latitude']
coordinates = input['coordinates']
window = input['window']
date1= input['date1']
date2= input['date2']
collections = input['collections']

if coordinates is not None and (longitude is not None or latitude is not None):
    # Mixed input: both coordinate table and long or lat are provided
    #raise ValueError("Error: Mixed input provided. Please provide either longitude and latitude or coordinate table.")
    biab_error_stop("Error: Mixed input provided. Please provide either longitude and latitude or coordinate table.")

if (longitude is None or latitude is None) and coordinates is None:
    # No coordinates table and not complete lat and long
    #raise ValueError("Error: No coordinates provided. Please provide either longitude and latitude or coordinates.")
    biab_error_stop("Error: No coordinates provided. Please provide either longitude and latitude or coordinates.")



#print("Retrieving data from openEO for a site with ", longitude, " longitude and ", latitude, " latitude", " and a window of ", window, " degrees", "between ", date1, " and ", date2)
#print(longitude)
#print(latitude)
#print(window)

####HARD CODE TEST#####
# longitude = 17.64  # Example longitude
# latitude = 59.86   # Example latitude
# window = 0.1       # Example window size (in degrees)
# date1 = "2020-08-01"  # Example start date
# date2 = "2020-08-02"  # Example end date
# refresh_token_file_path=os.path.join(os.getcwd(),"userdata", "refresh-tokens.json" )
# collections = ["SENTINEL2_L1C"]
# #coordinates = os.path.join(os.getcwd(),"userdata", "sample_example_short.csv")
# coordinates = None

###DEFINE FUNCTIONS###

#Function to retrieve the specific refresh token from the JSON file
def get_refresh_token(file_path):
      with open(file_path, 'r') as f:
        token_data = json.load(f)
        my_token = token_data.get('https://aai.egi.eu/auth/realms/egi', {}) \
            .get('vito-default-client', {}) \
            .get('refresh_token')
        if my_token is None:
            print("Refresh token not found in the JSON file.")
            sys.exit(1)
        return my_token

#Function to set the spatial extent to crop the data

def get_spatial_extent(longitude, latitude, window):
    """Calculate the bounding box for the given longitude, latitude, and window size."""
    # Calculate the bounds (min/max longitude/latitude for the window)
    west = longitude - window / 2
    east = longitude + window / 2
    south = latitude - window / 2
    north = latitude + window / 2
    return {"west": west, "south": south, "east": east, "north": north}

# Function to crop and download the data

def download_data(my_connection, collection_id, my_spatial_extent, my_temp_extent, my_band, output_path):
    sentinel2_cube=my_connection.load_collection(collection_id,
                                          spatial_extent=my_spatial_extent,
                                          temporal_extent=my_temp_extent, 
                                          bands=my_band)
    #remove the temporal extent
    sentinel2_cube_notime = sentinel2_cube.max_time()
   #Download data
    print("Downloading data for ", collection_id, "band", my_band)
    #1. Specify the path for the download
    #output_openeo=os.path.join("/","output", "biodiversity_dl_biab","get_openEO_data","sentinel2_cube.tif")
    #2. Specify the output of the download for the biab and download the data
    sentinel2_cube_notime.download(output_path)
   

###AUTHENTICATE###

print("Authenticating to openEO backend")
# Replace with the actual path to your JSON file containing the refresh token
refresh_token_file_path = input['refresh_token_file']  # e.g., "path/to/your/refresh_token.json"
#print(refresh_token_file_path)
my_token = get_refresh_token(refresh_token_file_path)
#Connect to the VITO backend
connection = openeo.connect("openeo.vito.be")
connection.authenticate_oidc_refresh_token(refresh_token=my_token)
#print("Authenticated with refresh token")

###RUN FUNCTIONS###

#For single point
if coordinates is None:
    #Create the empty list for the output paths
    file_paths = []
    # Create the spatial extent
    my_spatial_extent = get_spatial_extent(longitude, latitude, window)
    # Create the temporal extent
    my_temp_extent = [date1, date2]
    # set output path
    for collection in collections:
        #Select the first band of the collection ONLY FOR TESTING, BECAUSE WE HAVE LIMITED CREDITS
        all_bands=connection.describe_collection(collection)
        my_band = all_bands["cube:dimensions"]["bands"]["values"][3]
        file_name = collection + ".tif"
        output_openeo=os.path.join("/","output", "biodiversity_dl_biab","get_openEO_data", file_name)
        # Download the data
        download_data(connection, collection, my_spatial_extent, my_temp_extent, my_band, output_openeo)
        # Append the output path to the list
        file_paths.append(output_openeo)
    biab_output("openeo_raster", file_paths)


##For coordinates table
if coordinates is not None:
    #read the coordinates table
    coordinates_df = pd.read_csv(coordinates)
    #Check if the four first columns are the right ones
    colnames = list(coordinates_df.columns[0:5])
    if colnames != ['id','latitude', 'longitude', 'start_date', 'end_date']:
        biab_error_stop("Error: The first five columns must be named 'id', 'latitude', 'longitude', 'start_date' and 'end_date'.")
    #Create the empty list for the output paths
    file_paths = []
    #loop through the collections
    for collection in collections:
        #Select the first band of the collection ONLY FOR TESTING, BECAUSE WE HAVE LIMITED CREDITS
        all_bands=connection.describe_collection(collection)
        my_band = all_bands["cube:dimensions"]["bands"]["values"][3]
        # Loop through the coordinates table
        for row in coordinates_df.itertuples():
            sample_id = row[1]
            lat = row[2]
            lon = row[3]
            start_date = row[4]
            end_date = row[5]
            #spatial and temporal extent
            my_spatial_extent = get_spatial_extent(lon, lat, window)
            my_temp_extent = [start_date, end_date]
            #file name
            file_name = str(sample_id) + "_" +collection + ".tif"
            output_openeo=os.path.join("/","output", "biodiversity_dl_biab","get_openEO_data", file_name)
            #Download the data
            download_data(connection, collection, my_spatial_extent, my_temp_extent, my_band,output_openeo)
            # Append the output path to the list
            file_paths.append(output_openeo)
    
    biab_output("openeo_raster", file_paths)




# #Display the image
# import os
# import imageio.v2 as imageio
# import matplotlib.pyplot as plt

# parent_dir=os.getcwd()
# parent_dir=os.path.join(parent_dir,"output", "biodiversity_dl_biab","get_openEO_data")
# def convert_image(parent_dir, file):
#     tif_file= imageio.imread(os.path.join(parent_dir, file))
#     plt.imshow(tif_file)
#     plt.show()
#     return tif_file

# uppsala_wgs=convert_image(parent_dir, "STO_1_SENTINEL2_L1C.tif" )
# #Get shape of the image
# tif_file= imageio.imread(os.path.join(parent_dir,"sentinel2_cube.tif"))
# tif_file.shape