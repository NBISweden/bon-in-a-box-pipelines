import sys, json;
import os;
import pandas as pd;
import ee;
import geemap;


#import inputs
input = biab_inputs()
#select inputs

coordinates = input['coordinates']
collections = input['collections']
gee_account = input['gee_account']
gee_token_path = input['gee_token_path']
crs = input['output_crs']
input_crs = input['input_crs']
image_size = input['image_size']
meters_per_pixel = input['m_pixel']

#Check that the collection input is a list. Otherwise it will iterate by character
if not isinstance(collections, list):
    collections = [collections]  
    print("Single collection converted to a list for gee download")


####HARD CODE TEST#####
# longitude = 648029
# latitude = 6639925
# # window = 0.1       #  window size (in degrees)
# date1 = "2020-06-01"  #start date
# date2 = "2020-06-30"  #  end date
# # date1 = pd.to_datetime(date1)
# # date2 = pd.to_datetime(date2)

# gee_token_path=os.path.join(os.getcwd(),"userdata", "gee_keys","cnn-sdm-6a23e8ab235d.json")
# gee_account = "ma-service@cnn-sdm.iam.gserviceaccount.com"
# collections = ['COPERNICUS/S2_SR_HARMONIZED']
# # #coordinates = os.path.join(os.getcwd(),"userdata", "sample_example_short.csv")
#coordinates = os.path.join(os.getcwd(),"userdata", "sample_example_gee.csv")
# crs = 'EPSG:3006'


#set hard coded values FOR NOW

offset = (image_size / 2) * meters_per_pixel


if coordinates is  None:
    biab_error_stop("Error: No coordinates table provided.")


###DEFINE FUNCTIONS###

#initialize the Earth Engine API
def GEE_initialize(gee_account,json_path):
    # service_account = 'adrba-198@ee-adrianbaggstrom.iam.gserviceaccount.com'
    # json_path = '/Users/toban562/Desktop/ee-adrianbaggstrom-acc7ced721ce.json'
    credentials = ee.ServiceAccountCredentials(gee_account, json_path)
    ee.Initialize(credentials)

#Function to get the spatial extent
def get_window_coordinates(target_point, offset=150):
    coordinates = [target_point[0]-offset, target_point[1]-offset, target_point[0]+offset, target_point[1]+offset]
    return coordinates

#Function to download sentinel 2 data
def get_satellite_data_gee(start_date,end_date,gee_geometry,crs,gee_output_file,meters_per_pixel,target_channel='NDVI'):
    def maskS2clouds(image):
        # This function is copied from GEE examples library and converted from JavaScript to Python
        # Function to mask clouds using the Sentinel-2 QA band.
        qa = image.select('QA60')
        # Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = 1 << 10
        cirrusBitMask = 1 << 11
        # Both flags should be set to zero, indicating clear conditions.
        mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
        # Return the masked and scaled data, without the QA bands.
        return image.updateMask(mask).divide(10000).select("B.*").copyProperties(image, ["system:time_start"])
    s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterDate(start_date, end_date) \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70)) \
        .map(maskS2clouds)
    s2_collection = s2_collection.mean()
    if target_channel == 'NDVI':
        target_channel = s2_collection.normalizedDifference(['B8', 'B4']).rename("NDVI")
    else:
        target_channel = s2_collection.select(target_channel)
    geemap.ee_export_image(target_channel, gee_output_file, scale=meters_per_pixel, region=gee_geometry, crs=crs)

#Function to download human influence index data

def get_hii_data_gee(feat,gee_geometry,crs,gee_output_file,meters_per_pixel):
    hii_coll = ee.ImageCollection("projects/HII/v1/" + feat).filter(ee.Filter.calendarRange(2020, 2020, 'year'))
    hii_image = hii_coll.limit(1, 'system:time_start', False).first()
    # img_collection = ee.ImageCollection.fromImages([hii_image])
    # stacked_image = img_collection.toBands()
    geemap.ee_export_image(hii_image, gee_output_file, scale=meters_per_pixel, region=gee_geometry, crs=crs)

#Function to download elevation data

def get_elevation_data_gee(gee_geometry,crs,gee_output_file,meters_per_pixel,):
    dem_image = ee.Image('projects/sat-io/open-datasets/ASTER/GDEM')
    geemap.ee_export_image(dem_image, gee_output_file, scale=meters_per_pixel, region=gee_geometry, crs=crs)

#Function to download climate data

def get_climate_data_bioclim_gee(target_channel,gee_geometry,crs,gee_output_file,meters_per_pixel):
    bioclim_image = ee.Image('WORLDCLIM/V1/BIO')
    output_channel = bioclim_image.select(target_channel)
    geemap.ee_export_image(output_channel, gee_output_file, scale=meters_per_pixel, region=gee_geometry, crs=crs)


#Function to select the right collection and export the data

def download_collection(collection):
    if collection == 'COPERNICUS/S2_SR_HARMONIZED':
        get_satellite_data_gee(date1, date2, my_region, crs, output_gee, meters_per_pixel)
    elif collection.startswith('projects/HII/v1/'):
        get_hii_data_gee('hii', my_region, crs, output_gee, meters_per_pixel)
    elif collection.startswith('projects/sat-io/open-datasets/ASTER/GDEM'):
        get_elevation_data_gee(my_region, crs, output_gee, meters_per_pixel)
    elif collection.startswith('WORLDCLIM/V1/BIO'):
        get_climate_data_bioclim_gee('bio01', my_region, crs, output_gee, meters_per_pixel)

###AUTHENTICATE###

print("Authenticating to Google earth Engine...")
GEE_initialize(gee_account, gee_token_path)

###RUN FUNCTIONS###
   

##For coordinates table
if coordinates is not None:
    print("Retrieving data for coordinates table...")
    #read the coordinates table
    coordinates_df = pd.read_csv(coordinates)
    #Check if the four first columns are the right ones
    colnames = list(coordinates_df.columns[0:5])
    if colnames != ['id','latitude', 'longitude', 'start_date', 'end_date']:
        biab_error_stop("Error: The first five columns must be named 'id', 'latitude', 'longitude', 'start_date' and 'end_date'.")
    #Check if the coordinates have been transformed previously
    lat_transformed = "latitude_" + input_crs.replace(":", "_")
    lon_transformed = "longitude_" + input_crs.replace(":", "_")
    coordinates_df.columns = coordinates_df.columns.str.replace(":", "_")
    #Create the empty list for the output paths
    file_paths = []
    print("created file path")
    #loop through the collections
    for collection in collections:
        print("Processing collection: " + collection)
        # Loop through the coordinates table
        for row in coordinates_df.itertuples():
            sample_id = row[1]

            # Re-incorporating the check for pre-transformed columns
            # This logic now handles two scenarios
            if lat_transformed in coordinates_df.columns and lon_transformed in coordinates_df.columns:
                # If the coordinates have been transformed, use the transformed coordinates
                lon = getattr(row, lon_transformed)
                lat = getattr(row, lat_transformed)
                # The CRS for this point is the user-specified input_crs
                point_crs = input_crs
            else:
                # If the coordinates have not been transformed, use the original coordinates
                lon = row.longitude
                lat = row.latitude
                # The CRS for this point is the user-specified input_crs
                point_crs = input_crs

            # Create an ee.Geometry.Point from the input coordinates and the correct CRS
            point_geom = ee.Geometry.Point([lon, lat], point_crs)

            # Get the projected coordinates of the point in the desired output CRS
            # If the point_crs was already the output_crs, this will have no effect.
            projected_point = point_geom.transform(crs).getInfo()['coordinates']

            # Use the newly projected coordinates to define the window
            target_point = [projected_point[0], projected_point[1]]


            date1 = row.start_date
            date2 = row.end_date
            #get target point
            target_point = [lon, lat]
            #spatial  extent
            my_extent = get_window_coordinates(target_point, offset)
            my_region = ee.Geometry.Rectangle(my_extent, proj=crs, geodesic=False)
            #file name and path
            file_name = str(sample_id) + "_" +collection.replace("/", "_") + ".tif"
            #output_gee=os.path.join("/","output", "biodiversity_dl_biab","get_GEE_data", file_name)
            output_gee=os.path.join(output_folder, file_name)
            #output_gee=os.path.join(os.getcwd(),"output", "biodiversity_dl_biab","get_GEE_data", file_name)
            #Download the remote sensing data
            download_collection(collection)
            # Append the output path to the list
            file_paths.append(output_gee)
    biab_output("gee_raster", file_paths)

