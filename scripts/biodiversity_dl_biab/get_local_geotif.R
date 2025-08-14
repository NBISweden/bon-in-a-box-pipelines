#load packages
library(codetools)
library(rjson)
library(terra)
print("packages loaded succesfully")

###SET INPUTS###
input <- biab_inputs()

global_tif_path<-input$global_tifs
coordinates<-input$coordinates
window<-input$window


###SET ERROR MESSAGES FOR INCONSISTENT OR MISSING INPUTS###
#if no coordinates are provided we will print an error and stop


if (is.null(coordinates)) {
    #If no coordinates are provided, print an error and stop
    stop("Error: No coordinates provided. Please provide coordinate table in csv format.")
}

## SET OUTPUT FOR WHEN NO LOCAL RASTER IS PROVIDED

if (length(global_tif_path) == 0) {
    print("No local tiff files provided")
    biab_output("output_text", "No local tiff files: No local environmental source was provided")
}

###DEFINE FUNCTIONS###

#Function to extract a window of raster data based on coordinates and window size
extract_window<-function(raster_file, longitude, latitude, window){
window_extent <- ext(longitude - window/2, longitude + window/2, 
                        latitude - window/2, latitude + window/2)
cropped_raster <- crop(raster_file, window_extent)
#plot(cropped_raster)
return(cropped_raster)
}

#Function to export the cropped raster 
export_raster<-function(raster_file, output_file){
print("exporting cropped raster")
#print(output_cropped_raster)
writeRaster(raster_file, output_file, filetype="GTiff", overwrite=TRUE)
}

#Function to extract the name of the cropped raster
raster_name <- function(raster_path){
    # Get the file name without the path
    file_name <- basename(raster_path)
    # If I want to remove the file extension and the directory path
    # file_name_no_ext <- tools::file_path_sans_ext(file_name)
    return(file_name)
}


###RUN FUNCTIONS FOR COORDINATES TABLE

if (!is.null(coordinates) && length(global_tif_path) != 0) {
    #If coordinates are provided, use them
    print("Using coordinate file")
    coordinates<-read.csv(coordinates)
    head(coordinates)
    file_paths <- c()
    #iterate over the provided tif paths
    for (path_id in global_tif_path){
        global_tif <-rast(path_id)
        name <- raster_name(path_id)
    #Iterate over the provided coordinates
        for (i in 1:nrow(coordinates)) {
            tiff_name<-paste0("cropped_",coordinates[i,"id"],"_", name)
            print(tiff_name)
            #Because the local TIFFs are always in WGS84, we will use those crs (EPSG:4326)
            if ("latitude_EPSG_4326" %in% colnames(coordinates) && "longitude_EPSG_4326" %in% colnames(coordinates)) {
                print("Using latitude and longitude in EPSG:4326")
                latitude <- as.numeric(coordinates[i,"latitude_EPSG_4326"])
                longitude <- as.numeric(coordinates[i,"longitude_EPSG_4326"])
            } else {
                print("Using latitude and longitude without transformation")
                latitude <- as.numeric(coordinates[i,"latitude"])
                longitude <- as.numeric(coordinates[i,"longitude"])
            }
            cropped_raster<-extract_window(global_tif, longitude, latitude, window)
            #Export the cropped raster
            output_file<-file.path(outputFolder, tiff_name)
            export_raster(cropped_raster, output_file)
            file_paths <- cbind(file_paths, output_file)
            #biab_output("cropped_raster", output_file)
        }
    }
#Export the cropped rasters
biab_output("cropped_raster", file_paths)
}

