#load packages
library(codetools)
library(rjson)

print("packages loaded succesfully")

###SET INPUTS###
input <- biab_inputs()
print("Inputs: ")
print(input)
data_sources_config <- input$config_csv
choices <- input$user_choices

###HARD CODED INPUTS

#choices <- c("meteorology", "land cover", "sentinel2", "elevation")
#data_sources <- read.csv("userdata/datasources_config.csv",stringsAsFactors = FALSE)


#Import data sources config csv
data_sources <- read.csv(data_sources_config,stringsAsFactors = FALSE)
#subset local
local <- data_sources[data_sources[,1]=="local",c(2,3)]
#subset remote
remote <- data_sources[data_sources[,1]=="remote",c(2,3)]

#Select the names of tiff files or collections for the selected data sources
selection_local<-local[local[,1]%in%choices, 2]
selection_remote<-remote[remote[,1]%in%choices,2]
print(selection_local)
print(selection_remote)

# selection_remote_export<- c()
# for (remote in selection_remote){
# selection_remote_export <- cbind(selection_remote_export, remote)
# }


biab_output("selection_local", selection_local)
biab_output("selection_remote", selection_remote)

#I make a list of the sources that are local or remote
#The local sources will be sent to the crop raster script
#The remote sources will be sent to the download script
