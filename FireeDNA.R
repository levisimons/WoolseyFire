rm(list=ls())
require(raster)
require(rgdal)
require(sp)
require(elevatr)
require(rasterVis)
require(elsa)
require(ggplot2)
require(viridis)
require(sf)
require(nngeo)
require(reshape2)

wd <- "~/Desktop/CUBoulder/"
#wd <- "/home1/alsimons/CUBoulder"
setwd(wd)

#Run once.
#Generate a plot of the Moran's I value versus distance for environmental map layers.
#This is used to determine the resampling scale for the map layers.
if(length(list.files(pattern="MapMoran.txt",full.names=TRUE))==0){
  length(list.files(pattern="MapMoran.txt",full.names=TRUE))
  MapRasters <- list.files(pattern="(.*?).tif")
  MapMoran <- data.frame()
  for(MapRaster in MapRasters){
    Map <- raster(MapRaster)
    print(res(Map)[1])
    tmp <- correlogram(Map,res(Map)[1],7500)
    tmp <- tmp@correlogram
    tmp$Layer <- gsub(".tif","",MapRaster)
    MapMoran <- rbind(MapMoran,tmp)
  }
  write.table(MapMoran,"MapMoran.txt",quote=FALSE,sep="\t",row.names = FALSE)
}
#
MapMoran <- read.table("MapMoran.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
ggplot(MapMoran,aes(distance,moran,col=Layer))+geom_line()

#Read in study area boundary.
StudyArea <- sf::st_read("SMMNRA_study_area.shp")

#Create elevation raster.  This will get reprojected and clipped to the study area.
#Only run this once.
if(length(list.files(pattern="RawElevation.tif",full.names=TRUE))==0){
  elevation <- elevatr::get_elev_raster(StudyArea, z = 12)
  elevation[is.infinite(elevation)] <- -9999
  raster::NAvalue(elevation) <- -9999
  raster::writeRaster(elevation, filename="RawElevation.tif", options=c("COMPRESS=LZW", "TFW=YES"), format="GTiff", overwrite=TRUE)
}

#Read in eDNA sample coordinates.
SamplePoints <- read.table("Woolsey_metadata_for_phyloseq_May20.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#Only keep sediment samples.
SamplePoints <- SamplePoints[!is.na(SamplePoints$Latitude) & SamplePoints$Sample_type=="Sediment",]
#Convert to spatial data frame in the same project CRS.
SampleCoordinates <- sf::st_as_sf(SamplePoints,coords=c("Longitude","Latitude"),crs = st_crs(4326))
#Reproject the species observations into the same CRS as the western Santa Monica mountains.
SampleCoordinates <- sf::st_transform(SampleCoordinates,crs=st_crs(32611))

#Read in HUC 12 watershed boundaries.
watersheds <- sf::st_read("watersheds.shp")

#Get resampled raster lists.
env.files <- list.files(pattern="1000m.tif$",full.names=TRUE)
#Get raster layer names.
env.names <- gsub("./|.tif","",env.files)
#Stack environmental layers
env.data <- stack(c(env.files))

#Get summary statistics, per HUC12 watershed, for environmental layers used in this project.
#Get mean raster value per HUC 12 watershed.
watershedMeans <- sf::st_as_sf(extract(env.data, as_Spatial(watersheds), fun=mean, na.rm=TRUE, df=TRUE, sp=TRUE))
#Get standard deviation on the raster value per HUC 12 watershed.
watershedSDs <- sf::st_as_sf(extract(env.data, as_Spatial(watersheds), fun=sd, na.rm=TRUE, df=TRUE, sp=TRUE))
#Get the mode raster value per HUC 12 watershed.
watershedModes <- sf::st_as_sf(extract(env.data, as_Spatial(watersheds), fun=modal, na.rm=TRUE, df=TRUE, sp=TRUE))
#Create summary statistics layer for rasters values by HUC 12 watershed.
tmp <- sf::st_join(watershedMeans,watershedSDs,suffix=c(".Mean",".SD"),join=st_nn,k=1,maxdist=1000)
watershedSummary <- sf::st_join(tmp,watershedModes,join=st_nn,k=1,maxdist=1000)
tmp <- NULL
#Join sampling coordinates to watershed.
watershedSamplingIntersections <- sf::st_join(SampleCoordinates,watersheds,suffix=c("",".y"),join=st_nn,k=1,maxdist=1000)
#Join watersheds. summary statistics, and sampling data into a single spatial dataframe.
watershedsWithMetadata <- sf::st_join(watershedSamplingIntersections,watershedSummary,suffix=c(".y",""),join=st_nn,k=1,maxdist=1000)
#Get coordinates into a data frame.
xy <- as.data.frame(sf::st_coordinates(watershedsWithMetadata))
colnames(xy) <- c("longitude","latitude")
#Coerce merged spatial data frame to data frame.
watershedsWithMetadata <- as.data.frame(watershedsWithMetadata)
#Add into spatial information as latitude and longitude columns.
watershedsWithMetadata$geometry <- NULL
watershedsWithMetadata <- cbind(xy,watershedsWithMetadata)
#Remove duplicated columns
watershedsWithMetadata <- watershedsWithMetadata[!duplicated(as.list(watershedsWithMetadata))]

#Generate maps with watershed level summaries of environmental layers, along with eDNA sample points.
for(env.name in env.names){
  if(!(env.name %in% c("SoilProperties1000m"))){
    print(env.name)
    EnvVar <- paste(env.name,".Mean",sep="")
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_legend(paste(gsub("1000m","",env.name),"Mean")))
    TestMap <- TestMap+geom_sf(data=SampleCoordinates,color="red")
    TestMap
    ggsave(TestMap,file=paste(gsub("1000m","",env.name),"WatershedMeans.png",sep=""),width=7,height=7,units="in")
    EnvVar <- paste(env.name,".SD",sep="")
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_legend(paste(gsub("1000m","",env.name),"Mean")))
    TestMap <- TestMap+geom_sf(data=SampleCoordinates,color="red")
    TestMap
    ggsave(TestMap,file=paste(gsub("1000m","",env.name),"WatershedSDs.png",sep=""),width=7,height=7,units="in")
  }else{
    print(env.name)
    EnvVar <- env.name
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_legend(paste(gsub("1000m","",env.name),"Mode")))
    TestMap <- TestMap+geom_sf(data=SampleCoordinates,color="red")
    TestMap
    ggsave(TestMap,file=paste(gsub("1000m","",env.name),"WatershedModes.png",sep=""),width=7,height=7,units="in")
  }
}

#Plot Pearson correlation coefficients between environmental layers.
RasterCor <- raster::layerStats(env.data,'pearson',na.rm=T)
RasterCor <- as.matrix(RasterCor$`pearson correlation coefficient`)
RasterCor <- round(RasterCor,2)
RasterCor <- reshape2::melt(RasterCor)
ggplot(RasterCor,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)+
  scale_fill_viridis_c()+
  theme(axis.text.y=element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
