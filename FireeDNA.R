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
require(dplyr)

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
if(length(list.files(pattern="RawElevation.tif",full.names=TRUE))!=1){
  elevation <- elevatr::get_elev_raster(StudyArea, z = 12)
  elevation[is.infinite(elevation)] <- -9999
  raster::NAvalue(elevation) <- -9999
  raster::writeRaster(elevation, filename="RawElevation.tif", options=c("COMPRESS=LZW", "TFW=YES"), format="GTiff", overwrite=TRUE)
}

#Read in eDNA sample coordinates.
SamplePoints <- read.table("Woolsey_metadata_for_phyloseq_May20.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#Only keep sediment samples.
SamplePoints <- SamplePoints[!is.na(SamplePoints$Latitude) & SamplePoints$Sample_type=="Sediment",]
#Read in alpha diversity by sample data.
SampleDiversity <- read.table("AllPrimersAlphaDiversity.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#Create numerical variable for field season.
SampleDiversity$FieldSeasonNum <-sapply(as.character(SampleDiversity$Field_Season), switch, "Pre_burn_Aug2018" = 1, "Post_burn_NovDec2018" = 2, "Post_burn_Feb2019" = 3, "Post_burn_Jun2019" = 4, USE.NAMES = F)
#Only keep sediment samples
SampleDiversity <- SampleDiversity[SampleDiversity$sum.taxonomy %in% SamplePoints$sum.taxonomy,]
#Merge in diversity metrics.
SamplePoints <- dplyr::left_join(SamplePoints,SampleDiversity[,c("sum.taxonomy","Chao1","se.chao1","Shannon","Simpson","primer")])
#Remove missing primer data.
SamplePoints <- SamplePoints[!is.na(SamplePoints$primer),]
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

#Extract local environmental values at sampling locations.  Find the mean within a 2km buffer.
LocalBuffer <- sf::st_buffer(sf::st_as_sf(SampleCoordinates),2000)
LocalEnvironment <- exactextractr::exact_extract(env.data, LocalBuffer$geometry,"mean",stack_apply=T)
LocalEnvironment <- cbind(SamplePoints,LocalEnvironment)
names(LocalEnvironment) <- gsub(x = names(LocalEnvironment), pattern = "\\mean.", replacement = "")  
LocalEnvironment <- LocalEnvironment[!is.na(LocalEnvironment$primer),]

#Example map of environmental raster with sampling points and buffer.
plot(raster("RFAL20161000m.tif"))
plot(SampleCoordinates$geometry,add=T)
plot(LocalBuffer$geometry,add=T)

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
#Remove rows with missing primer data.
watershedsWithMetadata <- watershedsWithMetadata[!is.na(watershedsWithMetadata$primer),]
#Define variables to use for modeling
ModelVars <- c("AreaSqKm.y","aspect1000m.Mean","aspect1000m.SD","BurnSeverity1000m.Mean","BurnSeverity1000m.SD","DEF1000m.Mean","DEF1000m.SD","Elevation1000m.Mean","Elevation1000m.SD","PET1000m.Mean","PET1000m.SD","RFAL20161000m.Mean","RFAL20161000m.SD","slope1000m.Mean","slope1000m.SD","TDEW1000m.Mean","TDEW1000m.SD","TMAX1000m.Mean","TMAX1000m.SD","TMIN1000m.Mean","TMIN1000m.SD")
DiversityVars <- c("Chao1","Shannon","Simpson")
CategoricalVars <- c("Lagoon_Location","Field_Season","HUC12.y")
watershedsSubset <- watershedsWithMetadata[,c(ModelVars,CategoricalVars,DiversityVars,"se.chao1","primer")]
#Convert character columns to factors.
watershedsSubset[,CategoricalVars] <- lapply(watershedsSubset[,CategoricalVars],factor)
#Create numerical variable for field season.
watershedsSubset$FieldSeasonNum <-sapply(as.character(watershedsSubset$Field_Season), switch, "Pre_burn_Aug2018" = 1, "Post_burn_NovDec2018" = 2, "Post_burn_Feb2019" = 3, "Post_burn_Jun2019" = 4, USE.NAMES = F)

#Calculate F-statistics, and their significance, for linear models of diversity derived from watershed summary statistics.
WatershedTestSummary <- data.frame()
for(DiversityVar in DiversityVars){
  for(Primer in unique(watershedsSubset$primer)){
    for(Var in c("Lagoon_Location","HUC12.y",ModelVars)){
      ModelInput <- watershedsSubset[watershedsSubset$primer==Primer,c(Var,DiversityVar)]
      test <- anova(lm(ModelInput[,2]~ModelInput[,1]))
      TestOutput <- as.data.frame(t(unlist(na.omit(c(Primer,Var,DiversityVar,test$Df[1],test$Df[2],test$`F value`,test$`Pr(>F)`)))))
      colnames(TestOutput) <- c("Primer","Variable","DiversityMetric","Df","n","F","p")
      WatershedTestSummary <- rbind(WatershedTestSummary,TestOutput)
    }
  }
}

#Calculate F-statistics, and their significance, for linear models of diversity derived from local summary statistics
ModelVars <- c("aspect1000m","BurnSeverity1000m","DEF1000m","Elevation1000m","PET1000m","RFAL20161000m","slope1000m","TDEW1000m","TMAX1000m","TMIN1000m")
CategoricalVars <- c("Lagoon_Location","Field_Season")
DiversityVars <- c("Chao1","Shannon","Simpson")
LocalTestSummary <- data.frame()
for(DiversityVar in DiversityVars){
  for(Primer in unique(LocalEnvironment$primer)){
    for(Var in c(CategoricalVars,ModelVars)){
      ModelInput <- LocalEnvironment[LocalEnvironment$primer==Primer,c(Var,DiversityVar)]
      test <- anova(lm(ModelInput[,2]~ModelInput[,1]))
      TestOutput <- as.data.frame(t(unlist(na.omit(c(Primer,Var,DiversityVar,test$Df[1],test$Df[2],test$`F value`,test$`Pr(>F)`)))))
      colnames(TestOutput) <- c("Primer","Variable","DiversityMetric","Df","n","F","p")
      LocalTestSummary <- rbind(LocalTestSummary,TestOutput)
    }
  }
}

#Combine local and watershed model summary statistics.
FullSummary <- rbind(LocalTestSummary,WatershedTestSummary)
FullSummary <- FullSummary[!duplicated(FullSummary),]
FullSummary <- FullSummary[order(FullSummary$Primer,FullSummary$Variable,FullSummary$DiversityMetric,FullSummary$p),]
FullSummaryFiltered <- FullSummary[FullSummary$p <= 0.05,]

#Plot diversity metric distributions, split by primer, against sampling time.
SampleDiversity %>%
  ggplot(aes(x=FieldSeasonNum,y=Simpson,group=factor(FieldSeasonNum,1:4)))+
  guides(fill=guide_legend(title="Field season"))+
  geom_boxplot(aes(fill=factor(FieldSeasonNum)))+
  facet_grid(primer~.)

#Test how significant diversity distributions, split by primer, differ by time.
for(Primer in unique(SampleDiversity$primer)){
  tmp <- kruskal.test(Simpson~FieldSeasonNum,data=SampleDiversity[SampleDiversity$primer==Primer,])$p.value
  print(paste(Primer,tmp))
}

#Generate maps with watershed level summaries of environmental layers, along with eDNA sample points.
for(env.name in env.names){
  if(!(env.name %in% c("SoilProperties1000m"))){
    print(env.name)
    EnvVar <- paste(env.name,".Mean",sep="")
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_colorbar(paste(gsub("1000m","",env.name),"\nMean")))+theme(legend.position = "bottom")
    TestMap <- TestMap+geom_sf(data=SampleCoordinates,color="red")
    TestMap
    ggsave(TestMap,file=paste(gsub("1000m","",env.name),"WatershedMeans.png",sep=""),width=7,height=7,units="in")
    EnvVar <- paste(env.name,".SD",sep="")
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_colorbar(paste(gsub("1000m","",env.name),"\nStandard deviation")))+theme(legend.position = "bottom")
    TestMap <- TestMap+geom_sf(data=SampleCoordinates,color="red")
    TestMap
    ggsave(TestMap,file=paste(gsub("1000m","",env.name),"WatershedSDs.png",sep=""),width=7,height=7,units="in")
  }else{
    print(env.name)
    EnvVar <- env.name
    TestMap <- ggplot()+geom_sf(data=watershedSummary,aes(fill=!!sym(EnvVar)))+scale_fill_viridis_c()+guides(fill=guide_colorbar(paste(gsub("1000m","",env.name),"\nMode")))+theme(legend.position = "bottom")
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

