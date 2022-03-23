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
require(exactextractr)

wd <- "~/Desktop/CUBoulder/"
#wd <- "/home1/alsimons/CUBoulder"
setwd(wd)

#Only run this once.
#Create elevation raster.  This will get reprojected and clipped to the study area.
if(length(list.files(pattern="RawElevation.tif",full.names=TRUE))!=1){
  elevation <- elevatr::get_elev_raster(StudyArea, z = 12)
  elevation[is.infinite(elevation)] <- -9999
  raster::NAvalue(elevation) <- -9999
  raster::writeRaster(elevation, filename="RawElevation.tif", options=c("COMPRESS=LZW", "TFW=YES"), format="GTiff", overwrite=TRUE)
}

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
#Read in HUC 12 watershed boundaries.
watersheds <- sf::st_read("watersheds.shp")

#Get resampled raster lists.
env.files <- list.files(pattern="1000m.tif$",full.names=TRUE)
#Get raster layer names.
env.names <- gsub("./|.tif","",env.files)
#Stack environmental layers
env.data <- stack(c(env.files))

#Read in eDNA sample coordinates.
SamplePoints <- read.table("Woolsey_metadata_for_phyloseq_May20.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#Only keep sediment samples.
SamplePoints <- SamplePoints[!is.na(SamplePoints$Latitude) & SamplePoints$Sample_type=="Sediment",]
#Read in alpha diversity by sample data.
SampleDiversity <- read.table("AllPrimersAlphaDiversity.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE,quote = "\"", encoding = "UTF-8")
#Only keep sediment samples
SampleDiversity <- SampleDiversity[SampleDiversity$sum.taxonomy %in% SamplePoints$sum.taxonomy,]
#Create numerical variable for field season.
SampleDiversity$FieldSeasonNum <- as.numeric(unlist(sapply(as.character(SampleDiversity$Field_Season), switch, "Pre_burn_Aug2018" = 1, "Post_burn_NovDec2018" = 2, "Post_burn_Feb2019" = 3, "Post_burn_Jun2019" = 4, USE.NAMES = F)))
#Merge in diversity metrics.
SamplePoints <- dplyr::left_join(SamplePoints,SampleDiversity[,c("sum.taxonomy","Chao1","se.chao1","Shannon","Simpson","primer")])
#Remove missing primer data.
SamplePoints <- SamplePoints[!is.na(SamplePoints$primer),]
#Convert to spatial data frame in the same project CRS.
SampleCoordinates <- sf::st_as_sf(SamplePoints,coords=c("Longitude","Latitude"),crs = st_crs(4326))
#Reproject the species observations into the same CRS as the western Santa Monica mountains.
SampleCoordinates <- sf::st_transform(SampleCoordinates,crs=st_crs(32611))

#Extract local environmental values at sampling locations.  Find the mean within a 2km buffer.
LocalBuffer <- sf::st_buffer(sf::st_as_sf(SampleCoordinates),2000)
LocalEnvironment <- exactextractr::exact_extract(env.data, LocalBuffer$geometry,"mean",stack_apply=T)
LocalEnvironment <- cbind(SamplePoints,LocalEnvironment)
names(LocalEnvironment) <- gsub(x = names(LocalEnvironment), pattern = "\\mean.", replacement = "")  
LocalEnvironment <- LocalEnvironment[!is.na(LocalEnvironment$primer),]
#Extract the local mode of categorical variables, replace the mean values with it.
LocalEnvironment$CommunityType1000m <- exactextractr::exact_extract(env.data$CommunityType1000m, LocalBuffer$geometry,"mode",stack_apply=T)

#Join sampling coordinates to watershed.
watershedSamplingIntersections <- sf::st_join(SampleCoordinates,watersheds,suffix=c("",".y"),join=st_nn,k=1,maxdist=1000)

##Test if watershed environments are significantly different.
#Get all watersheds where sampling took place.
SamplingWatersheds <- watersheds[watersheds$HUC12 %in% watershedSamplingIntersections$HUC12,]
SamplingWatersheds <- SamplingWatersheds$HUC12
#Get all raster values per HUC12 watershed where sampling took place.
watershedValues <- data.frame()
for(SamplingWatershed in SamplingWatersheds){
  watershedValuesTmp <- exactextractr::exact_extract(env.data,watersheds[watersheds$HUC12==SamplingWatershed,"geometry"],default_value=0)
  watershedValuesTmp <- watershedValuesTmp[[1]]
  watershedValuesTmp$HUC12 <- SamplingWatershed
  watershedValuesTmp$coverage_fraction <- NULL
  watershedValues <- rbind(watershedValues,watershedValuesTmp)
}
#Generate boxplots of raster values based on sampled watersheds.
watershedValues <- reshape2::melt(watershedValues,id.var=c("HUC12"))
#Select environmental variable to plot
ggplot(data=watershedValues[watershedValues$variable==env.names[16],],aes(x=variable,y=value,fill=factor(HUC12)))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=T)+
  facet_grid(cols=vars(HUC12),labeller=label_both)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.position = "none")+
  labs(title="Distribution of mean annual minimum temperature by sampled watershed",y = "Degrees C")
#Map of environmental raster with sampling points and buffer.
plot(raster(env.files[16]))
title(main="1 km resolution map of mean annual minimum temperature",xlab="longitude",ylab="latitude")
plot(watersheds$geometry,add=T)
plot(SampleCoordinates$geometry,add=T)
plot(LocalBuffer$geometry,add=T)
#Test how significant watersheds differ by environmental metrics.
for(var in unique(watershedValues$variable)){
  tmp <- kruskal.test(value~HUC12,data=watershedValues[watershedValues$variable==var,])$p.value
  print(paste(var,tmp))
}
#Plot watersheds used in sampling.
ggplot(data=watersheds[watersheds$HUC12 %in% SamplingWatersheds,"HUC12"],aes(fill=factor(HUC12)))+
  labs(title="Sampled watersheds",x="longitude",y="latitude",fill='HUC 12\nwatersheds')+
  geom_sf()

#Get summary statistics, per HUC12 watershed, for environmental layers used in this project.
watershedSummary <- cbind(watersheds$HUC12,exactextractr::exact_extract(env.data,watersheds$geometry,fun=c("mean","stdev","mode"),default_value=0,stack_apply=T))
names(watershedSummary)[names(watershedSummary) == "watersheds$HUC12"] <- "HUC12"
#Merge summary statistics, per HUC12 watershed, for environmental layers used in this project
#along with local environmental measurements.
SampledWatersheds <- as.data.frame(watershedSamplingIntersections)
SampledWatersheds$geometry <- NULL
SampledWatersheds <- SampledWatersheds[,c("sum.taxonomy","HUC12","Chao1","Shannon","Simpson","primer")]
SampledWatershedsWithSummary <- dplyr::left_join(SampledWatersheds,watershedSummary,by=c("HUC12"))
SampledWatershedsWithSummary <- dplyr::left_join(SampledWatershedsWithSummary,LocalEnvironment[,c("sum.taxonomy",env.names)],by=c("sum.taxonomy"))
SampledWatershedsWithSummary <- dplyr::left_join(SampleDiversity[,c("sum.taxonomy","FieldSeasonNum")],SampledWatershedsWithSummary,by=c("sum.taxonomy"))

#Test how different alpha diversity values are compared to watershed summary statistics for static variables.
#Define diversity metrics.
DiversityVars <- c("Chao1","Shannon","Simpson")
#Select watershed summary variables.
ModelVars <- colnames(SampledWatershedsWithSummary)[!(colnames(SampledWatershedsWithSummary) %in% c("sum.taxonomy","FieldSeasonNum","HUC12","primer",DiversityVars))]
RemoveVars <- c("mean.CommunityType1000m","stdev.CommunityType1000m","mode.aspect1000m","mode.BurnSeverity1000m","mode.DEF1000m","mode.Elevation1000m","mode.Rainfall1Sample1000m","mode.Rainfall2Sample1000m","mode.Rainfall3Sample1000m","mode.Rainfall4Sample1000m","mode.RFAL20161000m","mode.slope1000m")
ModelVars <- ModelVars[!grepl(paste(RemoveVars,collapse="|"),ModelVars)]
ModelVars <- ModelVars[grepl(c("mean|stdev|mode"),ModelVars)]
StaticVars <- ModelVars[grepl(c("aspect|CommunityType|DEF|Elevation|RFAL2016|slope"),ModelVars)]
#Aggregate Kruskal-Wallis test results.
StaticWatershedTestSummary <- data.frame()
for(Primer in unique(SampledWatershedsWithSummary$primer)){
  watershedsSubset <- SampledWatershedsWithSummary[SampledWatershedsWithSummary$primer==Primer,]
  ModelInput <- watershedsSubset[!duplicated(watershedsSubset),]
  for(Var in StaticVars){
    for(DiversityVar in DiversityVars){
      StaticWatershedTest <- data.frame(matrix(nrow=1,ncol=6))
      colnames(StaticWatershedTest) <- c("Primer","DiversityVar","Var","ChiSquared","df","p")
      StaticWatershedTest[,"Primer"] <- Primer
      StaticWatershedTest[,"DiversityVar"] <- DiversityVar
      StaticWatershedTest[,"Var"] <- Var
      StaticWatershedTest[,"ChiSquared"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$statistic
      StaticWatershedTest[,"df"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$parameter
      StaticWatershedTest[,"p"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$p.value
      StaticWatershedTestSummary <- rbind(StaticWatershedTestSummary,StaticWatershedTest)
    }
  }
}
#Filter to only contain significant results.
StaticWatershedTestSummarySignificant <- StaticWatershedTestSummary[StaticWatershedTestSummary$p <= 0.05,]

#Test how different alpha diversity values are compared to local metrics for static variables.
#Define diversity metrics.
DiversityVars <- c("Chao1","Shannon","Simpson")
#Select watershed summary variables.
ModelVars <- colnames(SampledWatershedsWithSummary)[!(colnames(SampledWatershedsWithSummary) %in% c("sum.taxonomy","FieldSeasonNum","HUC12","primer",DiversityVars))]
RemoveVars <- c("mean.CommunityType1000m","stdev.CommunityType1000m","mode.aspect1000m","mode.BurnSeverity1000m","mode.DEF1000m","mode.Elevation1000m","mode.Rainfall1Sample1000m","mode.Rainfall2Sample1000m","mode.Rainfall3Sample1000m","mode.Rainfall4Sample1000m","mode.RFAL20161000m","mode.slope1000m")
ModelVars <- ModelVars[!grepl(paste(RemoveVars,collapse="|"),ModelVars)]
ModelVars <- ModelVars[!grepl(c("mean|stdev|mode"),ModelVars)]
StaticVars <- ModelVars[grepl(c("aspect|CommunityType|DEF|Elevation|RFAL2016|slope"),ModelVars)]
#Aggregate Kruskal-Wallis test results.
StaticLocalTestSummary <- data.frame()
for(Primer in unique(SampledWatershedsWithSummary$primer)){
  watershedsSubset <- SampledWatershedsWithSummary[SampledWatershedsWithSummary$primer==Primer,]
  ModelInput <- watershedsSubset[!duplicated(watershedsSubset),]
  for(Var in StaticVars){
    for(DiversityVar in DiversityVars){
      StaticLocalTest <- data.frame(matrix(nrow=1,ncol=6))
      colnames(StaticLocalTest) <- c("Primer","DiversityVar","Var","ChiSquared","df","p")
      StaticLocalTest[,"Primer"] <- Primer
      StaticLocalTest[,"DiversityVar"] <- DiversityVar
      StaticLocalTest[,"Var"] <- Var
      StaticLocalTest[,"ChiSquared"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$statistic
      StaticLocalTest[,"df"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$parameter
      StaticLocalTest[,"p"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$p.value
      StaticLocalTestSummary <- rbind(StaticLocalTestSummary,StaticLocalTest)
    }
  }
}
#Filter to only contain significant results.
StaticLocalTestSummarySignificant <- StaticLocalTestSummary[StaticLocalTestSummary$p <= 0.05,]

#Test how different alpha diversity values are compared to local or watershed summary metrics for static variables.
#Define diversity metrics.
DiversityVars <- c("Chao1","Shannon","Simpson")
#Select watershed summary variables.
ModelVars <- colnames(SampledWatershedsWithSummary)[!(colnames(SampledWatershedsWithSummary) %in% c("sum.taxonomy","FieldSeasonNum","HUC12","primer",DiversityVars))]
RemoveVars <- c("mean.CommunityType1000m","stdev.CommunityType1000m","mode.aspect1000m","mode.BurnSeverity1000m","mode.DEF1000m","mode.Elevation1000m","mode.Rainfall1Sample1000m","mode.Rainfall2Sample1000m","mode.Rainfall3Sample1000m","mode.Rainfall4Sample1000m","mode.RFAL20161000m","mode.slope1000m")
ModelVars <- ModelVars[!grepl(paste(RemoveVars,collapse="|"),ModelVars)]
DynamicVars <- ModelVars[grepl(c("Rainfall|BurnSeverity"),ModelVars)]
#Aggregate Kruskal-Wallis test results.
DynamicTestSummary <- data.frame()
for(Primer in unique(SampledWatershedsWithSummary$primer)){
  for(FieldSeason in unique(SampledWatershedsWithSummary$FieldSeasonNum)){
    watershedsSubset <- SampledWatershedsWithSummary[SampledWatershedsWithSummary$primer==Primer  & SampledWatershedsWithSummary$FieldSeasonNum==FieldSeason,]
    ModelInput <- watershedsSubset[!duplicated(watershedsSubset),]
    if(nrow(ModelInput)>1){
      for(Var in DynamicVars){
        for(DiversityVar in DiversityVars){
          DynamicTest <- data.frame(matrix(nrow=1,ncol=7))
          colnames(DynamicTest) <- c("Primer","FieldSeasonNum","DiversityVar","Var","ChiSquared","df","p")
          DynamicTest[,"Primer"] <- Primer
          DynamicTest[,"FieldSeasonNum"] <- FieldSeason
          DynamicTest[,"DiversityVar"] <- DiversityVar
          DynamicTest[,"Var"] <- Var
          DynamicTest[,"ChiSquared"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$statistic
          DynamicTest[,"df"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$parameter
          DynamicTest[,"p"] <- kruskal.test(ModelInput[,DiversityVar]~ModelInput[,Var])$p.value
          DynamicTestSummary <- rbind(DynamicTestSummary,DynamicTest)
        }
      }
    }
  }
}
#Filter to only contain comparisons between rainfall metrics gathered in the immediate months leading up to the sample collection.
#Or comparisons between burn severity and diversity metrics.
DynamicTestSummarySignificant <- DynamicTestSummary[DynamicTestSummary$FieldSeasonNum==as.numeric(gsub(".*Rainfall(.+)Sample1000m.*", "\\1", DynamicTestSummary$Var)),]
#DynamicTestSummarySignificant <- DynamicTestSummary
DynamicTestSummarySignificant <- DynamicTestSummarySignificant[complete.cases(DynamicTestSummarySignificant),]
DynamicTestSummarySignificant <- rbind(DynamicTestSummarySignificant,DynamicTestSummary[grepl("BurnSeverity",DynamicTestSummary$Var),])
#Filter to only contain significant results.
DynamicTestSummarySignificant <- DynamicTestSummarySignificant[DynamicTestSummarySignificant$p <= 0.05,]



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
SampleDiversity$primer <- gsub("\\_.*","",SampleDiversity$primer)
SampleDiversity %>%
  ggplot(aes(x=reorder(Field_Season,FieldSeasonNum),y=Simpson,group=factor(Field_Season)))+
  guides(fill=guide_legend(title="Field season"))+
  geom_boxplot(aes(fill=factor(Field_Season)))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")+
  facet_grid(primer~.,labeller=label_both)+
  labs(title="Diversity of samples, classified by primer, using a Simpson index",x="Field Season")

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
  scale_fill_viridis_c()+labs(title="Pearson correlation coefficients between environmental variables (sampled watersheds)")+
  theme(axis.text.y=element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))

