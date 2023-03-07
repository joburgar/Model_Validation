# Copyright 2022 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

#####################################################################################
# Mesocarnivore_Sampling_Design_GRTS.R
# script to create provincial GRTS sampling grid, with a mesocarnivore focus
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 18-Nov-2021
#####################################################################################

#####################################################################################
R_version <- paste0("R-",version$major,".",version$minor)
.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# LOAD PACKAGES ####
list.of.packages <- c("tidyverse","sf","PNWColors","nngeo", "bcdata", "bcmaps", "units")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################

bc_sf <- bc_bound(class = 'sf') %>% st_transform(3005) # to ensure an sf object in Albers

bc_sf <- bc_sf %>% 
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))

bc_buff <- st_buffer(bc_sf, dist=12000)
bc_sf$BC <- "BC" 

aoi_utm <- st_transform(bc_sf, crs=26910) # to have in metres for specifying grid cell size
aoi <- st_transform(bc_sf, crs=3005) # to have in metres for specifying grid cell size

# for 12x12 (wolverine)
aoi_grid <- st_make_grid(st_bbox(aoi), cellsize=12000, square=TRUE) #  grid for entire AOI (rectangle)
W_points <- st_point_on_surface(aoi_grid)  # if using portion of aoi
W_points <- st_sf(geom = W_points)
W_points$WID <- sample(nrow(W_points), replace=FALSE)

# for 6x6 (fisher)
aoi_grid <- st_make_grid(st_bbox(aoi), cellsize=6000, square=TRUE) #  grid for entire AOI (rectangle)
F_points <- st_point_on_surface(aoi_grid)  # if using portion of aoi
F_points <- st_sf(geom = F_points)
F_points$FID <- sample(nrow(F_points), replace=FALSE)

# for 3x3 (marten)
aoi_grid <- st_make_grid(st_bbox(aoi), cellsize=3000, square=TRUE) #  grid for entire AOI (rectangle)
rm(aoi_grid)
grid_9km2 = st_sf(geom = aoi_grid)
grid_9km2$Areakm2 <- st_area(grid_9km2)*1e-6
grid_9km2 <- drop_units(grid_9km2)
grid_9km2$Id <- as.numeric(rownames(grid_9km2))

# reduce to the buffered BC border for computation ease
BC_grid_9km2 <- st_intersection(grid_9km2, bc_buff)
ggplot()+
  geom_sf(data=BC_grid_9km2)

# attribute WID to each MID
WID.dist <- st_nn(BC_grid_9km2, W_points, k=1, returnDist = T)
str(WID.dist)

BC_grid_9km2$WID_dist <- unlist(WID.dist$dist)
BC_grid_9km2$WID_tmp <- unlist(WID.dist$nn)
BC_grid_9km2$WID <- W_points$WID[match(BC_grid_9km2$WID_tmp,rownames(W_points))]
BC_grid_9km2$WID_tmp <- NULL
as.data.frame(BC_grid_9km2 %>% filter(WID==8352))

ggplot()+
  geom_sf(data = BC_grid_9km2 %>% filter(WID==100))+
  geom_sf(data = W_points %>% filter(WID==100))

# attribute FID to each MID
FID.dist <- st_nn(BC_grid_9km2, F_points, k=1, returnDist = T)
str(FID.dist)

BC_grid_9km2$FID_dist <- unlist(FID.dist$dist)
BC_grid_9km2$FID_tmp <- unlist(FID.dist$nn)
BC_grid_9km2$FID <- F_points$FID[match(BC_grid_9km2$FID_tmp,rownames(F_points))]
BC_grid_9km2$FID_tmp <- NULL
as.data.frame(BC_grid_9km2 %>% filter(FID==16742))# will need to change this 

ggplot()+
  geom_sf(data = BC_grid_9km2 %>% filter(FID==16219))+
  geom_sf(data = F_points %>% filter(FID==16219))

st_write(BC_grid_9km2, paste0(getwd(),"/BC_grid_9km2.shp"), delete_layer = TRUE)

save(FID.dist, WID.dist, file = "FW_dist.RData")

#################################################################################
BC_grid_9km2 <- st_read(dsn="./data", layer="BC_grid_9km2")
FHZ <- st_read(dsn="./data", layer="FHE_zones")

ggplot()+
  geom_sf(data=BC_grid_9km2)+
  geom_sf(data=FHZ)

ggplot()+
  geom_sf(data=BC_grid_9km2 %>% filter(WID==8352), col="blue")+
  geom_sf(data=BC_grid_9km2 %>% filter(FID==16742))

BC_grid_9km2$MID <- BC_grid_9km2$Id
BC_grid_9km2$Id <- NULL
BC_grid_9km2$WID_dist <- NULL
BC_grid_9km2$FID_dist <- NULL

BC_grid_9km2$Areakm2 <- st_area(BC_grid_9km2)*1e-6
BC_grid_9km2 <- drop_units(BC_grid_9km2)
nrow(BC_grid_9km2 %>% filter(Areakm2!=9)) / 118005 # ~3% < 9km2

BC_grid_9km2 <- BC_grid_9km2 %>% rename(FID_6km=FID, WID_12km=WID, MID_3km=MID)
BC_grid_9km2 <- BC_grid_9km2 %>% arrange(WID_12km, FID_6km, MID_3km)

# create naming convention so that it's WID_Fquadrant_Mquadrant
# follow naming convention of quadrant where 1 is +x+y, 2 is -x+y, 3 is +x-y, and 4 is -x-y
# first find centroids for each cell
centre_points <- BC_grid_9km2 %>% st_transform(3005) %>% st_point_on_surface() %>% st_coordinates()
centre_points <- as.data.frame(centre_points)
BC_grid_9km2$centroidX <- centre_points$X
BC_grid_9km2$centroidY <- centre_points$Y

# work with dataframe for computation ease
BC_grid_df <- BC_grid_9km2 %>% st_drop_geometry() %>% as_tibble()

# create Fquad first (subset to WID and FID)
BC_WF_grid_df <- BC_grid_df %>% group_by(WID_12km) %>% count(FID_6km)
BC_WF_grid_df$centroidX <- BC_grid_df$centroidX[match(BC_WF_grid_df$FID_6km, BC_grid_df$FID_6km)]
BC_WF_grid_df$centroidY <- BC_grid_df$centroidY[match(BC_WF_grid_df$FID_6km, BC_grid_df$FID_6km)]

WID_min_centroid <- BC_WF_grid_df %>% group_by(WID_12km) %>% summarise(WIDminX = min(centroidX), WIDminY = min(centroidY))
BC_WF_grid_df <- left_join(BC_WF_grid_df, WID_min_centroid)

BC_WF_grid_df <- BC_WF_grid_df %>% mutate(Fquad = case_when(centroidX>WIDminX & centroidY>WIDminY ~ 1,
                                                            centroidX==WIDminX & centroidY>WIDminY ~ 2,
                                                            centroidX==WIDminX & centroidY==WIDminY ~ 3,
                                                            centroidX>WIDminX & centroidY==WIDminY ~ 4))
# add the Fquad back to the spatial file and plot to check
BC_grid_9km2 <- left_join(BC_grid_9km2, BC_WF_grid_df %>% ungroup() %>% select(FID_6km, Fquad))

ggplot()+
  geom_sf(data=BC_grid_9km2 %>% filter(WID_12km==1), aes(fill=Fquad))


# create Mquad next (subset to FID and MID)
BC_FM_grid_df <- BC_grid_df %>% select(-Areakm2, -WID_12km)

FID_min_centroid <- BC_FM_grid_df %>% group_by(FID_6km) %>% summarise(FIDminX = min(centroidX), FIDminY = min(centroidY))
BC_FM_grid_df <- left_join(BC_FM_grid_df, FID_min_centroid)

BC_FM_grid_df <- BC_FM_grid_df %>% mutate(Mquad = case_when(centroidX>FIDminX & centroidY>FIDminY ~ 1,
                                                            centroidX==FIDminX & centroidY>FIDminY ~ 2,
                                                            centroidX==FIDminX & centroidY==FIDminY ~ 3,
                                                            centroidX>FIDminX & centroidY==FIDminY ~ 4))


# add the Mquad back to the spatial file and plot to check
BC_grid_9km2 <- left_join(BC_grid_9km2, BC_FM_grid_df %>% ungroup() %>% select(MID_3km, Mquad))

ggplot()+
  geom_sf(data=BC_grid_9km2 %>% filter(WID_12km==7), aes(fill=Mquad))


### subset to just BC
BC_grid_9km2 <- BC_grid_9km2 %>% arrange(WID_12km, Fquad, Mquad)
BC_meso_grid <- BC_grid_9km2 %>% select(-Areakm2, -centroidX, -centroidY)

BC_meso_grid <- st_intersection(BC_meso_grid, bc_sf)

BC_meso_grid <- BC_meso_grid[c("MesoCell","Wgrid","Fquad","Mquad","WID_12km","FID_6km","MID_3km")]


# create new WID in same order as original WID but no missing values 
uniqueWID <- as.data.frame(unique(BC_meso_grid$WID_12km))
colnames(uniqueWID)[1] <- "origWID"
uniqueWID$newWID <- rownames(uniqueWID)
tail(uniqueWID)


# create labels and re-order sf object
BC_meso_grid$Wgrid <- as.numeric(uniqueWID$newWID[match(BC_meso_grid$WID_12km, uniqueWID$origWID)])
BC_meso_grid$MesoCell <- paste0("W", str_pad(BC_meso_grid$Wgrid, 4, pad = "0"),
                                "_F", BC_meso_grid$Fquad,
                                "_M", BC_meso_grid$Mquad)

tail(BC_meso_grid)
# write shapefile
st_write(BC_meso_grid, paste0(getwd(),"/data/BC_meso_grid.shp"), delete_layer = TRUE)

st_bbox(BC_meso_grid)
# xmin      ymin      xmax      ymax 
# 275942.4  367537.4 1867409.8 1735251.6 
# MescoCell = character = Wgrid_Fquad_Mquad
# Wgrid = numeric = 1:7601 where lower number = higher priority for random provincial sampling
# Fquad = numeric = 1:4 where number refers to quadrant location
# Mquad = numeric = 1:4 where number refers to quadrant location
# WID_12km = original random ordering for wolverine sized cells (same order but without missing values as Wgrid)
# FID_6km = random ordering for fisher sized cells (if fisher focused study)
# MID_3km = random ordering for marten sized cells (if marten focused study)


#################################################################################
# View shapefile
BC_meso_grid <- st_read(dsn=paste0(getwd(),"/data"), layer="BC_meso_grid")

glimpse(BC_meso_grid)

ggplot()+
    geom_sf(data=BC_meso_grid %>% filter(Wgrid==2), aes(fill=Mquad))

ggplot()+
  geom_sf(data=BC_meso_grid %>% filter(Wgrid==2), aes(fill=Fquad))+
  geom_sf(data=BC_meso_grid %>% filter(FID_6km==29599), aes(fill=Mquad))


# BC_meso_Wgrid <- BC_meso_grid %>% group_by(Wgrid) %>%
#   summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
#   summarise(across(geometry, ~ st_combine(.)))
# 
# st_write(BC_meso_Wgrid, paste0(getwd(),"/data/BC_meso_Wgrid.shp"), delete_layer = TRUE)

BC_meso_Wgrid <- st_read(dsn=paste0(getwd(),"/data"), layer="BC_meso_Wgrid")

# BC_meso_Fquad <- BC_meso_grid %>% group_by(FID_6km) %>%
#   summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
#   summarise(across(geometry, ~ st_combine(.)))
# 
# st_write(BC_meso_Fquad, paste0(getwd(),"/data/BC_meso_Fquad.shp"), delete_layer = TRUE)

BC_meso_Fquad <- st_read(dsn=paste0(getwd(),"/data"), layer="BC_meso_Fquad")

ggplot()+
  geom_sf(data=bc_sf) +
  geom_sf(data=BC_meso_Wgrid)

#################################################################################
###--- FUNCTIONS ---###

###--- Bring in study area - clip grid to buffered area and export
SA_meso <- function(input_sf = input_sf, buff_dist = 12000){
  sa_input_buff <- st_buffer(input_sf, dist = buff_dist)
  find_grid_cells <- sf::st_within(BC_meso_grid, st_bbox(sa_input_buff) %>% st_as_sfc(.))
  filt_grid_cells <- BC_meso_grid[which(lengths(find_grid_cells) != 0), ]
  SA_meso <- BC_meso_grid %>% filter(Wgrid %in% filt_grid_cells$Wgrid)
  return(SA_meso)
}

###--- Retrieve data from the BC data catalogue
retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  return(aoi.geodata)
}

#################################################################################
# # for Leah Anderson (DECAR_Wildlife_LSA)
# input_sf <- st_read(dsn=paste0(getwd(),"/data/Decar_Wildlife_LSA"), layer="Decar_Wildlife_LSA") %>% st_transform(crs=3005)
# 
# DECAR_Wildlife_LSA <- SA_meso(input_sf = input_sf)
# 
# ggplot()+
#   geom_sf(data=DECAR_Wildlife_LSA, aes(fill=Mquad))+
#   geom_sf(data=input_sf, fill=NA)
# 
# st_write(DECAR_Wildlife_LSA, paste0(getwd(),"./data/Decar_Wildlife_LSA/meso_grid_DECAR_Wildlife_LSA.shp"), delete_layer = TRUE)

#################################################################################
# # for Enterprise Telemetry (DRAFT_Enterprise_telemetry)
GIS_Dir <- "//Sfp.idir.bcgov/s140/S40203/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Enterprise fisher telemetry/GIS_data"
input_sf <- st_read(dsn=GIS_Dir, layer="DRAFT_Enterprise_telemetry") %>% st_transform(crs=3005)

Enterprise_SA <- SA_meso(input_sf = input_sf, buff_dist = 48000)

ggplot()+
  geom_sf(data=Enterprise_SA, aes(fill=Wgrid))+
  geom_sf(data=input_sf, fill=NA, lwd = 2)

st_write(Enterprise_SA, paste0(getwd(),"./data/meso_grid_Enterprise_SA.shp"), delete_layer = TRUE)

# Enterprise_SA %>% count(Wgrid) %>% st_drop_geometry()
# 
# Enterprise_Fquad <- Enterprise_SA %>% group_by(FID_6km) %>%
#   summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
#   summarise(across(geometry, ~ st_combine(.)))
# 
# aoi <- Enterprise_SA %>%
#   summarise(across(geometry, ~ st_union(.))) %>%
#   summarise(across(geometry, ~ st_combine(.)))
# 
# bcdc_search("BEC", res_format = "wms")
# aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
# 
# # bcdc_search("indian", res_format = "wms")
# # 1: Indian Reserves - Administrative Boundaries (multiple, wms, kml)
# # ID: 8efe9193-80d2-4fdf-a18c-d531a94196ad
# # Name: indian-reserves-administrative-boundaries
# aoi.IR <- retrieve_geodata_aoi((ID = "8efe9193-80d2-4fdf-a18c-d531a94196ad"))
# 
# # bcdc_search("city", res_format = "wms")
# # 2: BC Major Cities Points 1:2,000,000 (Digital Baseline Mapping) (multiple, wms, kml)
# # ID: b678c432-c5c1-4341-88db-0d6befa0c7f8
# aoi.city <- retrieve_geodata_aoi(ID = "b678c432-c5c1-4341-88db-0d6befa0c7f8")
# 
# bcdc_search("road", res_format = "wms")
# # 1: Digital Road Atlas (DRA) - Master Partially-Attributed Roads (multiple, zip, wms, pdf, kml)
# # ID: bb060417-b6e6-4548-b837-f9060d94743e
# # Name: digital-road-atlas-dra-master-partially-attributed-roads
# aoi.road <- retrieve_geodata_aoi(ID = "bb060417-b6e6-4548-b837-f9060d94743e")
# 
# 
# ggplot()+
#   geom_sf(data=aoi)+
#   geom_sf(data=aoi.BEC, aes(fill=ZONE))+
#   geom_sf(data=Enterprise_Fquad, fill=NA, lwd=0.3, col="lightgrey")+
#   geom_sf(data=input_sf, fill=NA, lwd=1.5)+
#   geom_sf(data=aoi.city, aes(col=NAME), cex=5)+
#   geom_sf(data=aoi.IR %>% filter(grepl("CANOE|ALKALI", ENGLISH_NAME)), aes(fill=NA))
# 
# 
# as.data.frame(aoi.IR %>% count(ENGLISH_NAME))
#################################################################################
# for Gitanyow Study Area (DRAFT_Enterprise_telemetry)
# input_sf <- st_read("./data/Gitanyow/Gitanyow_draft_study_area.kml") %>% st_transform(crs=3005)
# 
# Gitanyow_Boundary <- st_read(dsn="./data/Gitanyow", layer="Gitanyow_House_Territory_Outside_Boundary") %>% st_transform(crs=3005)
# 
# ggplot()+
#   geom_sf(data=Gitanyow_Boundary)+
#   geom_sf(data=input_sf)
#   
# Gitanyow_SA <- SA_meso(input_sf = input_sf, buff_dist = 10000)
# 
# ggplot()+
#   geom_sf(data=Gitanyow_SA, aes(fill=Wgrid))+
#   geom_sf(data=input_sf, fill=NA, lwd = 2)
# 
# st_write(Gitanyow_SA, paste0(getwd(),"./data/meso_grid_Gitanyow_SA.shp"), delete_layer = TRUE)
# 
# Gitanyow_SA %>% count(Wgrid) %>% st_drop_geometry()
# 
# aoi <- Gitanyow_SA %>%
#   summarise(across(geometry, ~ st_union(.))) %>%
#   summarise(across(geometry, ~ st_combine(.)))
# 
# # bcdc_search("BEC", res_format = "wms")
# aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
# 
# # bcdc_search("city", res_format = "wms")
# # # 2: BC Major Cities Points 1:2,000,000 (Digital Baseline Mapping) (multiple, wms, kml)
# # # ID: b678c432-c5c1-4341-88db-0d6befa0c7f8
# # aoi.city <- retrieve_geodata_aoi(ID = "b678c432-c5c1-4341-88db-0d6befa0c7f8")
# 
# ggplot()+
#   geom_sf(data=aoi)+
#   geom_sf(data=aoi.BEC, aes(fill=MAP_LABEL))+
#   geom_sf(data=Gitanyow_Fquad, fill=NA, lwd=0.3, col="lightgrey")+
#   geom_sf(data=input_sf, fill=NA, lwd=1.5)
#   # geom_sf(data=Gitanyow_Boundary, fill=NA, lwd=0.3, col="black")
#   
# ggplot()+
#   geom_sf(data=Gitanyow_Fquad, aes(fill=FID_6km))+
#   # geom_sf(data=Gitanyow_SA, aes(fill=FID_6km))+
#   geom_sf(data=input_sf, fill=NA, lwd=1.5)
# # geom_sf(data=Gitanyow_Boundary, fill=NA, lwd=0.3, col="black")
# 
# 
# Gitanyow_Fquad <- Gitanyow_SA %>% group_by(FID_6km) %>%
#   summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
#   summarise(across(geometry, ~ st_combine(.)))


#################################################################################
# for Esk'etemc Study Area ()
input_sf <- st_read(dsn = "./data/Esketemc", layer="SOI_Esket2020")

ggplot()+
  geom_sf(data=bc_sf)+
  geom_sf(data=input_sf, col="blue", fill="blue")

Esketemc_SA <- SA_meso(input_sf = input_sf, buff_dist = 3000)

ggplot()+
  geom_sf(data=Esketemc_SA)+
  geom_sf(data=input_sf, col="blue", fill="blue")

Esketemc_Fquad <- Esketemc_SA %>% group_by(FID_6km) %>%
  summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
  summarise(across(geometry, ~ st_combine(.)))

Esketemc_Fquad$Wgrid <- Esketemc_SA$Wgrid[match(Esketemc_Fquad$FID_6km, Esketemc_SA$FID_6km)]
Esketemc_Fquad$Fquad <- Esketemc_SA$Fquad[match(Esketemc_Fquad$FID_6km, Esketemc_SA$FID_6km)]

st_write(Esketemc_Fquad,paste0(getwd(),"/data/Esketemc/Esketemc_Fquad.shp"), delete_layer = TRUE)


Esketemc_Wgrid <- Esketemc_SA %>% group_by(Wgrid) %>%
  summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
  summarise(across(geometry, ~ st_combine(.)))

aoi <- Esketemc_SA %>%
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))

# bcdc_search("BEC", res_format = "wms")
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")

# bcdc_search("indian", res_format = "wms")
# 1: Indian Reserves - Administrative Boundaries (multiple, wms, kml)
# ID: 8efe9193-80d2-4fdf-a18c-d531a94196ad
# Name: indian-reserves-administrative-boundaries
aoi.IR <- retrieve_geodata_aoi((ID = "8efe9193-80d2-4fdf-a18c-d531a94196ad"))
aoi.IR$FN <- word(aoi.IR$ENGLISH_NAME)
aoi.EsCC <- aoi.IR %>% filter(grepl("DOG|ALKALI", FN))

# bcdc_search("city", res_format = "wms")
# 2: BC Major Cities Points 1:2,000,000 (Digital Baseline Mapping) (multiple, wms, kml)
# ID: b678c432-c5c1-4341-88db-0d6befa0c7f8
aoi.city <- retrieve_geodata_aoi(ID = "b678c432-c5c1-4341-88db-0d6befa0c7f8")

# bcdc_search("road", res_format = "wms")
# 1: Digital Road Atlas (DRA) - Master Partially-Attributed Roads (multiple, zip, wms, pdf, kml)
# ID: bb060417-b6e6-4548-b837-f9060d94743e
# Name: digital-road-atlas-dra-master-partially-attributed-roads
aoi.road <- retrieve_geodata_aoi(ID = "bb060417-b6e6-4548-b837-f9060d94743e")
aoi.road %>% count(ROAD_SURFACE) %>% st_drop_geometry()
as.data.frame(aoi.road %>% count(ROAD_CLASS) %>% st_drop_geometry())

ggplot()+
  theme_minimal()+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.BEC %>% filter(grepl("SBS|IDF|SBPS", ZONE)), aes(fill=ZONE))+
  geom_sf(data=Esketemc_Wgrid, fill=NA)+
  # geom_sf(data=Esketemc_Fquad, fill=NA, lwd=0.3, col="lightgrey")+
  geom_sf(data=input_sf, fill=NA, col="black",lwd=1.2)+
  # geom_sf(data=aoi.road %>% filter(FEATURE_TYPE=="Road"))+
  # geom_sf(data=aoi.road %>% filter(ROAD_SURFACE==c("paved","loose")))+
  # geom_sf(data=aoi.road %>% filter(ROAD_CLASS==c("highway","local","service")))+
  geom_sf(data=aoi.EsCC, fill="black", col="black",lwd=1.5)


temp.grid <- st_join(aoi.BEC %>% filter(grepl("SBS|IDF|SBPS", ZONE)), Esketemc_Wgrid)
Esketemc_Fquad$Fquad <- as.factor(Esketemc_Fquad$Fquad)

tmp.Wgrid <- Esketemc_Fquad %>% filter(Wgrid %in% temp.grid$Wgrid) %>% st_drop_geometry()
tmp.Wgrid %>% arrange(Wgrid) %>% count(Wgrid) #1599, 1673
rdm.points <- st_centroid(Esketemc_Fquad %>% filter(Wgrid %in% unique(tmp.Wgrid$Wgrid)))

st_write(rdm.points,paste0(getwd(),"/data/Esketemc/RndmCamPnts.shp"), delete_layer = TRUE)

rdm.points %>% st_transform(crs=4326)

ggplot()+
  theme_minimal()+
  geom_sf(data=aoi.BEC %>% filter(grepl("SBS|IDF|SBPS", ZONE)), aes(fill=ZONE))+
  geom_sf(data=input_sf, fill=NA, col="black",lwd=1.2)+
  geom_sf(data=Esketemc_Fquad %>% filter(Wgrid %in% c("1599","1673")), fill=NA)+
  geom_sf(data=rdm.points, col="blue")+
  # geom_sf(data=Esketemc_Fquad %>% filter(Wgrid %in% c("1599","1673")), aes(fill=Fquad))+
  # geom_sf(data=Esketemc_Wgrid %>% filter(Wgrid %in% temp.grid$Wgrid), col="black", fill=NA, lwd=1.2)+
  # geom_sf(data=aoi.road %>% filter(FEATURE_TYPE=="Road"))+
  # geom_sf(data=aoi.road %>% filter(ROAD_SURFACE==c("paved","loose")))+
  # geom_sf(data=aoi.road %>% filter(ROAD_CLASS==c("highway","local","service")))+
  geom_sf(data=aoi.EsCC, fill="black", col="black",lwd=1.5)

pal <- pnw_palette(name="Shuksan2",n=4, type="discrete")


ggplot()+
  geom_sf(data=Esketemc_Fquad, aes(fill=Fquad))+
  scale_color_manual(values = pnw_palette("Shuksan2",4))+
  # geom_sf(data=Esketemc_Wgrid, col="darkblue", lwd=1, fill=NA)+
  # geom_sf(data=input_sf, col="black", lwd=2, fill=NA)+
  geom_sf(data=aoi.EsCC, lwd=0)


Esketemc_Fquad_inner <- Esketemc_Fquad %>% st_intersection(input_sf)
Esketemc_Fquad_inner <- st_join(Esketemc_Fquad_inner, aoi.BEC, largest=TRUE)
Esketemc_Fquad_inner$use <- ifelse(grepl("SBS|IDF|SBPS", Esketemc_Fquad_inner$ZONE),1,0)
Esketemc_Fquad_inner %>% group_by(use) %>% count(ZONE)
Esketemc_Fquad_inner <- Esketemc_Fquad_inner %>% arrange(!use, FID_6km)
Esketemc_Fquad_inner$priority <- seq(1:nrow(Esketemc_Fquad_inner))

st_write(Esketemc_Fquad_inner %>% select(Wgrid, Fquad, FID_6km, ZONE, priority),paste0(getwd(),"/data/Esketemc/Priority_Grid_Cells.shp"), delete_layer = TRUE)



ggplot()+
  geom_sf(data=Esketemc_Fquad_inner)
rlang::last_error()
  
#################################################################################
FHZ <- st_read(dsn="./data", layer="FHE_zones")
st_simplify(FHZ)

getwd()
# fgdb = "/TTML_PlantInventory/p12/bcsarga.gdb"
fgdb = "/M-PEPE-SPI/M-PEPE_SpatialObservations_Aug2022.gdb"

# List  feature classes in a file geodatabase
st_layers(paste0(getwd(),fgdb))

# for the 06-Aug-2021 output, the p20 gdb was imported
# Available layers:
#                          layer_name geometry_type features fields          crs_name
# 1     SPI_SURVEY_OBS_M_PEPE_AUG2022         Point      153    133 NAD83 / BC Albers
# 2 SPI_INCIDENTAL_OBS_M_PEPE_AUG2022         Point       51     65 NAD83 / BC Albers
# 3  SPI_TELEMETRY_OBS_M_PEPE_AUG2022         Point     3169     39 NAD83 / BC Albers

# just importing one layer
SPI_SURVEY <- st_read(dsn=paste0(getwd(),fgdb),layer="SPI_SURVEY_OBS_M_PEPE_AUG2022") %>%
    st_transform(crs=3005)%>% select("SURVEY_OBS","SPECIES_EN","OBSERVED_N","OBSERVATIO","OBSERVAT_1","LATITUDE","LONGITUDE")
as.data.frame(SPI_SURVEY %>% count(OBSERVAT_1) %>% st_drop_geometry())

SPI_INC <- st_read(dsn=paste0(getwd(),fgdb),layer="SPI_INCIDENTAL_OBS_M_PEPE_AUG2022") %>%
  st_transform(crs=3005)
as.data.frame(SPI_INC %>% count(OBSERVAT_1) %>% st_drop_geometry())

SPI_TELEM <- st_read(dsn=paste0(getwd(),fgdb),layer="SPI_TELEMETRY_OBS_M_PEPE_AUG2022") %>%
  st_transform(crs=3005)
as.data.frame(SPI_TELEM %>% count(OBSERVAT_1) %>% st_drop_geometry())

ggplot()+
  geom_sf(data=bc_sf)+
  # geom_sf(data=FHZ)+
  geom_sf(data=SPI_SURVEY) +
  geom_sf(data=SPI_INC, col="red") +
  geom_sf(data=SPI_TELEM, col="blue")+
  geom_sf(data=input_sf, fill=NA)






