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
# BC_meso_grid <- st_read(dsn="./data", layer="BC_meso_grid")
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
  geom_sf(data=BC_grid_9km2 %>% filter(WID_12km==7), aes(fill=as.factor(Mquad)))


BC_grid_9km2 <- BC_grid_9km2 %>% arrange(WID_12km, Fquad, Mquad)

# create new WID in same order as original WID but no missing values 
uniqueWID <- as.data.frame(unique(BC_grid_9km2$WID_12km))
colnames(uniqueWID)[1] <- "origWID"
uniqueWID$newWID <- rownames(uniqueWID)
tail(uniqueWID)


# create labels and re-order sf object
BC_grid_9km2$Wgrid <- as.numeric(uniqueWID$newWID[match(BC_grid_9km2$WID_12km, uniqueWID$origWID)])
BC_grid_9km2$MesoCell <- paste0("W", str_pad(BC_grid_9km2$WID_12km, 4, pad = "0"),
                                "_F", BC_grid_9km2$Fquad,
                                "_M", BC_grid_9km2$Mquad)

BC_meso_grid <- BC_grid_9km2 %>% select(-Areakm2, -centroidX, -centroidY)
BC_meso_grid <- BC_meso_grid[c("MesoCell","Wgrid","Fquad","Mquad","WID_12km","FID_6km","MID_3km")]

# write shapefile
st_write(BC_meso_grid, paste0(getwd(),"/data/BC_meso_grid.shp"), delete_layer = TRUE)
# MescoCell = character = Wgrid_Fquad_Mquad
# Wgrid = numeric = 1:7601 where lower number = higher priority for random provincial sampling
# Fquad = numeric = 1:4 where number refers to quadrant location
# Mquad = numeric = 1:4 where number refers to quadrant location
# WID_12km = original random ordering for wolverine sized cells (same order but without missing values as Wgrid)
# FID_6km = random ordering for fisher sized cells (if fisher focused study)
# MID_3km = random ordering for marten sized cells (if marten focused study)

