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


# ###--- older code
# # assign GridId to fisher cells
# 
# grid_9km2 <- st_join(grid_36km2, grid_144km2 %>% select(GridId), largest=TRUE)
# grid_36km2 <- grid_36km2 %>% arrange(GridId)
# 
# ggplot()+
#   geom_sf(data = grid_36km2 %>% filter(GridId==14992), aes(fill=GridId))
# 
# tmp <- grid_36km2 %>% count(GridId) %>% st_drop_geometry()
# tmp4cells <- tmp %>% filter(n==4)
# 
# grid_36km2$NumCells <- tmp$n[match(grid_36km2$GridId, tmp$GridId)]
# 
# grid_36km2_4cells <- grid_36km2 %>% filter(NumCells==4)
# nrow(grid_36km2_4cells)/4
# grid_36km2_4cells$FID <- rep(1:4, times=15104)
# 
# grid_36km2_4cells$WID <- grid_36km2_4cells$GridId
# max(grid_36km2_4cells$WID)
# grid_36km2_4cells$GridId <- paste("Id",str_pad(grid_36km2_4cells$WID, 5, pad = "0"),grid_36km2_4cells$FID, sep="_")
# 
# # assign GridId to marten cells
# grid_9km2 <- st_join(grid_9km2, grid_36km2_4cells %>% select(GridId,FID,WID), largest=TRUE)
# grid_9km2 <- grid_9km2 %>% arrange(WID, FID)
# 
# 
# tmp <- grid_9km2 %>% count(GridId) %>% st_drop_geometry()
# tmp %>% count(n)
# tmp4cells <- tmp %>% filter(n==4)
# 
# grid_9km2$NumCells <- tmp$n[match(grid_9km2$GridId, tmp$GridId)]
# 
# grid_9km2_4cells <- grid_9km2 %>% filter(NumCells==4)
# nrow(grid_9km2_4cells)/4
# grid_9km2_4cells <- grid_9km2_4cells %>% arrange(WID, FID)
# grid_9km2_4cells$MID <- rep(1:4, times=59987)
# 
# grid_9km2_4cells$GridId <- paste(grid_9km2_4cells$GridId,grid_9km2_4cells$MID, sep="_")
# 
# ggplot()+
#   geom_sf(data = grid_9km2_4cells %>% filter(WID==1), aes(fill=MID))
# 
# 
# BC_meso_grid <- st_intersection(grid_9km2_4cells %>% st_transform(3005),
#                                 bc_sf %>% st_transform(3005))
# 
# st_write(BC_meso_grid, paste0(getwd(),"/BC_meso_grid.shp"), delete_layer = TRUE)
# 
# ###################################################################################

BC_meso_grid <- st_read(dsn="./data", layer="BC_meso_grid")
FHZ <- st_read(dsn="./data", layer="FHE_zones")

ggplot()+
  geom_sf(data=BC_meso_grid)+
  geom_sf(data=FHZ)

