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

#####################################################################################
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

###--- Clip mesocarnivore grid and NABat grid to study area and produce shapefiles

clip_meso_NABat <- function(dsn=dsn, layer=layer){
  
  input_sf <- st_read(dsn = dsn, layer=layer)
  
  # ggplot()+
  #   geom_sf(data=bc_sf)+
  #   geom_sf(data=input_sf, col="blue", fill="blue")

  SA <- SA_meso(input_sf = input_sf, buff_dist = 3000)
  
  # ggplot()+
  #   geom_sf(data=SA)+
  #   geom_sf(data=input_sf, col="blue", fill="blue")
  
  Fquad <- SA %>% group_by(FID_6km) %>%
    summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
    summarise(across(geometry, ~ st_combine(.)))
  
  Fquad$Wgrid <- SA$Wgrid[match(Fquad$FID_6km, SA$FID_6km)]
  Fquad$Fquad <- SA$Fquad[match(Fquad$FID_6km, SA$FID_6km)]
  
  st_write(Fquad,paste0(getwd(),"/",dsn,"/Fquad.shp"), delete_layer = TRUE)
  
  Wgrid <- SA %>% group_by(Wgrid) %>%
    summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
    summarise(across(geometry, ~ st_combine(.)))
  
  aoi <- SA %>%
    summarise(across(geometry, ~ st_union(.))) %>%
    summarise(across(geometry, ~ st_combine(.)))
  
  # bcdc_search("BEC", res_format = "wms")
  aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
  
  # bcdc_search("indian", res_format = "wms")
  # 1: Indian Reserves - Administrative Boundaries (multiple, wms, kml)
  # ID: 8efe9193-80d2-4fdf-a18c-d531a94196ad
  # Name: indian-reserves-administrative-boundaries
  aoi.IR <- retrieve_geodata_aoi((ID = "8efe9193-80d2-4fdf-a18c-d531a94196ad"))
  
  # ggplot()+
  #   theme_minimal()+
  #   geom_sf(data=aoi)+
  #   geom_sf(data=aoi.BEC %>% filter(grepl("SBS|IDF|SBPS", ZONE)), aes(fill=ZONE))+
  #   geom_sf(data=Wgrid, fill=NA)+
  #   geom_sf(data=input_sf, fill=NA, col="black",lwd=1.2)+
  #   geom_sf(data=aoi.IR, fill="black", col="black",lwd=1.5)
  # 
  
  temp.grid <- st_join(aoi.BEC %>% filter(grepl("SBS|IDF|SBPS", ZONE)), Wgrid)
  Fquad$Fquad <- as.factor(Fquad$Fquad)
  
  tmp.Wgrid <- Fquad %>% filter(Wgrid %in% temp.grid$Wgrid) %>% st_drop_geometry()
  tmp.Wgrid %>% arrange(Wgrid) %>% count(Wgrid)
  rdm.points <- st_centroid(Fquad %>% filter(Wgrid %in% unique(tmp.Wgrid$Wgrid)))
  
  st_write(rdm.points,paste0(getwd(),"/",dsn,"/RndmCamPnts.shp"), delete_layer = TRUE)
  
  Fquad_all <- st_join(Fquad, aoi.BEC, largest=TRUE)
  Fquad_all$use <- ifelse(grepl("SBS|IDF|SBPS", Fquad_all$ZONE),1,0)
  Fquad_all %>% group_by(use) %>% count(ZONE)
  st_write(Fquad_all %>% filter(use==1) %>% select(Wgrid, Fquad, FID_6km, ZONE),
           paste0(getwd(),"/",dsn,"/Priority_Grid_Cells_use.shp"), delete_layer = TRUE)
  
  SA_NABat_all <- BC_NABat %>% st_intersection(aoi)
  SA_NABat_all <- SA_NABat_all %>% filter(!GRTS_ID %in% BC_NABat_2021$GRTS_ID) # exclude the NABat cells already surveyed
  
  # ggplot()+
  #   theme_minimal()+
  #   geom_sf(data=SA_NABat)
  
  st_write(SA_NABat_all %>% select(GRTS_ID, AKCAN5KM_I),paste0(getwd(),"/",dsn,"/SA_NABat.shp"), delete_layer = TRUE)
  SA_NABat <- BC_NABat %>% st_intersection(input_sf)
  
  return(list(SA_NABat=SA_NABat, SA_Fquad=Fquad_all))
  
}

################################################################################
###--- Input files
### NABat
BC_NABat_2021 <- st_read(dsn = "./data/NABat_Total_Sampled_Cells_2021", 
                         layer="NABat_Total_Sampled_Cells") %>% st_transform(crs=3005)

BC_NABat <- st_read(dsn = "data/british_columbia_master_sample",
                    layer="master_sample_british_columbia") %>% st_transform(crs=3005)

### Mesocarnivore grid
BC_meso_grid <- st_read(dsn=paste0(getwd(),"/data"), layer="BC_meso_grid") %>% 
  st_transform(crs=3005)

### BC boundaries
bc_sf <- bc_bound(class = 'sf') %>% st_transform(3005) # to ensure an sf object in Albers

bc_sf <- bc_sf %>% 
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))

bc_buff <- st_buffer(bc_sf, dist=12000)
bc_sf$BC <- "BC" 

### Fisher habitat boundaries
# FHZ <- st_read(dsn="./data", layer="FHE_zones")
# FHZ <- st_simplify(FHZ)


#################################################################################
### Natural Resource Districts
# 10: Natural Resource (NR) Districts (multiple, wms, oracle_sde)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
NRD <- bcdc_get_data(record="0bc73892-e41f-41d0-8d8e-828c16139337")
NRD

###--- for the NR of choice
aoi <- NRD %>% filter(REGION_ORG_UNIT=="RCB") %>% st_transform(3005)

ggplot()+
  geom_sf(data=NRD)+
  geom_sf(data=aoi, col="red")

#################################################################################

# for each Study Area
dsn="data/Esketemc"
layer="SOI_Esket2020"

Esket <- clip_meso_NABat(dsn="data/Esketemc", layer="SOI_Esket2020")

ggplot()+
  # geom_sf(data=bc_sf)+
  geom_sf(data=Esket$SA_NABat, col="red")+
  geom_sf(data=Esket$SA_Fquad %>% filter(use==1), aes(fill=NA, col="blue"))


