###############################################################################
#  Prepare plant layers for optimization with 'Prioritizr' package.
#  April 23, 2018; last updated April 27, 2018
#  Script to project and resample plant SDM layers to match butterfly and vertebrate layers 
#   for use in first round of prioritizations.
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages.
# =============================================================================
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(sf)

	# ----------------------	
	#  Function to normalize rasters to 0 to 1 scale.
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	
	
	
# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------	
	#  Stack all .tif files of endemic butterfly ranges.
	setwd("C:/Users/saraw/Desktop/SDM_Butterflies/Endemics")
	end_fly_ls <- list.files(pattern='\\.tif$')
	end_fly_stack <- stack(end_fly_ls)
	
	# ----------------------	
	#  Stack all .tif files of non-endemic butterfly ranges.
	setwd("C:/Users/saraw/Desktop/SDM_Butterflies/Non-endemics")
	non_end_fly_ls <- list.files(pattern='\\.tif$')
	non_end_fly_stack <- stack(non_end_fly_ls)
	
	# ----------------------	
	#  Stack all .tif files of endemic plant ranges.
	setwd("C:/Users/saraw/Desktop/SDM_Plants/Endemics")
	end_plant_ls <- list.files(pattern='\\.asc$')
	end_plant_stack <- stack(end_plant_ls)
	
	# ----------------------	
	#  Stack all .tif files of non-endemic plant ranges.
	setwd("C:/Users/saraw/Desktop/SDM_Plants/Non-endemics")
	non_end_plant_ls <- list.files(pattern='\\.asc$')
	non_end_plant_stack <- stack(non_end_plant_ls)
	
	# ----------------------	
	# Rare plant spp.
	rare_end_plant_sp <- shapefile("C:/Users/saraw/Desktop/Locked_rare_spp/Endemic/AllEndemicPlants_SquareBuffers.shp")
	rare_non_end_plant_sp <- shapefile("C:/Users/saraw/Desktop/Locked_rare_spp/Non-endemic/Dipterocarp_500mSquareBufferNonEndemic.shp")
	
	
	
# =============================================================================
#  Generate raster template
# =============================================================================		
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	temp <- end_fly_stack[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 
	

	
	
# =============================================================================
#  Project and resample plants SDM rasters to match other inputs.
# =============================================================================
	
	# ----------------------	
	#  Set projection for plant raster stacks.
	projection(end_plant_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
	projection(non_end_plant_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------	
	# Project to CRS we are using.
	newproj <- "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
	end_plant_stack_proj <- projectRaster(end_plant_stack, crs = newproj)
	non_end_plant_stack_proj <- projectRaster(non_end_plant_stack, crs = newproj)
	
	# ----------------------	
	# Resample according to butterfly rasters.
	end_plant_stack_rs <- raster::resample(end_plant_stack_proj, temp, method = 'bilinear')
	non_end_plant_stack_rs <- raster::resample(non_end_plant_stack_proj, temp, method = 'bilinear')
	
	
	
# =============================================================================
#  Create raster input layers for rare plants species (the species that did not have SDMs run).
# =============================================================================	
	
	# ----------------------	
	# Convert to sf object and make a raster value column.
	rare_end_plant_sf <- st_as_sf(rare_end_plant_sp) %>%
		mutate(raster_val = 1)
		
	rare_non_end_plant_sf <- st_as_sf(rare_non_end_plant_sp) %>%
		mutate(Orchid = 0) %>%
		mutate(EnTreesShr = 0) %>%
		mutate(raster_val = 1)	
		
	# ----------------------	
	# Generate empty stack to hole outputs.	
	rare_end_plant_stack <- stack()
	rare_non_end_plant_stack <- stack()
	
	# ----------------------	
	# For each polygon of rare species, convert to sp object and the to raster.
	for(i in 1:nrow(rare_end_plant_sf)){
		sf_tmp <- rare_end_plant_sf[i,]
		sp_tmp <- as(sf_tmp, 'Spatial')
		r_tmp <- rasterize(sp_tmp, r_template, field = sp_tmp$raster_val)
		rare_end_plant_stack <- stack(rare_end_plant_stack, r_tmp)
		}
	
	for(i in 1:nrow(rare_non_end_plant_sf)){
		sf_tmp <- rare_non_end_plant_sf[i,]
		sp_tmp <- as(sf_tmp, 'Spatial')
		r_tmp <- rasterize(sp_tmp, r_template, field = sp_tmp$raster_val)
		rare_non_end_plant_stack <- stack(rare_non_end_plant_stack, r_tmp)
		}
	
	# ----------------------	
	# Add the rare plant stacks to the other plant raster layer stacks.
	all_end_plant_stack <- stack(end_plant_stack_rs, rare_end_plant_stack)
	all_non_end_plant_stack <- stack(non_end_plant_stack_rs, rare_non_end_plant_stack)
	
	
# =============================================================================
#  Resample ACD input to match other feature inputs.
# =============================================================================	
	
	# ----------------------
	#  Load unmasked carbon map from Asner et al.
	acd <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/raw_spat_data/CAO_ACD_30m_unmasked.tif")
	x_min <- acd@data@min
	x_max <- acd@data@max
	acd_tmp <- raster::calc(acd, range01)

	# ----------------------
	#  Resample raster template.
	acd_feat_in <- raster::resample(acd_tmp, temp, method = 'bilinear')
	
	
	
# =============================================================================
#  Resample elevation connectivity input to match other feature inputs.
# =============================================================================	
	
	# ----------------------
	#  Load raster stack from Hodgsen et al.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/stackedresultsreweight.Rdata")
	elev_conn <- condstack2$reweight
	x_min <- elev_conn@data@min
	x_max <- elev_conn@data@max
	elev_conn_tmp <- raster::calc(elev_conn, range01)
	
	# ----------------------
	#  Resample to raster template.
	elev_conn_feat_in <- raster::resample(elev_conn_tmp, temp, method = 'bilinear')
	
	
	
# =============================================================================
#  Save outputs.
# =============================================================================	
	
	# ----------------------	
	# Save.
	#writeRaster(all_end_plant_stack, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/all_end_plant_stack.grd")
	#writeRaster(all_non_end_plant_stack, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/all_non_end_plant_stack.grd")
	#writeRaster(end_fly_stack, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/end_fly_stack.grd")
	#writeRaster(non_end_fly_stack, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/non_end_fly_stack.grd")
	#writeRaster(acd_feat_in, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/acd_feat_in.grd")
	#writeRaster(elev_conn_feat_in, file = "C:/Users/saraw/Desktop/5_5_18/elev_conn_feat_in.grd")
	

	
# =============================================================================	
###############################################################################