###############################################################################
#  Prepare plant layers for optimization with 'Prioritizr' package.
#  April 23, 2018; last updated October 1, 2018
#  Script to prepare all conservation inputs to matching resolution, extent and projection
#   so they can be stacked and input into prioritization. These all have the areas of 
#   Sabah's TPAs removed.
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

	
	
# =============================================================================
#  Load raster data.
# =============================================================================		

	# ----------------------	
	#  Stack all .tif files of endemic butterfly ranges.
	setwd("C:/Users/saraw/Documents/Prioritization/feature_prep/SDM_Butterflies/Endemic_15")
	end_fly_ls <- list.files(pattern='\\.tif$')
	end_fly_stack <- stack(end_fly_ls)
	
	# ----------------------	
	#  Stack all .tif files of non-endemic butterfly ranges.
	setwd("C:/Users/saraw/Documents/Prioritization/feature_prep/SDM_Butterflies/NonEndemic_62")
	non_end_fly_ls <- list.files(pattern='\\.tif$')
	non_end_fly_stack <- stack(non_end_fly_ls)
	
	# ----------------------	
	#  Stack all .tif files of endemic plant ranges.
	setwd("C:/Users/saraw/Documents/Prioritization/feature_prep/SDM_Plants/Endemics")
	end_plant_ls <- list.files(pattern='\\.asc$')
	end_plant_stack <- stack(end_plant_ls)
	
	# ----------------------	
	#  Stack all .tif files of non-endemic plant ranges.
	setwd("C:/Users/saraw/Documents/Prioritization/feature_prep/SDM_Plants/Non-endemics")
	non_end_plant_ls <- list.files(pattern='\\.asc$')
	non_end_plant_stack <- stack(non_end_plant_ls)
	
	

# =============================================================================
#  Load shapefile data.
# =============================================================================		
	
	setwd("C:/Users/saraw/Documents/Prioritization/")
	
	# ----------------------
	#  Mammal spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_mammals_CR_EN_VU"
	load(file = "vert_range_data/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_mammals_CR_EN_VU.Rdata")
	
	# ----------------------
	#  Amphibian spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_amphibians_CR_EN_VU"
	load(file = "vert_range_data/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_amphibians_CR_EN_VU.Rdata")
	
	# ----------------------
	#  Bird spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_birds_CR_EN_VU"
	load(file = "vert_range_data/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_birds_CR_EN_VU.Rdata")

	# ----------------------	
	# Rare plant spp.
	rare_end_plant_sp <- shapefile("feature_prep/Locked_rare_spp/Endemic/AllEndemicPlants_SquareBuffers.shp")
	rare_non_end_plant_sp <- shapefile("feature_prep/Locked_rare_spp/Non-endemic/Dipterocarp_500mSquareBufferNonEndemic.shp")
	
	# ----------------------
	# Existing TPA 
	load("study_area_boundaries_and_pas/sabah_tpa.Rdata")
	
	# ----------------------
	#  Load Sabah mainland border
	load("study_area_boundaries_and_pas/main_sabah_sf.Rdata")
	
	
	
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
	
	# ----------------------	
	#  Function to normalize rasters to 0 to 1 scale.
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	
	# ----------------------
	#  Create sp obect of Sabah mainland border
	border_sp <- as(main_sabah_sf, "Spatial")
	
	
	
# =============================================================================
#  PLANTS
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
	
	# ----------------------	
	# Convert to rare plants (no SDMs) to sf object and make a raster value column.
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
	# For each polygon of rare species, convert to sp object and then to raster.
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
	
	# ----------------------
	#  Combine all plant layers (endemic and non-endemic) that occur outside PAs
	plant_all <- stack(all_end_plant_stack, all_non_end_plant_stack)
	plant_all[is.na(plant_all)] <- 0
	
	# ----------------------
	#  Generate weighting vector.
	plant_rep_weight_w_pa_range <- c(rep(2, nlayers(all_end_plant_stack)), rep(1, nlayers(all_non_end_plant_stack)))
	
	# ----------------------
	#  Save outputs.
	#writeRaster(plant_all, file = "Prioritization/feature_inputs/plant_all_w_pa_range.grd")
	#save(plant_rep_weight_w_pa_range, file = "Prioritization/feature_inputs/plant_rep_weight_w_pa_range.Rdata")

	

# =============================================================================
#  BUTTERFLIES
# =============================================================================
  	
	# ----------------------
	#  Combine all butterfly layers (endemic and non-endemic)
 	fly_all <- stack(end_fly_stack, non_end_fly_stack)
	fly_all[is.na(fly_all)] <- 0
	
 	# ----------------------
	#  Generate weighting vector.
	fly_rep_weight_w_pa_range <- c(rep(2, nlayers(end_fly_stack)), rep(1, nlayers(non_end_fly_stack)))
	
	# ----------------------
	#  Save outputs.
	#writeRaster(fly_all, file = "feature_inputs/fly_all_w_pa_range.grd")
	#save(fly_rep_weight_w_pa_range, file = "feature_inputs/fly_rep_weight_w_pa_range.Rdata")
	
	

# =============================================================================
#  VERTEBRATES
# =============================================================================	
	
	# ----------------------
	#  Prep sf objects for mammals
	sf_m_tmp <- sf_adj_df_mammals_CR_EN_VU %>%
		mutate(ras_val = 1) %>%
		mutate(threat_wt = ifelse(RL_code == "VU", 2, 
			ifelse(RL_code == "EN", 3, 4))) %>%
		dplyr::filter(genus !=  "Dicerorhinus") %>% # remove rhino
		dplyr::select(id_no, binomial,  RL_code,  ras_val , threat_wt, order, family, 
			forest_dep, elev_min, elev_max) 
	df_m_tmp <- sf_m_tmp
	st_geometry(df_m_tmp) <- NULL
	df_m <- droplevels(df_m_tmp)
	geom_m <- sf_m_tmp$geometry
	sf_m <- st_as_sf(cbind(df_m, geom_m))
	
	# ----------------------
	#  Rasterize and stack each mammal range
	feat_m <- stack()
	for(i in 1:nrow(sf_m)){
		tmp_sf <- sf_m %>%
			dplyr::filter(row_number() == i)
		tmp_sp <- as(tmp_sf, "Spatial")
		tmp_r <- rasterize(tmp_sp, r_template, value = tmp_sp$ras_val)
		feat_m <- stack(feat_m, tmp_r)
		}
	names(feat_m) <- sf_m$binomial
	
	# ----------------------
	#  Prep sf objects for amphibians
	sf_a_tmp <-  sf_adj_df_amphibians_CR_EN_VU %>%
		mutate(ras_val = 1) %>%
		mutate(threat_wt = ifelse(RL_code == "VU", 2, 
			ifelse(RL_code == "EN", 3, 4))) %>%
		dplyr::select(id_no, binomial, RL_code,  ras_val , threat_wt, order, family, 
			forest_dep, elev_min, elev_max)
	df_a_tmp <- sf_a_tmp
	st_geometry(df_a_tmp) <- NULL
	df_a <- droplevels(df_a_tmp)
	df_a <- df_a[1:23,] #  remove species that only occurs in Sarawak
	geom_a <- sf_a_tmp$geometry
	sf_a <- st_as_sf(cbind(df_a, geom_a))
	
	# ----------------------
	#  Rasterize and stack each amphibian range
	feat_a <- stack()
	for(i in 1:nrow(sf_a)){
		tmp_sf <- sf_a %>%
			dplyr::filter(row_number() == i)
		tmp_sp <- as(tmp_sf, "Spatial")
		tmp_r <- rasterize(tmp_sp, r_template, value = tmp_sp$ras_val)
		feat_a <- stack(feat_a, tmp_r)
		}
	names(feat_a) <- sf_a$binomial
	
	# ----------------------
	#  Prep sf objects for birds
	sf_b_tmp <-  sf_adj_df_birds_CR_EN_VU %>%
		mutate(ras_val = 1) %>%
		mutate(threat_wt = ifelse(RL_code == "VU", 2, 
			ifelse(RL_code == "EN", 3, 4))) %>%
		dplyr::filter(binomial !=  "Fregata andrewsi") %>%
		dplyr::select(id_no, binomial, RL_code,  ras_val , threat_wt, order, family, 
			forest_dep, elev_min, elev_max)
	df_b_tmp <- sf_b_tmp
	st_geometry(df_b_tmp) <- NULL
	df_b <- droplevels(df_b_tmp)
	geom_b <- sf_b_tmp$geometry
	sf_b <- st_as_sf(cbind(df_b, geom_b))
	
	# ----------------------
	#  Rasterize and stack each bird range
	feat_b <- stack()
	for(i in 1:nrow(sf_b)){
		tmp_sf <- sf_b %>%
			dplyr::filter(row_number() == i)
		tmp_sp <- as(tmp_sf, "Spatial")
		tmp_r <- rasterize(tmp_sp, r_template, value = tmp_sp$ras_val)
		feat_b <- stack(feat_b, tmp_r)
		}
	names(feat_b) <- sf_b$binomial
	
	# ----------------------
	#  Combine all taxa in a single stack.
	sf_all <- rbind(sf_m, sf_a, sf_b)
	vert_all <- stack(feat_m, feat_a, feat_b)
	vert_all[is.na(vert_all)] <- 0
	
	vert_rep_weight_w_pa_range <- c(df_m$threat_wt, df_a$threat_wt, df_b$threat_wt)
	
	# ----------------------
	#  Save outputs.
	#writeRaster(vert_all, file = "Prioritization/feature_inputs/vert_all_w_pa_range.grd")
	#save(vert_rep_weight_w_pa_range, file = "Prioritization/feature_inputs/vert_rep_weight_w_pa_range.Rdata")
	

	
# =============================================================================
#  Resample ACD input to match other feature inputs.
# =============================================================================	
	
	setwd("C:/Users/saraw/Documents/")
	
	# ----------------------
	#  Process unmasked carbon map from Asner et al.
	acd <- raster("Spatial_data/raw_spat_data/CAO_ACD_30m_unmasked.tif")
	acd_tmp1 <- raster::resample(acd, r_template, method = 'bilinear')
	
	x_min <- cellStats(acd_tmp1, min)
	x_max <- cellStats(acd_tmp1, max)
	acd_feat_in <- raster::calc(acd_tmp1, range01)

	#writeRaster(acd_feat_in, "Prioritization/feature_inputs/acd_feat_in_w_pa_range.grd")
	
	
	
# =============================================================================
#  Resample elevation connectivity input to match other feature inputs.
# =============================================================================	
	
	# ----------------------
	#  Load raster stack from Hodgsen et al.
	load("Prioritization/feature_prep/stackedresultsreweight.Rdata")
	elev_conn <- condstack2$reweight
	elev_conn_tmp1 <- raster::resample(elev_conn, r_template, method = 'bilinear')
	elev_conn_tmp2 <- raster::mask(elev_conn_tmp1, border_sp, updateValue = 0)
		
	x_min <- elev_conn_tmp2@data@min
	x_max <- elev_conn_tmp2@data@max
	elev_conn_feat_in <- raster::calc(elev_conn_tmp2, range01)

	#writeRaster(elev_conn_feat_in,"Prioritization/feature_inputs/elev_conn_feat_in_w_pa_range.grd")

	
	
# =============================================================================
#  Resample corridor connectivity input to match other feature inputs.
# =============================================================================		
	
	# ----------------------
	#  Resample to raster template.
	corr_r  <- raster("Prioritization/feature_prep/corr_sm.grd")
	corr_conn_tmp1 <- raster::resample(corr_r, r_template, method = 'bilinear')
	
	x_min <- corr_conn_tmp1@data@min
	x_max <- corr_conn_tmp1@data@max
	corr_conn_feat_in <- raster::calc(corr_conn_tmp1, range01)

	#writeRaster(corr_conn_feat_in,"Prioritization/feature_inputs/corr_conn_feat_in_w_pa_range.grd")
	
	
	
# =============================================================================
#  Forest Formations
# =============================================================================	
	
	# ----------------------	
	#  Read in forest formations raw layer (polygon) and removes mangroves (as we are not 
	#   interested in protecting mangroves via this project.
	for_form_tmp <- shapefile("Spatial_data/ForestFormations/ForestFormations_Sabah_extended.shp")
	for_form_sf <- st_as_sf(for_form_tmp) %>%
		mutate(raster_val = 1)		
	beach_mangrove_for <- for_form_sf %>%
		dplyr::filter(MAJOR_FORE == "Mangrove Forest" | MAJOR_FORE == "Beach Forest")

	for_form_red <- for_form_sf %>%
		dplyr::filter(!(MAJOR_FORE == "Mangrove Forest" | MAJOR_FORE == "Beach Forest"))
		
	# ----------------------	
	# Generate empty stack to hole outputs.	
	for_form_stack <- stack()
	
	# ----------------------	
	# For each polygon of forest formation, convert to sp object and then to raster.
	for(i in 1:nrow(for_form_red)){
		sf_tmp <- for_form_red[i,]
		sp_tmp <- as(sf_tmp, 'Spatial')
		r_tmp <- rasterize(sp_tmp, r_template, field = sp_tmp$raster_val)
		for_form_stack <- stack(for_form_stack, r_tmp)
		}
	
	for_form_feat_in <- for_form_stack
	
	#writeRaster(for_form_feat_in, "Prioritization/feature_inputs/forest_form_feat_in_w_pa_range.grd")
	
	

# =============================================================================	
###############################################################################