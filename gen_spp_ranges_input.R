###############################################################################
#  Connectivity through movement simulation.
#  September 18, 2017; last updated December 14, 2017
#  Script to generate spp ranges inputs for use in prioritizr package. In this packages' 
#   nomenclature, each spp range input is called a 'conservation feature' and is a 
#   RasterLayer within a RasterStack object. Cell values are 1 (spp range does occur) and 
#   0 (spp range does not occur).
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages.
# =============================================================================
library(raster)
library(sf)
library(dplyr)
library(ggplot2)



# =============================================================================
#  Load input data: for vertebrates, this is the output from adjusting the spp ranges based on 
#   elevation constraints, which was done using the script 'adjust_spp_ranges.R'. For butterflies
#   and plants, these are the SDM outputs from Sarah Scriven.
# =============================================================================

	# ----------------------
	#  Boundaries
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sabah_pa_sf.Rdata")

	
	# ----------------------
	#  Mammal spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_mammals_CR_EN_VU"
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_mammals_CR_EN_VU.Rdata")
	
	# ----------------------
	#  Amphibian spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_amphibians_CR_EN_VU"
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_amphibians_CR_EN_VU.Rdata")
	
	# ----------------------
	#  Bird spp ranges and associated spp data: 
	#   'sf' obejct named 'sf_adj_df_birds_CR_EN_VU"
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/spp_ranges_w_constraints/ranges_and_df/sf_adj_df_birds_CR_EN_VU.Rdata")

	# ----------------------
	#  Prepped butterfly SDM raster layers
	fly_end <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/end_fly_stack.grd")
	fly_non_end <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/non_end_fly_stack.grd")

	# ----------------------
	#  Prepped plant raster layers
	plants_end <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/all_end_plant_stack.grd")
	plants_non_end <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/all_non_end_plant_stack.grd")


# =============================================================================
#  Load vertebrate species range Excel data.
# =============================================================================	
	
	# ----------------------
	#  Data frame holding species information (spp name, status, etc.)
	#df_m <- read.csv("C:/Users/saraw/Documents/SEARRP_Analyses/processed_excel_data/spp_data_w_constraint_info/df_mammals_CR_EN_VU.csv")
	#df_a <- read.csv("C:/Users/saraw/Documents/SEARRP_Analyses/processed_excel_data/spp_data_w_constraint_info/df_amphibians_CR_EN_VU.csv")
	#df_b <- read.csv("C:/Users/saraw/Documents/SEARRP_Analyses/processed_excel_data/spp_data_w_constraint_info/df_birds_CR_EN_VU.csv")

	
	
# =============================================================================
#  Prep species ranges sf objects for rasterization.
# =============================================================================	
	
	# ----------------------
	#  Mammals
	sf_m_tmp <-  sf_adj_df_mammals_CR_EN_VU %>%
		mutate(ras_val = 1) %>%
		mutate(threat_wt = ifelse(RL_code == "VU", 2, 
			ifelse(RL_code == "EN", 3, 4))) %>%
		dplyr::filter(genus !=  "Dicerorhinus") %>% # remove rhino
		dplyr::select(id_no, binomial,  RL_code,  ras_val , threat_wt, order, family, 
			forest_dep, elev_min, elev_max) 
	#sf_m_tmp2 <- st_erase(sf_m_tmp1, st_buffer(sabah_pa_sf, 1))
	df_m_tmp <- sf_m_tmp
	st_geometry(df_m_tmp) <- NULL
	df_m <- droplevels(df_m_tmp)
	geom_m <- sf_m_tmp$geometry
	sf_m <- st_as_sf(cbind(df_m, geom_m))
	
	# ----------------------
	#  Amphibians
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
	#  Birds
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
	
	

# =============================================================================
#  Generate needed template for rasterizing species ranges.
# =============================================================================

	# ----------------------
	#  Load template raster.
	temp <- fly_end[[1]]
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Create existing PA raster.
	pa_sp <- as(sabah_pa_sf, 'Spatial')
	pa_r <- rasterize(pa_sp, r_template, field = pa_sp$locked_in_val) 
	
	
	
# =============================================================================
#  Rasterize each spp and add to a raster stack of input features using the raster template.
# =============================================================================

	# ----------------------
	#  Mammals
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
	#  Amphibians
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
	#  Birds
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
	#  All taxa in a single stack.
	sf_all <- rbind(sf_m, sf_a, sf_b)
	vert_all_tmp <- stack(feat_m, feat_a, feat_b)
	vert_all <- raster::mask(vert_all_tmp, pa_r, updateValue = 0, inverse = TRUE)
	vert_all[is.na(vert_all)] <- 0
	
	# ----------------------
	#  Determine which species have range ONLY within EXISTING PA and remove.
	cellStats(vert_all, max) 
	#  Those with 0 should be removed: layer 52, 54, 63 (2 amphibians, 1 bird).
	vert_all <- dropLayer(vert_all, c(52, 54, 63))
	
	
# =============================================================================
#  Weighting input features for objective in optimization - VERTS
# =============================================================================
 
 	# ----------------------
	#  Generate weighting vector.
	#   First remove rows for species that were removed from raster stack above.
	df_a_rem <- df_a %>%
		dplyr::filter(binomial != "Philautus gunungensis") %>%
		dplyr::filter(binomial != "Philautus saueri")
	df_b_rem <- df_b %>%
		dplyr::filter(binomial != "Ducula pickeringii")
	vert_rep_weight <- c(df_m$threat_wt, df_a_rem$threat_wt, df_b_rem$threat_wt)
	
	# ----------------------
	#  Using repeated rasters for spp depending on weight.
	vert_feat_rep <- stack()
	for(i in 1:nlayers(vert_all)){
		rep_stack <- stack(replicate(vert_rep_weight[[i]], vert_all[[i]]))
		vert_feat_rep <- stack(vert_feat_rep, rep_stack)
		}

	# ----------------------
	#  Save outputs.
	#writeRaster(vert_feat_rep, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_feat_rep.grd")
	#writeRaster(vert_all, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_all.grd")
	#save(vert_rep_weight, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_rep_weight.Rdata")
	
	
# =============================================================================
#  Weighting input features for objective in optimization - BUTTERFLIES
# =============================================================================
  	
	# ----------------------
	#  Combine all butterfly layers (endemic and non-endemic)
 	fly_all <- stack(fly_end, fly_non_end)
	
 	# ----------------------
	#  Generate weighting vector.
	fly_rep_weight <- c(rep(2, nlayers(fly_end)), rep(1, nlayers(fly_non_end)))
	
	# ----------------------
	#  Using repeated rasters for spp depending on weight.
	fly_feat_rep <- stack()
	for(i in 1:nlayers(fly_all)){
		rep_stack <- stack(replicate(fly_rep_weight[[i]], fly_all[[i]]))
		fly_feat_rep <- stack(fly_feat_rep, rep_stack)
		}

	# ----------------------
	#  Save outputs.
	#writeRaster(fly_feat_rep, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_feat_rep.grd")
	#writeRaster(fly_all, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_all.grd")
	#save(fly_rep_weight, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_rep_weight.Rdata")
	
	
# =============================================================================
#  Weighting input features for objective in optimization - PLANTS
# =============================================================================
  	
	# ----------------------
	#  Determine which species have range ONLY within EXISTING PA and remove.
 	plants_end_rem_tmp <- plants_end
	plants_end_rem <- raster::mask(plants_end_rem_tmp, pa_r, updateValue = 0, inverse = TRUE)
	plants_end_rem[is.na(plants_end_rem)] <- 0
	lay_max <- cellStats(plants_end_rem, max) 
	#  Those with 0 should be removed: layer 
	rem_lay <- which(lay_max %in% 0)
	plants_end_rem <- dropLayer(plants_end_rem, rem_lay)
	
	# ----------------------
	#  Determine which species have range ONLY within EXISTING PA and remove.
 	plants_non_end_rem_tmp <- plants_non_end
	plants_non_end_rem <- raster::mask(plants_non_end_rem_tmp, pa_r, updateValue = 0, inverse = TRUE)
	plants_non_end_rem[is.na(plants_non_end_rem)] <- 0
	lay_max <- cellStats(plants_non_end_rem, max) 
	#  Those with 0 should be removed: layer 
	rem_lay <- which(lay_max %in% 0)
	plants_non_end_rem <- dropLayer(plants_non_end_rem, rem_lay)
	
	# ----------------------
	#  Combine all plant layers (endemic and non-endemic)
	plant_all <- stack(plants_end_rem, plants_non_end_rem)
	
	# ----------------------
	#  Generate weighting vector.
	plant_rep_weight <- c(rep(2, nlayers(plants_end_rem)), rep(1, nlayers(plants_non_end_rem)))
	
	# ----------------------
	#  Using repeated rasters for spp depending on weight.
	plants_feat_rep <- stack()
	for(i in 1:nlayers(plant_all)){
		rep_stack <- stack(replicate(plant_rep_weight[[i]], plant_all[[i]]))
		plants_feat_rep <- stack(plants_feat_rep, rep_stack)
		}
	
	# ----------------------
	#  Save outputs.
	#writeRaster(plants_feat_rep, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plants_feat_rep.grd")
	#writeRaster(plants_all, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plants_all.grd")
	#save(plant_rep_weight, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plant_rep_weight.Rdata")

	
	
# =============================================================================
#  spp richness
# =============================================================================		
	
	vert_all <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_all.grd")
	fly_all <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_all.grd")
	plant_all <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plant_all.grd")	
	
	vert_sum <- raster::calc(vert_all, sum, na.rm = TRUE)
	x_min <- vert_sum@data@min
	x_max <- vert_sum@data@max
	vert_sum_in <- raster::calc(vert_sum, range01)

	fly_sum <- raster::calc(fly_all, sum, na.rm = TRUE)
	x_min <- fly_sum@data@min
	x_max <- fly_sum@data@max
	fly_sum_in <- raster::calc(fly_sum, range01)
	
	plant_sum <- raster::calc(plant_all, sum, na.rm = TRUE)
	x_min <- plant_sum@data@min
	x_max <- plant_sum@data@max
	plant_sum_in <- raster::calc(plant_sum, range01)
	
	writeRaster(vert_sum_in, "C:/Users/saraw/Desktop/tmp/vert_sum_in.grd")
	writeRaster(fly_sum_in, "C:/Users/saraw/Desktop/tmp/fly_sum_in.grd")
	writeRaster(plant_sum_in, "C:/Users/saraw/Desktop/tmp/plant_sum_in.grd")
	
	
	
# =============================================================================
#  Plots for output
# =============================================================================	
	
	# ----------------------
	#  Boundaries
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sarawak_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_kali_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf.Rdata")
	
	# ----------------------
	#  Boundaries
	#   Loop over each spp (row)  within each taxon.
	for(s in 1:nrow(sf_b)){

		dat <- sf_b %>%
			slice(s)
		nam <- dat$binomial
		
		#  Plot
		range_p <- ggplot() +
			geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey80") +
			geom_sf(data = ssk_pa_near_sf, colour = "darkorange2", fill = "darkorange2") +
			geom_sf(data = dat, colour = "transparent", fill = "darkred", alpha = 0.3) +
			geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
			geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey70") +
			coord_sf(crs = st_crs(32650)) +
			xlab("Longitude") +
			ylab("Latitude") +
			xlim(315000, 755000) +
			ylim(455000, 815000) +
			theme_bw()
		
		#  Save to pdf
		ggsave(filename = paste("C:/Users/saraw/Desktop/SEARRP_Meetings_Materials/new_spp_plots/", 
			nam, ".jpeg", sep = ""), plot = range_p)
		}   


	
# =============================================================================	
###############################################################################