###############################################################################
#  Prioritization: baseline scenario with chunk species inputs.
#  October 22, 2018; last updated November 26, 2018
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
library(gurobi)
library(prioritizr)
library(ggplot2)


setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Load species ranges (1 layer per species).
	vert_feat_in <- stack("chunks/old_version/prior_vert_feat_r_w_pa_range.grd")
	fly_feat_in <- stack("chunks/old_version/prior_fly_feat_r_w_pa_range.grd")	
	plant_feat_in <- stack("chunks/old_version/prior_plant_feat_r_w_pa_range.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("feature_inputs/vert_rep_weight_w_pa_range.Rdata")
	load("feature_inputs/fly_rep_weight_w_pa_range.Rdata")
	load("feature_inputs/plant_rep_weight_w_pa_range.Rdata")
	
	# ----------------------
	#  Forest formations inputs - 	raster 13, 18 and 19 (in order) in the stack have no 
	#   area outside of existing PAs so remove for this analysis.
	forest_form_feat_in <- stack("feature_inputs/forest_form_feat_in_w_pa_range.grd")
	#forest_form_feat_in <- forest_form_feat_in[[-13]]
	#forest_form_feat_in <- forest_form_feat_in[[-18]]
	#forest_form_feat_in <- forest_form_feat_in[[-19]]
	forest_form_feat_in[is.na(forest_form_feat_in)] <- 0
	
	# ----------------------
	#  Set up problem for connectivity and carbon input layers
	elev_conn_feat_in <- raster("feature_inputs/elev_conn_feat_in_w_pa_range.grd")
	corr_conn_feat_in <- raster("feature_inputs/corr_conn_feat_in_w_pa_range.grd")
	acd_feat_in <- raster("feature_inputs/acd_feat_in_w_pa_range.grd")
	
	# ----------------------
	#  Create raster following template of study area 
	temp <- fly_feat_in[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	
	st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
	
	
	
# =============================================================================
#  Set up planning units and locked in areas.
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	pu_r <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 

	area_ha <- cellStats(pu_r, sum)
	
	# ----------------------
	#  Get area of tpas
	locked_in_pas <- tpa_r
	# locked_in_pas_tmp <- stack(tpa_r, pu_r)
	# locked_in_pas_sum <- raster::calc(locked_in_pas_tmp, sum, na.rm = TRUE)
	# locked_in_pas <- locked_in_pas_sum
	# locked_in_pas[locked_in_pas < 200] <- 0
	# locked_in_pas[locked_in_pas > 0] <- 1
	
	area_tpa_ha <- cellStats(locked_in_pas, sum) 
	
	prior_area_ha <- 410000 + area_tpa_ha
	# 2,069,500 ha without TPAs that are not forested
	# 2,272,200

	
# =============================================================================
#  Run prioritzation with no BLM constraint
# =============================================================================	
	

	
# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert <- solve(p_vert)
	#cellStats(s_vert, sum)
	#plot(s_vert)
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly <- solve(p_fly)
	#cellStats(s_fly, sum)
	#plot(s_fly)
	
	# ----------------------
	#  Set up first round for plants
	p_plant <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_plant <- solve(p_plant)
	#cellStats(s_plant, sum)
	#plot(s_plant)
	
	# ----------------------
	#  Set up first round for plants
	p_form <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.475) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_form <- solve(p_form)
	#cellStats(s_form, sum)
	#plot(s_form)
	

	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	

	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	p2 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.55) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2 <- solve(p2)
	cellStats(s2, sum) * 100
	plot(s2)

	
	# ----------------------
	#  First round processing
	s2_sp <- rasterToPolygons(s2, fun = function(x){x>0}, dissolve = TRUE)
	s2_sf <- st_as_sf(s2_sp)
	area_thresh <- units::set_units(100, km^2)
	s2_sf_drop <- s2_sf %>%
		drop_crumbs(area_thresh)
	s2_sf_fill <- s2_sf_drop %>%
		fill_holes(area_thresh)
	s2_new_sp <- as(s2_sf_fill, "Spatial")
	s2_new <- rasterize(s2_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s2_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	#locked_in_new <- s2
	locked_in_new <- locked_in_sum
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.2) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.0000000003, 0.5) %>%
		#add_boundary_penalties(0.00000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	cellStats(s3, sum)
	plot(s3)
	
	
	
	
	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	p2 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.55) %>%
		add_boundary_penalties(0.00000000001, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2 <- solve(p2)
	cellStats(s2, sum) * 100
	plot(s2)

	
	# ----------------------
	#  First round processing
	s2_sp <- rasterToPolygons(s2, fun = function(x){x>0}, dissolve = TRUE)
	s2_sf <- st_as_sf(s2_sp)
	area_thresh <- units::set_units(100, km^2)
	s2_sf_drop <- s2_sf %>%
		drop_crumbs(area_thresh)
	s2_sf_fill <- s2_sf_drop %>%
		fill_holes(area_thresh)
	s2_new_sp <- as(s2_sf_fill, "Spatial")
	s2_new <- rasterize(s2_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s2_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	#locked_in_new <- s2
	locked_in_new <- locked_in_sum
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.2) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.0000000003, 0.5) %>%
		add_boundary_penalties(0.00000000003, 0.5) %>%
		#add_boundary_penalties(0.00000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	cellStats(s3, sum)
	plot(s3)
	
	# ----------------------
	#  Second round processing
	s3_sp <- rasterToPolygons(s3, fun = function(x){x>0}, dissolve = TRUE)
	s3_sf <- st_as_sf(s3_sp)
	area_thresh <- units::set_units(100, km^2)
	s3_sf_drop <- s3_sf %>%
		drop_crumbs(area_thresh)
	 #s3_sf_fill <- s3_sf_drop %>2%
		#fill_holes(area_thresh)
	s3_new_sp <- as(s3_sf_drop, "Spatial")
	s3_new <- rasterize(s3_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s3_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
		
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	#locked_in_new <- s2
	locked_in_new <- locked_in_sum2
	
	# ----------------------
	#  Set up problem for round 3.
	p4 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.0000000005, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s4 <- solve(p4)
	cellStats(s4, sum)
	plot(s4)
	
	# ----------------------
	#  Secpnd round processing
	s4_sp <- rasterToPolygons(s4, fun = function(x){x>0}, dissolve = TRUE)
	s4_sf <- st_as_sf(s4_sp)
	area_thresh <- units::set_units(100, km^2)
	s4_sf_drop <- s4_sf %>%
		drop_crumbs(area_thresh)
	 #s4_sf_fill <- s4_sf_drop %>%
		#fill_holes(area_thresh)
	s4_new_sp <- as(s4_sf_drop, "Spatial")
	s4_new <- rasterize(s4_new_sp, pu_r)
	final_selection_stack_tmp <- stack(locked_in_pas, s4_new)
	final_selection <- raster::calc(final_selection_stack_tmp, sum, na.rm = TRUE)
	final_selection[final_selection > 0] <- 1
	
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen1_w_pa_r <- s3
	writeRaster(scen1_w_pa_r, "scenario_outputs/scen1_w_pa_r.tif")

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
# =============================================================================
#  Run prioritzation with BLM constraint
# =============================================================================	
	


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert_blm <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		#add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000001, 0.5) %>%
		#add_boundary_penalties(0.00000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert_blm <- solve(p_vert_blm)
	cellStats(s_vert_blm, sum)
	plot(s_vert_blm)
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly_blm <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		#add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly_blm <- solve(p_fly_blm)
	cellStats(s_fly_blm, sum)
	plot(s_fly_blm)
	
	# ----------------------
	#  Set up first round for plants
	p_plant_blm <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		#add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 2.
	s_plant_blm <- solve(p_plant_blm)
	cellStats(s_plant_blm, sum)
	plot(s_plant_blm)
	
	# ----------------------
	#  Set up first round for plants
	p_form_blm <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.475) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		#add_boundary_penalties(0.00000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_form_blm <- solve(p_form_blm)
	cellStats(s_form_blm, sum)
	plot(s_form_blm)
	
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	#corr_conn_feat_in[corr_conn_feat_in < 0.35] <- 0
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in <- stack(vert_feat_in, fly_feat_in, plant_feat_in, s_form_blm,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	p2_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.2) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.0000000001, 0.5) %>%
		#add_boundary_penalties(0.00000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2_blm <- solve(p2_blm)
	cellStats(s2_blm, sum)
	plot(s2_blm)
	
	# ----------------------
	#  First round processing
	s2_blm_sp <- rasterToPolygons(s2_blm, fun = function(x){x>0}, dissolve = TRUE)
	s2_blm_sf <- st_as_sf(s2_blm_sp)
	area_thresh <- units::set_units(1000, ha)
	s2_blm_sf_drop <- s2_blm_sf %>%
		drop_crumbs(area_thresh)
	s2_blm_sf_fill <- s2_blm_sf_drop %>%
		fill_holes(area_thresh)
	s2_blm_new_sp <- as(s2_blm_sf_fill, "Spatial")
	s2_blm_new <- rasterize(s2_blm_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s2_blm_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	#locked_in_new_blm <- s2_blm
	locked_in_new_blm <- locked_in_sum
	
	# ----------------------
	#  Set up problem for round 3.
	p3_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.2) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.00000000025, 0.5) %>%
		#add_boundary_penalties(0.00000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3_blm <- solve(p3_blm)
	cellStats(s3_blm, sum)
	plot(s3_blm)
	
	# ----------------------
	#  Second round processing
	s3_blm_sp <- rasterToPolygons(s3_blm, fun = function(x){x>0}, dissolve = TRUE)
	s3_blm_sf <- st_as_sf(s3_blm_sp)
	area_thresh <- units::set_units(1000, ha)
	s3_blm_sf_drop <- s3_blm_sf %>%
		drop_crumbs(area_thresh)
	 #s3_blm_sf_fill <- s3_blm_sf_drop %>2%
		#fill_holes(area_thresh)
	s3_blm_new_sp <- as(s3_blm_sf_drop, "Spatial")
	s3_blm_new <- rasterize(s3_blm_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s3_blm_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	#locked_in_new_blm <- s2_blm
	locked_in_new_blm <- locked_in_sum2
	
	# ----------------------
	#  Set up problem for round 3.
	p4_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.0000000005, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s4_blm <- solve(p4_blm)
	cellStats(s4_blm, sum)
	plot(s4_blm)
	
	# ----------------------
	#  Secpnd round processing
	s4_blm_sp <- rasterToPolygons(s4_blm, fun = function(x){x>0}, dissolve = TRUE)
	s4_blm_sf <- st_as_sf(s4_blm_sp)
	area_thresh <- units::set_units(1000, ha)
	s4_blm_sf_drop <- s4_blm_sf %>%
		drop_crumbs(area_thresh)
	 #s4_blm_sf_fill <- s4_blm_sf_drop %>%
		#fill_holes(area_thresh)
	s4_blm_new_sp <- as(s4_blm_sf_drop, "Spatial")
	s4_blm_new <- rasterize(s4_blm_new_sp, pu_r)
	
	final_selection_stack_tmp <- stack(locked_in_pas, s4_blm)
	final_selection <- raster::calc(final_selection_stack_tmp, sum, na.rm = TRUE)
	final_selection[final_selection > 0] <- 1
	cellStats(final_selection, sum)
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	# scen1_w_pa_blm_r <- s3_blm
	# writeRaster(scen1_w_pa_blm_r, "scenario_outputs/scen1_w_pa_blm_r.tif")

	
	
	
	
	## Key
	myColors <- c('palegreen4',  'firebrick', 'dodgerblue4')
	myKey <- list(text=list(lab = c("Planning Units", "Selected Units", "Existing Protected Area")),
				  rectangles=list(col = myColors, border = FALSE),
				  space='bottom', columns=3)

	p0 <- levelplot(pu_r, margin = FALSE, col.regions=c('palegreen4'))
	p1 <- levelplot(final_selection, margin = FALSE, col.regions=c('transparent', 'firebrick'), 
		colorkey=FALSE, key = myKey, par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
	p1 + as.layer(p0, under = TRUE) +  layer(sp.polygons(main_sabah_tpa_sp, fill='dodgerblue4', alpha=0.8))
	
	
	
	