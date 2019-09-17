###############################################################################
#  Prioritization: baseline scenario with chunk species inputs.
#  October 22, 2018; last updated November 13, 2018
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
	vert_feat_in <- stack("feature_inputs/vert_all.grd")
	fly_feat_in <- stack("feature_inputs/fly_all.grd")	
	plant_feat_in <- stack("feature_inputs/plant_all.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("feature_inputs/vert_rep_weight.Rdata")
	load("feature_inputs/fly_rep_weight.Rdata")
	load("feature_inputs/plant_rep_weight.Rdata")
	
	# ----------------------
	# Connectivity and carbon input layers
	elev_conn_feat_in <- raster("feature_inputs/elev_conn_feat_in.grd")
	corr_conn_feat_in <- raster("feature_inputs/corr_conn_feat_in.grd")
	acd_feat_in <- raster("feature_inputs/acd_feat_in.grd")
	
	# ----------------------
	#  Forest formations inputs - 	raster 13, 18 and 19 (in order) in the stack have no 
	#   area outside of existing PAs so remove for this analysis.
	forest_form_feat_in <- stack("feature_inputs/forest_form_feat_in.grd")
	forest_form_feat_in <- forest_form_feat_in[[-13]]
	forest_form_feat_in <- forest_form_feat_in[[-18]]
	forest_form_feat_in <- forest_form_feat_in[[-19]]
	forest_form_feat_in[is.na(forest_form_feat_in)] <- 0
	
	# ----------------------
	#  Create raster following template of study area 
	temp <- fly_feat_in[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	
# =============================================================================
#  Set up planning units WITHOUT existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	pu_r_tmp <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 
	pu_r <- raster::mask(pu_r_tmp, tpa_r, inverse = TRUE)
	
	area_ha <- cellStats(pu_r, sum)
	
	
# =============================================================================
#  Run prioritzation with no BLM constraint
# =============================================================================	
	

	
# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert <- solve(p_vert)
	cellStats(s_vert, sum)
	plot(s_vert)
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>% 
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly <- solve(p_fly)
	cellStats(s_fly, sum)
	plot(s_fly)
	
	# ----------------------
	#  Set up first round for plants
	p_plant <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_plant <- solve(p_plant)
	cellStats(s_plant, sum)
	plot(s_plant)
	
	# ----------------------
	#  Set up first round for forest formations
	p_form <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.25) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_form <- solve(p_form)
	cellStats(s_form, sum)
	plot(s_form)
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, 
		elev_conn_feat_in, corr_conn_feat_in, acd_feat_in)
	#writeRaster(all_feat_in, "all_feat_in_scen1.tif")
		
	# ----------------------
	#  Set up problem for round 2.
	p2 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.59) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2 <- solve(p2)

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new <- s2
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(410000) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)

	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen1_r <- s3
	writeRaster(scen1_r, "scenario_outputs/scen1_r.tif")

	
	
	
	
	
	
	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in <- stack(s_vert, s_fly, s_plant, s_form,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)
	
	# ----------------------
	#  Set up problem for round 2.
	p2_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
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
	#area_thresh <- units::set_units(100, ha)
	s2_blm_sf_drop <- s2_blm_sf %>%	
		drop_crumbs(area_thresh)
	s2_blm_sf_fill <- s2_blm_sf_drop %>%
		fill_holes(area_thresh)
	s2_blm_new_sp <- as(s2_blm_sf_fill, "Spatial")
	s2_blm_new <- rasterize(s2_blm_new_sp, pu_r)

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- s2_blm_new
	cellStats(locked_in_new_blm, sum)
	plot(locked_in_new_blm)
	
	# ----------------------
	#  Set up problem for round 3.
	p3_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.225) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.000000000001, 0.5) %>%
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
	#area_thresh <- units::set_units(100, ha)	
	s3_blm_sf_drop <- s3_blm_sf %>%
		drop_crumbs(area_thresh)
	s3_blm_new_sp <- as(s3_blm_sf_drop, "Spatial")
	s3_blm_new <- rasterize(s3_blm_new_sp, pu_r)
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- s3_blm_new
	cellStats(locked_in_new_blm, sum)
	plot(locked_in_new_blm)
	
	# ----------------------
	#  Set up problem for round 3.
	p4_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.25) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.00000000001, 0.5) %>%
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
	#area_thresh <- units::set_units(100, ha)
	s4_blm_sf_drop <- s4_blm_sf %>%
		drop_crumbs(area_thresh)
	s4_blm_new_sp <- as(s4_blm_sf_drop, "Spatial")
	s4_blm_new <- rasterize(s4_blm_new_sp, pu_r)
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- s4_blm_new
	cellStats(locked_in_new_blm, sum)
	plot(locked_in_new_blm)
	
 	# ----------------------
	#  Set up problem for round 3.
	p5_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.275) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s5_blm <- solve(p5_blm)
	cellStats(s5_blm, sum)
	plot(s5_blm)
	
	# ----------------------
	#  Secpnd round processing
	s5_blm_sp <- rasterToPolygons(s5_blm, fun = function(x){x>0}, dissolve = TRUE)
	s5_blm_sf <- st_as_sf(s5_blm_sp)
	area_thresh <- units::set_units(1000, ha)
	#area_thresh <- units::set_units(500, ha)
	s5_blm_sf_drop <- s5_blm_sf %>%
		drop_crumbs(area_thresh)
	s5_blm_new_sp <- as(s5_blm_sf_drop, "Spatial")
	s5_blm_new <- rasterize(s5_blm_new_sp, pu_r)
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- s5_blm_new
	cellStats(locked_in_new_blm, sum)
	plot(locked_in_new_blm)
	
	# ----------------------
	#  Finalize raster output
	final_selection <- locked_in_new_blm
	final_selection[final_selection > 0] <- 1
	cellStats(final_selection, sum)
	
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen1_w_pa_blm_r <- final_selection
	writeRaster(scen1_w_pa_blm_r, "scenario_outputs/11_26_18/scen1_w_pa_blm_r.tif")

	
	
	
	
	## Key
	myColors <- c('palegreen4',  'firebrick', 'dodgerblue4')
	myKey <- list(text=list(lab = c("Planning Units", "Selected Units", "Existing Protected Area")),
				  rectangles=list(col = myColors, border = FALSE),
				  space='bottom', columns=3)

	p0 <- levelplot(pu_r, margin = FALSE, col.regions=c('palegreen4'))
	p1 <- levelplot(locked_in_new_blm, margin = FALSE, col.regions=c('firebrick'), 
		colorkey=FALSE, key = myKey, par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
	#p2 <- levelplot(locked_in_pas, margin = FALSE, col.regions
	p1 + as.layer(p0, under = TRUE) +  layer(sp.polygons(main_sabah_tpa_sp, fill='dodgerblue4', alpha=0.8))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
# =============================================================================
#  Run prioritzation with BLM constraint
# =============================================================================	
	


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert_blm <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>%
		add_feature_weights(vert_rep_weight) %>% 
		#add_boundary_penalties(0.000001, 0.5) %>%
		add_boundary_penalties(0.00000025, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert_blm <- solve(p_vert_blm)
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly_blm <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>% 
		add_feature_weights(fly_rep_weight) %>% 
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly_blm <- solve(p_fly_blm)
	
	# ----------------------
	#  Set up first round for plants
	p_plant_blm <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.175) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_plant_blm <- solve(p_plant_blm)
	
	# ----------------------
	#  Set up first round for forest formations
	p_form_blm <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.176) %>%
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_form_blm <- solve(p_form_blm)
	
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in <- stack(s_vert, s_fly, s_plant, s_form,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	p2_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(.625) %>%
		#add_boundary_penalties(0.0000000005, 0.5) %>%
		add_boundary_penalties(0.00000000025, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 2.
	s2_blm <- solve(p2_blm)
	cellStats(s2_blm, sum)
	plot(s2_blm)
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- s2_blm
	
	# ----------------------
	#  Set up problem for round 3.
	p3_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_utility_objective(410000) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3_blm <- solve(p3_blm)
	cellStats(s3_blm, sum)
	plot(s3_blm)
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen1_blm_r <- s3_blm
	writeRaster(scen1_blm_r, "scenario_outputs/scen1_blm_r.tif")

	