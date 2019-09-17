

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
library(smoothr)


setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Load species ranges (1 layer per species).
	vert_feat_in <- stack("feature_inputs/vert_all_w_pa_range.grd")
	fly_feat_in <- stack("feature_inputs/fly_all_w_pa_range.grd")	
	plant_feat_in <- stack("feature_inputs/plant_all_w_pa_range.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("feature_inputs/vert_rep_weight_w_pa_range.Rdata")
	load("feature_inputs/fly_rep_weight_w_pa_range.Rdata")
	load("feature_inputs/plant_rep_weight_w_pa_range.Rdata")
	
	# ----------------------
	#  Forest formations inputs - 	raster 13, 18 and 19 (in order) in the stack have no 
	#   area outside of existing PAs so remove for this analysis.
	forest_form_feat_in <- stack("feature_inputs/forest_form_feat_in_w_pa_range.grd")
	forest_form_feat_in <- forest_form_feat_in[[-13]]
	forest_form_feat_in <- forest_form_feat_in[[-18]]
	forest_form_feat_in <- forest_form_feat_in[[-19]]
	forest_form_feat_in[is.na(forest_form_feat_in)] <- 0
	
	# ----------------------
	#  Set up problem for connectivity and carbon input layers
	elev_conn_feat_in <- raster("feature_inputs/elev_conn_feat_in_w_pa_range.grd")
	corr_conn_feat_in <- raster("feature_inputs/corr_conn_feat_in.grd")
	acd_feat_in <- raster("feature_inputs/acd_feat_in_w_pa_range.grd")
	
	# ----------------------
	#  Planning unit grids
	pu_r <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 
	area_ha <- cellStats(pu_r, sum)
	
	# ----------------------
	#  Get area of tpas
	locked_in_pas <- tpa_r
	locked_in_pas_tmp <- stack(tpa_r, pu_r)
	locked_in_pas_sum <- raster::calc(locked_in_pas_tmp, sum, na.rm = TRUE)
	locked_in_pas <- locked_in_pas_sum
	locked_in_pas[locked_in_pas < 200] <- 0
	locked_in_pas[locked_in_pas > 0] <- 100
	locked_in_pas[locked_in_pas == 0] <- NA
	area_tpa_ha <- cellStats(locked_in_pas, sum) 
	
	prior_area_ha <- 410000 + area_tpa_ha
	# 2,069,500 ha without TPAs that are not forested
	# 2,272,200
	
	
# =============================================================================
#  Set up and solve first round prioritization for species ranges, forest 
#   foramtions, ACD and connectivity layers for whole of Sabah
# =============================================================================

	# ----------------------
	#  Vertebrates
	p_vert <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	s_vert <- solve(p_vert)
	cellStats(s_vert, sum)

	#plot(s_vert)
	#s_vert_pa_m <- raster::mask(s_vert, locked_in_pas)
	#vert_ach <- cellStats(s_vert_pa_m, sum)/20695
	
	
	# ----------------------
	#  Butterflies
	p_fly <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	s_fly <- solve(p_fly)
	cellStats(s_fly, sum)
	
	#plot(s_fly)
	#s_fly_pa_m <- raster::mask(s_fly, locked_in_pas)
	#fly_ach <- cellStats(s_fly_pa_m, sum)/20695
	
	
	# ----------------------
	#  Plants
	p_plant <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	s_plant <- solve(p_plant)
	cellStats(s_plant, sum)
	
	#plot(s_plant)
	#s_plant_pa_m <- raster::mask(s_plant, locked_in_pas)
	#plant_ach <- cellStats(s_plant_pa_m, sum)/20695
	
	
	# ----------------------
	#  Forest formations
	p_form <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	s_form <- solve(p_form)
	cellStats(s_form, sum)
	
	#plot(s_form)
	#s_form_pa_m <- raster::mask(s_form, locked_in_pas)
	#form_ach <- cellStats(s_form_pa_m, sum)/20695
	
	
	# ----------------------
	#  Elevational connectivity
	p_cond <- problem(x = pu_r, features = elev_conn_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.31) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	s_cond <- solve(p_cond)
	cellStats(s_cond, sum)
	
	#plot(s_cond)


	# ----------------------
	#  Corridors
	p_corr <- problem(x = pu_r, features = corr_conn_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.999) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	s_corr <- solve(p_corr)
	cellStats(s_corr, sum)

	#plot(s_corr)

	
	# ----------------------
	#  ACD
	p_acd <- problem(x = pu_r, features = acd_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.5) %>%
		add_relative_targets(0.525) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	s_acd <- solve(p_acd)
	cellStats(s_acd, sum)

	#plot(s_acd)

	
	
	#writeRaster(s_vert, "s_vert.tif")
	#writeRaster(s_fly, "s_fly.tif")
	#writeRaster(s_plant, "s_plant.tif")
	#writeRaster(s_form, "s_form.tif")
	#writeRaster(s_cond, "s_cond.tif")
	#writeRaster(s_corr, "s_corr.tif")
	#writeRaster(s_acd, "s_acd.tif")
	
	

# =============================================================================
#  Set up and solve final prioritization using prioritized inputs from each category
# =============================================================================
	
	s_vert <-  raster("first_round_outputs/s_vert.tif")
	s_fly <-  raster("first_round_outputs/s_fly.tif")
	s_plant <-  raster("first_round_outputs/s_plant.tif")
	s_form <-  raster("first_round_outputs/s_form.tif")
	s_cond <-  raster("first_round_outputs/s_cond.tif")
	s_acd <-  raster("first_round_outputs/s_acd.tif")
	s_corr <-  raster("first_round_outputs/s_corr.tif")

	
	
	
	# ----------------------
	# Stack all prioritized areas for each category of input.
	#all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_cond, s_acd, s_corr)
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_acd) # for no connectivity outputs
	
	# ----------------------
	#  Set up problem with all prioritized inputs
	p_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		#add_relative_targets(0.9) %>% # for no connectivity outputs
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve.
	s_all <- solve(p_all)
	cellStats(s_all, sum)
	plot(s_all)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))

	sol <- s_all
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol

	# ----------------------
	#  Process to remove tiny "speckled" selections
	s_all_sp <- rasterToPolygons(s_all, fun = function(x){x>0}, dissolve = TRUE)
	s_all_sf <- st_as_sf(s_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s_all_sf_drop <- s_all_sf %>%
		drop_crumbs(area_thresh)
	s_all_sf_fill <- s_all_sf_drop %>%
		fill_holes(area_thresh)
	s_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s_all_new <- rasterize(s_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in areas (previous solution) for next round of selection.
	locked_in_new <- locked_in_sum
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Select additional area to fulfill budget with increased BLM
	p2_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		#add_relative_targets(0.9) %>% # for no connectivity outputs
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s2_all <- solve(p2_all)
	cellStats(s2_all, sum)
	plot(s2_all)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))

	sol <- s2_all
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	# ----------------------
	#  Process to remove tiny "speckled" selections again.
	s2_all_sp <- rasterToPolygons(s2_all, fun = function(x){x>0}, dissolve = TRUE)
	s2_all_sf <- st_as_sf(s2_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s2_all_sf_drop <- s2_all_sf %>%
		drop_crumbs(area_thresh)
	s2_all_new_sp <- as(s2_all_sf_drop, "Spatial")
	s2_all_new <- rasterize(s2_all_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s2_all_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
		
	# ----------------------
	#  Set up locked in areas (previous solution) for final round of selection
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)

	# ----------------------
	#  Select additional area to fulfill final remaining land area. 
	p3_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		#add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.9) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve.
	s3_all <- solve(p3_all)
	cellStats(s3_all, sum)
	#plot(s3_all)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Rename final output.
	final_sol <- s3_all
	no_conn_sol <- final_sol # for no connectivity outputs
	
	
	
	
	# ----------------------
	#  Check final coverage of each category of input
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_in))

	sol <- final_sol
	
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_prior_input_sol[i] <- cov
		}
	
	prop_cov_prior_input_sol
	
	

	
# =============================================================================
#  Sequentially remove each input and resolve.
# =============================================================================
	
	
	# sol_no_corr <- s3_all
	# sol_no_acd <- s3_all
	# sol_no_cond <- s3_all	
	
	# sol_no_form <- s3_all
	# sol_no_plant <- s3_all
	# sol_no_fly <- s3_all
	# sol_no_vert <- s3_all
	
	# ----------------------
	#  Save output - all input solution and sequential removal in a raster stack
	all_sol_stack <- stack(final_sol, sol_no_vert, sol_no_fly, sol_no_plant, sol_no_form, 
		sol_no_cond, sol_no_acd, sol_no_corr)
	# writeRaster(all_sol_stack, file = "all_sol_stack.grd")
	
	# ----------------------
	#  Save output  from no connectivity layers output
	# writeRaster(no_conn_sol, file = "no_conn_sol.tif")
	
	

# =============================================================================	
###############################################################################
	
	
	
	
	
	
	
	
	
	
	
	

# =============================================================================
#  Set up and solve final prioritization using prioritized inputs from each category
# =============================================================================
	
	# ----------------------
	# Stack all prioritized areas for each category of input.
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_acd)
	
	
	# ----------------------
	#  Set up problem with all prioritized inputs
	p_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_connectivity_penalties(1, data = d_matrix) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
			
	# ----------------------
	#  Solve.
	s_all <- solve(p_all)
	cellStats(s_all, sum)
	plot(s_all)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))

	sol <- s_all
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol

	# ----------------------
	#  Process to remove tiny "speckled" selections
	s_all_sp <- rasterToPolygons(s_all, fun = function(x){x>0}, dissolve = TRUE)
	s_all_sf <- st_as_sf(s_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s_all_sf_drop <- s_all_sf %>%
		drop_crumbs(area_thresh)
	s_all_sf_fill <- s_all_sf_drop %>%
		fill_holes(area_thresh)
	s_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s_all_new <- rasterize(s_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in areas (previous solution) for next round of selection.
	locked_in_new <- locked_in_sum
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Select additional area to fulfill budget with increased BLM
	p2_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s2_all <- solve(p2_all)
	cellStats(s2_all, sum)
	plot(s2_all)
	plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))

	sol <- s2_all
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	# ----------------------
	#  Process to remove tiny "speckled" selections again.
	s2_all_sp <- rasterToPolygons(s2_all, fun = function(x){x>0}, dissolve = TRUE)
	s2_all_sf <- st_as_sf(s2_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s2_all_sf_drop <- s2_all_sf %>%
		drop_crumbs(area_thresh)
	s2_all_new_sp <- as(s2_all_sf_drop, "Spatial")
	s2_all_new <- rasterize(s2_all_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s2_all_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
		
	# ----------------------
	#  Set up locked in areas (previous solution) for final round of selection
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)

	# ----------------------
	#  Select additional area to fulfill final remaining land area. 
	p3_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		#add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.9) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve.
	s3_all <- solve(p3_all)
	cellStats(s3_all, sum)
	#plot(s3_all)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Rename final output.
	final_sol <- s3_all
	
	# ----------------------
	#  Check final coverage of each category of input
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_in))

	sol <- final_sol
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_in_m <- raster::mask(all_feat_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_prior_input_sol[i] <- cov
		}
	
	prop_cov_prior_input_sol
	
	
		p3 <- levelplot(final_sol, margin = FALSE, col.regions=c('transparent', 'dodgerblue4'),
			xlab="", ylab="", scales=list(draw=FALSE), colorkey=FALSE, key = myKey,
			par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
		p3 + layer(sp.polygons(main_sabah_tpa_sp, lwd=0.8, col = "transparent", fill = "red", alpha = 0.4)) + layer(sp.polygons(main_sabah_sp, lwd=0.8))
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	all_feat <- stack(vert_feat_in, fly_feat_in, plant_feat_in, forest_form_feat_in, elev_conn_feat_in,
		acd_feat_in, corr_conn_feat_in)

	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- as(st_read(
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_m <- raster::mask(all_feat, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat)){
		range_cov_feat_sum <- cellStats(all_feat[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	vert_cov <- prop_cov_all_feat_sol[1:81]
	fly_cov <- prop_cov_all_feat_sol[82:158]
	plant_cov <- prop_cov_all_feat_sol[159:555]
	form_cov <- prop_cov_all_feat_sol[556:574]
	elev_conn_cov <- prop_cov_all_feat_sol[575]
	acd_cov <- prop_cov_all_feat_sol[576]
	corr_conn_cov <- prop_cov_all_feat_sol[577]
	
					# # ----------------------
					# #  Combine output of proportion covered across all input layers
					# #   into single DF
					# all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
						# acd_cov, corr_conn_cov)
					# cats <- c(rep("VERT", length(vert_cov)), 
						# rep("FLY", length(fly_cov)), 
						# rep("PLANT", length(plant_cov)),
						# rep("FORM", length(form_cov)),
						# rep("ELEV", length(elev_conn_cov)),
						# rep("ACD", length(acd_cov)),
						# rep("CORR", length(corr_conn_cov)))
					# weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, plant_rep_weight_w_pa_range,
						# rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
					# all_cov_df <- as.data.frame(cats) %>%
						# cbind(all_cov, weights)
					# names(all_cov_df) <- c("Category", "Proportion", "Weight")	
					# all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
					# all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
					# all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
					# all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
					# all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
					# all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
					# all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
					# all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
					# all_cov_df <- all_cov_df %>%
						# dplyr::filter(!is.na(Proportion)) 
					
						
					# all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
						# "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))
					
						
					# # ----------------------
					# # Box plot - for current scenario
					# cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
						# geom_point(aes(colour = Category, shape = Weight)) +
						# scale_shape_manual(name = "Weight", values = c(1, 2, 8, 19, 15, 5)) +
						# #geom_boxplot(aes(colour = Category)) +
						# scale_colour_manual(name = "Features of Conservation Interest",
							# values = c("mediumorchid1", "goldenrod2", "darkgreen",
								# "darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
							# breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
							# labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
						# ylab("Proportion of Feature Protected \n") +
						# ylim(0, 1) +
						# geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
						# theme_bw()
					# cov_p


	
	
	# Remove rare plants
	plant_cov1 <- prop_cov_all_feat_sol[159:206]
	plant_cov2 <- prop_cov_all_feat_sol[512:529]
	plant_cov <- c(plant_cov1, plant_cov2)
	
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, rep(1, length(plant_cov)),
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
		
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))
	
	cat_cov <- as.data.frame(cbind(c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"), 
		prop_cov_prior_input_sol))
	colnames(cat_cov) <- c("Category", "Proportion")
	cat_cov$Proportion <- as.numeric(as.character(cat_cov$Proportion))
	
	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		geom_point(data = cat_cov, aes(x = Category, y = Proportion, colour = Category)) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p
		
		
		
		
		
		
	

	prop_cov_all_feat_tpa <- numeric(length = nlayers(all_feat))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- locked_in_pas
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_m <- raster::mask(all_feat, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat)){
		range_cov_feat_sum <- cellStats(all_feat[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_tpa[i] <- cov
		}
	
	add_cov <- prop_cov_all_feat_sol - prop_cov_all_feat_tpa
	
	vert_cov <- add_cov[1:81]
	fly_cov <- add_cov[82:158]
	plant_cov <- add_cov[159:555]
	form_cov <- add_cov[556:574]
	elev_conn_cov <- add_cov[575]
	acd_cov <- add_cov[576]
	corr_conn_cov <- add_cov[577]
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, plant_rep_weight_w_pa_range,
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
		
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))

		
	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		#geom_point(aes(colour = Category, shape = Weight)) +
		#scale_shape_manual(name = "Weight", values = c(1, 2, 8, 19, 15, 5)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p

	
	
	# Remove rare plants
	plant_cov1 <- add_cov[159:206]
	plant_cov2 <- add_cov[512:529]
	plant_cov <- c(plant_cov1, plant_cov2)
	
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, rep(1, length(plant_cov)),
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
		
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))
	
	
	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		#geom_point(aes(colour = Category, shape = Weight)) +
		#scale_shape_manual(name = "Weight", values = c(1, 2, 8, 19, 15, 5)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	p4b <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		#add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.9) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	