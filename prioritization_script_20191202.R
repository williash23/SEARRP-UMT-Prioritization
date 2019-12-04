

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
#  Load feature input data
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



# =============================================================================
#  Load study area data
# =============================================================================
	
	# ---------------------
	#  Planning unit grids
	pu_r <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 
	area_ha <- cellStats(pu_r, sum)
	


# =============================================================================
#  Calculate areas and budgets
# =============================================================================

	# ----------------------
	#  Get total area of existing protected areas (tpas)
	locked_in_pas <- tpa_r
	locked_in_pas_tmp <- stack(tpa_r, pu_r)
	locked_in_pas_sum <- raster::calc(locked_in_pas_tmp, sum, na.rm = TRUE)
	locked_in_pas <- locked_in_pas_sum
	locked_in_pas[locked_in_pas < 200] <- 0
	locked_in_pas[locked_in_pas > 0] <- 100
	locked_in_pas[locked_in_pas == 0] <- NA
	area_tpa_ha <- cellStats(locked_in_pas, sum) 
	
	# ----------------------
	#  Calculate area to be prioritized (ha) including area of existing protected areas
	prior_area_ha <- 410000 + area_tpa_ha
	# 2,069,500 ha without TPAs that are not forested
	# 2,272,200
	
	
	
# =============================================================================
#  Set up and solve first round of prioritization for species ranges, forest 
#   foramtions, ACD and connectivity layers for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up conersvation problem for vertebrate ranges
	#   Input features are the 81 vertebrate ranges
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target was chosen manually to maximize coverage across all inputs 
	#   Feature weights based on IUCN conservation status (loaded above)
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_vert <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for vertebrate ranges
	s_vert <- solve(p_vert)
	
	# ----------------------
	#  The following 3 lines were used to assess coverage of all inputs layers for manual adjustment
	#   of targets
	#plot(s_vert)
	#cellStats(s_vert, sum)
	#s_vert_pa_m <- raster::mask(s_vert, locked_in_pas)
	#vert_ach <- cellStats(s_vert_pa_m, sum)/20695
	
	
	
	# ----------------------
	#  Set up conersvation problem for invertebrate ranges
	#   Input features are the 77 invertebrate ranges
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target was chosen manually to maximize coverage across all inputs 
	#   Feature weights based on endemism (loaded above)
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_fly <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for invertebrate ranges
	s_fly <- solve(p_fly)

	# ----------------------
	#  The following 3 lines were used to assess coverage of all inputs layers for manual adjustment
	#   of targets
	#plot(s_fly)
	#cellStats(s_fly, sum)
	#s_fly_pa_m <- raster::mask(s_fly, locked_in_pas)
	#fly_ach <- cellStats(s_fly_pa_m, sum)/20695
	
	
	
	# ----------------------
	#  Set up conersvation problem for plant ranges
	#   Input features are the 397 plant ranges
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target was chosen manually to maximize coverage across all inputs 
	#   Feature weights based on endemism (loaded above)
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_plant <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for plant ranges
	s_plant <- solve(p_plant)
	
	# ----------------------
	#  The following 3 lines were used to assess coverage of all inputs layers for manual adjustment
	#   of targets
	#plot(s_plant)
	#cellStats(s_plant, sum)
	#s_plant_pa_m <- raster::mask(s_plant, locked_in_pas)
	#plant_ach <- cellStats(s_plant_pa_m, sum)/20695
	
	
	
	# ----------------------
	#  Set up conersvation problem for forest formation distributions
	#   Input features are the 19 forest formation distributions
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target was chosen manually to maximize coverage across all inputs 
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_form <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for forest formations
	s_form <- solve(p_form)

	# ----------------------
	#  The following 3 lines were used to assess coverage of all inputs layers for manual adjustment
	#   of targets
	#plot(s_form)
	#	cellStats(s_form, sum)
	#s_form_pa_m <- raster::mask(s_form, locked_in_pas)
	#form_ach <- cellStats(s_form_pa_m, sum)/20695
	
	
	
	# ----------------------
	#  Set up conersvation problem for elevational connectivity
	#   Input feature is single elevational connectivity layer (output from Condatis analysis)
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target adjusted manually to reach full land area buget (Prioritizing over a sinlge layer 
	#    so simply selecting area with highest values)
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_cond <- problem(x = pu_r, features = elev_conn_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.31) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for elevational connectivity
	s_cond <- solve(p_cond)
	
	# ----------------------
	#  The following 2 lines were used to assess total area of solution for manual adjustment
	#   of target
	#plot(s_cond)
	#cellStats(s_cond, sum)
	
	

	# ----------------------
	#  Set up conersvation problem for corridor connectivity
	#   Input feature is single corridors layer
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target adjusted manually to reach full land area buget (Prioritizing over a sinlge layer 
	#    so simply selecting area with highest values)
	#   Existing proctected areas (TPAs) locked into solution
	#   Include very small boundary length modifier penalty to promote contiguous corridors
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_corr <- problem(x = pu_r, features = corr_conn_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.999) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.0000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve conservation problem for corridor connectivity
	s_corr <- solve(p_corr)
	
	# ----------------------
	#  The following 2 lines were used to assess total area of solution for manual adjustment
	#   of target
	#plot(s_corr)
	#cellStats(s_corr, sum)

	
	
	# ----------------------
	#  Set up conersvation problem for above ground carbon (ACD)
	#   Input feature is single ACD layer
	#   Land area budget is 410,000 ha plus area of tpas
	#   Target adjusted manually to reach full land area buget (Prioritizing over a sinlge layer 
	#    so simply selecting area with highest values)
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p_acd <- problem(x = pu_r, features = acd_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.525) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve conservation problem for aboveground carbon (ACD)
	s_acd <- solve(p_acd)
	
	# ----------------------
	#  The following 2 lines were used to assess total area of solution for manual adjustment
	#   of target
	#plot(s_acd)
	#cellStats(s_acd, sum)

	

# =============================================================================
#  Save all outputs from first round of prioritization
# =============================================================================

	#writeRaster(s_vert, "s_vert.tif")
	#writeRaster(s_fly, "s_fly.tif")
	#writeRaster(s_plant, "s_plant.tif")
	#writeRaster(s_form, "s_form.tif")
	#writeRaster(s_cond, "s_cond.tif")
	#writeRaster(s_corr, "s_corr.tif")
	#writeRaster(s_acd, "s_acd.tif")
	
	

# =============================================================================
#  Set up and solve second round of prioritization for using outputs from first round of prioritization
#   Outputs of first round gerenated a single layer per category (7 layers total) that will be 
#   prioritized in second roud.
#  Perform *TWO* separate scenarios in second round: 
#   1. Using all 7 inputs
#   2. Exlcuding connectivity inputs (using 5 inputs)
#  Additionally, the second round of prioritization is performed in 3 iterations such that very small 
#   and isolated areas that were selected by the priorititzation were removed from the selection. 
#   The remaining selected area was then locked in, and the process was repeated two more times.
# =============================================================================

	# ----------------------
	#  Reload outputs from first round of prioritization if needed
	#s_vert <-  raster("first_round_outputs/s_vert.tif")
	#s_fly <-  raster("first_round_outputs/s_fly.tif")
	#s_plant <-  raster("first_round_outputs/s_plant.tif")
	#s_form <-  raster("first_round_outputs/s_form.tif")
	#s_cond <-  raster("first_round_outputs/s_cond.tif")
	#s_acd <-  raster("first_round_outputs/s_acd.tif")
	#s_corr <-  raster("first_round_outputs/s_corr.tif")



# =============================================================================
# FIRST SCENARIO OF SECOND ROUND OF PRIORITIZATION: ALL INPUTS
# =============================================================================

	# ----------------------
	#  Stack all prioritized areas for each category of input.
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_cond, s_acd, s_corr)
	
	# ----------------------
	#  Set up first iteration of conservation problem
	#   Input features are ALL outputs from first round of prioritization (7 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include very small boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p1_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve problem
	s1_all <- solve(p1_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s1_all, sum)
	plot(s1_all)
	plot(locked_in_pas, add = T)
	
	sol <- s1_all
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
	
	

	# ----------------------
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
	s1_all_sp <- rasterToPolygons(s1_all, fun = function(x){x>0}, dissolve = TRUE)
	s1_all_sf <- st_as_sf(s1_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s1_all_sf_drop <- s1_all_sf %>%
		drop_crumbs(area_thresh)
	s1_all_sf_fill <- s1_all_sf_drop %>%
		fill_holes(area_thresh)
	s1_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s1_all_new <- rasterize(s1_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s1_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	

	# ----------------------
	#  Set up second iteration of conservation problem: selecting additional area to fulfill land areas 
	#   budget with increased BLM
	#   Input features are ALL outputs from first round of prioritization (7 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include slightly larger boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p2_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s2_all <- solve(p2_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s2_all, sum)
	plot(s2_all)
	plot(locked_in_pas, add = T)
	
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
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
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
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)



	# ----------------------
	#  Set up third iteration of conservation problem: selecting additional area to fulfill land area
	#   budget 
	#   Input features are ALL outputs from first round of prioritization (7 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize utility objective  does not need target
	#   Target adjusted manually to reach full land area buget
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p3_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		add_locked_in_constraints(locked_in_new) %>%
				add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s3_all <- solve(p3_all)
	
	# ----------------------
	#  Check final coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s3_all, sum)
	plot(s3_all)
	plot(locked_in_pas, add = T)
	
	sol <- s3_all
	
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
#  Save output from first scenario second round of prioritization for ALL INPUTS
# =============================================================================
	
	# ----------------------
	#  Rename final output.
	final_sol <- s3_all
	
	# ----------------------
	#  Save output - all input solution and sequential removal in a raster stack
	# writeRaster(final_sol, file = "final_sol.tif")
	
	
	
	

# =============================================================================
# SECOND SCENARIO OF SECOND ROUND OF PRIORITIZATION: NO CONNECTIVITY INPUTS
#  Note that target is higher in this second prioritization because there are fewer inputs 
# =============================================================================

	# ----------------------
	#  Stack all prioritized areas for each category of input.
	all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_acd) #  Note no connectivity inputs included
	
	# ----------------------
	#  Set up first iteration of conservation problem
	#   Input features are species ranges, forest formations and ACD outputs from first round of prioritization
	#	  (5 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include very small boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p1_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.90) %>% 
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve problem
	s1_all <- solve(p1_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s1_all, sum)
	plot(s1_all)
	plot(locked_in_pas, add = T)
	
	sol <- s1_all
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
	
	

	# ----------------------
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
	s1_all_sp <- rasterToPolygons(s1_all, fun = function(x){x>0}, dissolve = TRUE)
	s1_all_sf <- st_as_sf(s1_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s1_all_sf_drop <- s1_all_sf %>%
		drop_crumbs(area_thresh)
	s1_all_sf_fill <- s1_all_sf_drop %>%
		fill_holes(area_thresh)
	s1_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s1_all_new <- rasterize(s1_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s1_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	

	# ----------------------
	#  Set up second iteration of conservation problem: selecting additional area to fulfill land areas 
	#   budget with increased BLM
	#   Input features are species ranges, forest formations and ACD outputs from first round of prioritization
	#	  (5 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include slightly larger boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p2_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.90) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s2_all <- solve(p2_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s2_all, sum)
	plot(s2_all)
	plot(locked_in_pas, add = T)
	
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
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
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
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)



	# ----------------------
	#  Set up third iteration of conservation problem: selecting additional area to fulfill land area
	#   budget 
	#   Input features are species ranges, forest formations and ACD outputs from first round of prioritization
	#	  (5 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize utility objective does not need target
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p3_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s3_all <- solve(p3_all)
	
	# ----------------------
	#  Check final coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s3_all, sum)
	plot(s3_all)
	plot(locked_in_pas, add = T)
	
	sol <- s3_all
	
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
#  Save output from second scenario from second round of prioritization EXCLUDING 
#   CONNECTIVITY INPUTS
# =============================================================================
	
	# ----------------------
	#  Rename final output.
	final_sol_no_conn <- s3_all
	
	# ----------------------
	#  Save output - all input solution and sequential removal in a raster stack
	# writeRaster(final_sol_no_conn, file = "final_sol_no_conn.tif")
	
	


	

# =============================================================================
#  Lastly, to perform the sensitivity analysis, repeat the first scenario of the second round of 
#   prioritization and sequentially remove ONE of the inputs. Repeat this with removal of each input.
# =============================================================================


	# ----------------------
	#  Stack all prioritized areas for each category of input and remove a single input sequentially
	#all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_cond, s_acd)
	#all_feat_in <- stack(s_vert, s_fly, s_plant, s_form, s_cond, s_corr)
	#all_feat_in <- stack(s_vert, s_fly, s_plant, s_form,  s_acd, s_corr)
	#all_feat_in <- stack(s_vert, s_fly, s_plant, s_cond, s_acd, s_corr)
	#all_feat_in <- stack(s_vert, s_fly, s_form, s_cond, s_acd, s_corr)
	#all_feat_in <- stack(s_vert, s_plant, s_form, s_cond, s_acd, s_corr)
	#all_feat_in <- stack(s_fly, s_plant, s_form, s_cond, s_acd, s_corr)
	
	# ----------------------
	#  Set up first iteration of conservation problem
	#   Input features are ALL BUT ONE of the outputs from first round of prioritization 
	#   (6 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include very small boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p1_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve problem
	s1_all <- solve(p1_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s1_all, sum)
	plot(s1_all)
	plot(locked_in_pas, add = T)
	
	sol <- s1_all
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
	
	

	# ----------------------
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
	s1_all_sp <- rasterToPolygons(s1_all, fun = function(x){x>0}, dissolve = TRUE)
	s1_all_sf <- st_as_sf(s1_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s1_all_sf_drop <- s1_all_sf %>%
		drop_crumbs(area_thresh)
	s1_all_sf_fill <- s1_all_sf_drop %>%
		fill_holes(area_thresh)
	s1_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s1_all_new <- rasterize(s1_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s1_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	

	# ----------------------
	#  Set up second iteration of conservation problem: selecting additional area to fulfill land areas 
	#   budget with increased BLM
	#   Input features are ALL BUT ONE of the outputs from first round of prioritization 
	#   (6 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize features objective
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Include slightly larger boundary length modifier penalty to promote "clumping"
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p2_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s2_all <- solve(p2_all)
	
	# ----------------------
	#  Check coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s2_all, sum)
	plot(s2_all)
	plot(locked_in_pas, add = T)
	
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
	#  Remove areas of prioritized selection that are very small (less than 1,000 ha) and isolated
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
	#  Make remaining area selected in prioritzation (minus the removed small, isolated areas) locked in 
	#   for next iteration of selection.
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)



	# ----------------------
	#  Set up third iteration of conservation problem: selecting additional area to fulfill land area
	#   budget 
	#   Input features are ALL BUT ONE of the outputs from first round of prioritization 
	#   (6 input layers total)
	#   Land area budget is 410,000 ha plus area of tpas - use maximize utility objective  does not need target
	#   Target adjusted manually to reach full land area buget 
	#   Existing proctected areas (TPAs) locked into solution
	#   Require bindary outputs (0/1 = no/yes for inclusion in solution)
	p3_all <- problem(x = pu_r, features = all_feat_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		add_locked_in_constraints(locked_in_new) %>%
				add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve problem
	s3_all <- solve(p3_all)
	
	# ----------------------
	#  Check final coverage of each input and assess eveness of coverage (if one input has much higher 
	#   coverage, adjust target) 
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_in))
	cellStats(s3_all, sum)
	plot(s3_all)
	plot(locked_in_pas, add = T)
	
	sol <- s3_all
	
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
	
	

	# ----------------------
	#  Rename solution according to which single input was removed
	# sol_no_corr <- s3_all
	# sol_no_acd <- s3_all
	# sol_no_cond <- s3_all	
	
	# sol_no_form <- s3_all
	# sol_no_plant <- s3_all
	# sol_no_fly <- s3_all
	# sol_no_vert <- s3_all
	
	

# =============================================================================
#  Save all outputs from sequentially removal in  a raster stack.
# =============================================================================
	
	
	# ----------------------
	#  Save output - solution with all 7 inputs and sequential removal outputs in a raster stack
	seq_rem_sol_stack <- stack(final_sol, sol_no_vert, sol_no_fly, sol_no_plant, sol_no_form, 
		sol_no_cond, sol_no_acd, sol_no_corr)
	# writeRaster(seq_rem_sol_stack, file = "seq_rem_sol_stack.tif")
	
