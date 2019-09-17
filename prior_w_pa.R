	
	

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
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert_blm <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert_blm <- solve(p_vert_blm)
	cellStats(s_vert_blm, sum)
	plot(s_vert_blm)
	
	s_vert_pa_m <- raster::mask(s_vert_blm, locked_in_pas)
	vert_ach <- cellStats(s_vert_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly_blm <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly_blm <- solve(p_fly_blm)
	cellStats(s_fly_blm, sum)
	plot(s_fly_blm)
	
	s_fly_pa_m <- raster::mask(s_fly_blm, locked_in_pas)
	fly_ach <- cellStats(s_fly_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for plants
	p_plant_blm <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 2.
	s_plant_blm <- solve(p_plant_blm)
	cellStats(s_plant_blm, sum)
	plot(s_plant_blm)
	
	s_plant_pa_m <- raster::mask(s_plant_blm, locked_in_pas)
	plant_ach <- cellStats(s_plant_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for plants
	p_form_blm <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_form_blm <- solve(p_form_blm)
	cellStats(s_form_blm, sum)
	plot(s_form_blm)

	s_form_pa_m <- raster::mask(s_form_blm, locked_in_pas)
	form_ach <- cellStats(s_form_pa_m, sum)/20695
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_form_blm,
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in_sm2)
	
	corr_pa_m <-  raster::mask(corr_conn_feat_in_sm2, locked_in_pas) 
	corr_pa_ach <- cellStats(corr_pa_m, sum)/cellStats(corr_conn_feat_in_sm2, sum)
	acd_pa_m <-  raster::mask(acd_feat_in, locked_in_pas) 
	acd_pa_ach <- cellStats(acd_pa_m, sum)/cellStats(acd_feat_in, sum)
	
	# ----------------------
	#  Set up problem for round 2.
	p_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(c(rep(0.95, 4), 0.2, 0.6, 0.95)) %>%
		add_boundary_penalties(0.000000001, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_blm <- solve(p_blm)
	cellStats(s_blm, sum)
	plot(s_blm)
	
	# ----------------------
	#  First round processing
	s_sp <- rasterToPolygons(s_blm, fun = function(x){x>0}, dissolve = TRUE)
	s_sf <- st_as_sf(s_sp)
	area_thresh <- units::set_units(1000, ha)
	s_sf_drop <- s_sf %>%
		drop_crumbs(area_thresh)
	s_sf_fill <- s_sf_drop %>%
		fill_holes(area_thresh)
	s_new_sp <- as(s_sf_fill, "Spatial")
	s_new <- rasterize(s_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new <- locked_in_sum
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	tst <- stack(s_vert, s_fly, s_plant, s_form, 
		elev_conn_feat_in, acd_feat_in, conn)
	
	
	p2 <- levelplot(tst, layout=c(4, 2), margin = FALSE, col.regions=c('transparent', 'dodgerblue4'),
		xlab="", ylab="", scales=list(draw=FALSE), colorkey=FALSE, key = myKey,
		names.attr = c("Prioritized Verts", "Prioritized Butterflies", "Prioritized Plants", "Prioritized Forest Forms",
			"Raw Elev Conn", "Raw ACD", "Raw Corridors"), 
		par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
	p2 + layer(sp.polygons(main_sabah_tpa_sp, lwd=0.8, col = "transparent", fill = "red", alpha = 0.6))
	
	
	
	
	
########################################################################################
	
	
	
	
	
	
	

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
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up first round for vertebrates
	p_vert_blm <- problem(x = pu_r, features = vert_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(vert_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve
	s_vert_blm <- solve(p_vert_blm)
	cellStats(s_vert_blm, sum)
	plot(s_vert_blm)
	
	s_vert_pa_m <- raster::mask(s_vert_blm, locked_in_pas)
	vert_ach <- cellStats(s_vert_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for butterflies
	p_fly_blm <- problem(x = pu_r, features = fly_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(fly_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s_fly_blm <- solve(p_fly_blm)
	cellStats(s_fly_blm, sum)
	plot(s_fly_blm)
	
	s_fly_pa_m <- raster::mask(s_fly_blm, locked_in_pas)
	fly_ach <- cellStats(s_fly_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for plants
	p_plant_blm <- problem(x = pu_r, features = plant_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.435) %>%
		add_feature_weights(plant_rep_weight_w_pa_range) %>% 
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 2.
	s_plant_blm <- solve(p_plant_blm)
	cellStats(s_plant_blm, sum)
	plot(s_plant_blm)
	
	s_plant_pa_m <- raster::mask(s_plant_blm, locked_in_pas)
	plant_ach <- cellStats(s_plant_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for plants
	p_form_blm <- problem(x = pu_r, features = forest_form_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_form_blm <- solve(p_form_blm)
	cellStats(s_form_blm, sum)
	plot(s_form_blm)

	s_form_pa_m <- raster::mask(s_form_blm, locked_in_pas)
	form_ach <- cellStats(s_form_pa_m, sum)/20695
	
	
	
	# ----------------------
	#  Set up first round for elevational connectivity
	p_cond_blm <- problem(x = pu_r, features = elev_conn_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.00000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_cond_blm <- solve(p_cond_blm)
	cellStats(s_cond_blm, sum)
	plot(s_cond_blm)
	
	s_cond_pa_m <- raster::mask(s_cond_blm, locked_in_pas)
	cond_ach <- cellStats(s_cond_pa_m, sum)/20695
	
	# ----------------------
	#  Set up first round for plants
	p_corr_blm <- problem(x = pu_r, features = corr_conn_feat_in_sm5) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.7) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_corr_blm <- solve(p_corr_blm)
	cellStats(s_corr_blm, sum)
	plot(s_corr_blm)

	# ----------------------
	#  Set up first round for plants
	p_acd_blm <- problem(x = pu_r, features = acd_feat_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		#add_boundary_penalties(0.000000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s_acd_blm <- solve(p_acd_blm)
	cellStats(s_acd_blm, sum)
	plot(s_acd_blm)
	
	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_form_blm,
		s_cond_blm, s_acd_blm, s_corr_blm)
	
	
	
	# ----------------------
	#  Set up problem for round 2.
	p2_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2_blm <- solve(p2_blm)
	cellStats(s2_blm, sum)
	plot(s2_blm)
	plot(locked_in_pas, add = T)
	

	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- s2_blm
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
		

		
		
		
		
	#s2_rem <- raster::mask(s2_blm, locked_in_pas, inverse = TRUE)
		
	# ----------------------
	#  First round processing
	s2_sp <- rasterToPolygons(s2_blm, fun = function(x){x>0}, dissolve = TRUE)
	s2_sf <- st_as_sf(s2_sp)
	area_thresh <- units::set_units(1000, ha)
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
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)
	
	
	
	
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	cellStats(s3, sum)
	plot(s3)
	plot(locked_in_pas, add = T)
	
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- s3
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	





	
	#s3_rem <- raster::mask(s3, locked_in_pas, inverse = TRUE)
	
	# ----------------------
	#  Second round processing
	s3_sp <- rasterToPolygons(s3, fun = function(x){x>0}, dissolve = TRUE)
	s3_sf <- st_as_sf(s3_sp)
	area_thresh <- units::set_units(1000, ha)
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
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)

	
	
	
	
	
	# ----------------------
	#  Set up problem for round 3.
	p4b <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		#add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.9) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s4b <- solve(p4b)
	cellStats(s4b, sum)
	plot(s4b)
	plot(locked_in_pas, add = T)
	
	s4b_sp <- rasterToPolygons(s4b, fun = function(x){x>0}, dissolve = TRUE)
	s4b_sf <- st_as_sf(s4b_sp)
	area_thresh <- units::set_units(5000, ha)
	s4b_sf_drop <- s4b_sf %>%
		drop_crumbs(area_thresh)
		
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- s4b
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	
	
	
	
	
	# ----------------------
	#  Set up problem for round 3.
	p4 <- problem(x = pu_r, features = all_feat_blm_in) %>%
		#add_max_utility_objective(prior_area_ha) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.925) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.0000000005, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s4 <- solve(p4)
	cellStats(s4, sum)
	plot(s4)
	plot(locked_in_pas, add = T)
		
		
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- s4
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
		
		
		
	s4_rem <- raster::mask(s4, locked_in_pas, inverse = TRUE)
	
	# ----------------------
	#  Second round processing
	s4_sp <- rasterToPolygons(s4_rem, fun = function(x){x>0}, dissolve = TRUE)
	s4_sf <- st_as_sf(s4_sp)
	area_thresh <- units::set_units(5000, ha)
	s4_sf_drop <- s4_sf %>%
		drop_crumbs(area_thresh)
	 #s4_sf_fill <- s4_sf_drop %>2%
		#fill_holes(area_thresh)
	s4_new_sp <- as(s4_sf_drop, "Spatial")
	s4_new <- rasterize(s4_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s4_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
		
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new <- locked_in_sum2
	cellStats(locked_in_new, sum)
	plot(locked_in_new)
	plot(locked_in_pas, add = T)

		

	
	myColors <- c('darkred',  'dodgerblue4')
	myKey <- list(text=list(lab = c("Selected Units", "PAs")),
		  rectangles=list(col = myColors, border = FALSE),
		  space='bottom', columns=2)
		
	p1 <- levelplot(s4b, margin = FALSE, col.regions=c('transparent', 'darkred'),
		xlab="", ylab="", scales=list(draw=FALSE), colorkey=FALSE, key = myKey,
		par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))	
	p0 <- levelplot(locked_in_pas, margin = FALSE, col.regions = rep("dodgerblue3",
		each = 8), xlab="", ylab="", scales=list(draw=FALSE), colorkey = FALSE)
	p1 + as.layer(p0, under = FALSE) + layer(sp.polygons(main_sabah_sp))
	
	
	p2 <- levelplot(tst, layout=c(4, 2), margin = FALSE, col.regions=c('transparent', 'dodgerblue4'),
		xlab="", ylab="", scales=list(draw=FALSE), colorkey=FALSE, key = myKey,
		names.attr = c("Prioritized Verts", "Prioritized Butterflies", "Prioritized Plants", "Prioritized Forest Forms",
			"Prioritized Elev Conn", "Prioritized ACD", "Prioritized Corridors"), 
		par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
	p2 + layer(sp.polygons(main_sabah_tpa_sp, lwd=0.8, col = "transparent", fill = "red", alpha = 0.6))
	
	
	

	p <- ggplot() +
		geom_sf(data = main_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = s4b_sf, colour = "darkred", fill = "darkred") +
		geom_sf(data = main_sabah_tpa, fill = "goldenrod3", colour = "transparent", alpha = 0.75) +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	p
	
	
	