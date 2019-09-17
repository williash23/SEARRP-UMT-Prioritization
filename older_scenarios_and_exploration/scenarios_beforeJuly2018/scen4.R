###############################################################################
#  Prioritization scenario - 4
#  Sara Williams
###############################################################################


#  Conditions:
#     1. SFI/Idris:
#        - SFI- 50,000 ha
#        - Idris - 100,000 ha
#        - 6 individual inputs prioritized separately for each area then prioritized together for each
#          area for a total ha in each
#     2. Together, SFI, Idris and Deramakot became locked in for overall prioritization
#     3. Segaliud Lokan Forest Reserve locked in and added to locked in area for overall prioritization
#     4. Area between Crocker Range and Kinabalu Park is locked out   
#     5. Total area for overall prioritization 410,000 ha
#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
#          for 410,000 ha across all Sabah

#  Ran above for both:
#     1. No BLM constraints
#     2. BLM constraint

setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
# install.packages("C:/gurobi751/win64/R/gurobi_7.5-1.zip", repos = NULL)
library(gurobi)
#  devtools::install_github("prioritizr/priortizr") # can use official CRAN version
library(prioritizr)
library(ggplot2)



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Load species ranges (1 layer per species).
	vert_feat_in_single <- stack("feature_inputs/vert_all.grd")
	fly_feat_in_single <- stack("feature_inputs/fly_all.grd")	
	plant_feat_in_single <- stack("feature_inputs/plant_all.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("feature_inputs/vert_rep_weight.Rdata")
	load("feature_inputs/fly_rep_weight.Rdata")
	load("feature_inputs/plant_rep_weight.Rdata")
	
	# ----------------------
	#  Set up problem for connectivity and carbon input layers
	elev_conn_feat_r <- raster("feature_inputs/elev_conn_feat_in.grd")
	corr_feat_r <- raster("feature_inputs/corr_feat_in.grd")
	acd_feat_r <- raster("feature_inputs/acd_feat_in.grd")
	
	# ----------------------
	#  Create raster following template of study area 
	fly_feat_in_single <- stack("feature_inputs/fly_all.grd")
	temp <- fly_feat_in_single[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	
	
# =============================================================================
#  Set up planning units WITHOUT existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	load(file = "planning_unit_grids/sa_grid.Rdata")
	load(file = "planning_unit_grids/sfi_grid.Rdata")
	load(file = "planning_unit_grids/idris_grid.Rdata")
	load(file = "planning_unit_grids/deram_grid.Rdata")
	load(file = "planning_unit_grids/crock_kina_grid.Rdata")
	load(file = "planning_unit_grids/sega_grid.Rdata")
	

	const_cost_h <- 500
	pu_sf_tmp <- sa_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_sf_tmp$area_h_tmp <- as.numeric(pu_sf_tmp$area_h_tmp)
	pu_sf <- pu_sf_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100)
	pu_in <- as(pu_sf, "Spatial")
	
	pu_sfi_tmp <- sfi_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_sfi_tmp$area_h_tmp <- as.numeric(pu_sfi_tmp$area_h_tmp)
	pu_sfi <- pu_sfi_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100)
	pu_sfi_in <- as(pu_sfi, "Spatial")
	
	pu_idris_tmp <- idris_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_idris_tmp$area_h_tmp <- as.numeric(pu_idris_tmp$area_h_tmp)
	pu_idris <- pu_idris_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100)
	pu_idris_in <- as(pu_idris, "Spatial")
	
	pu_deram_tmp <- deram_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_deram_tmp$area_h_tmp <- as.numeric(pu_deram_tmp$area_h_tmp)
	pu_deram <- pu_deram_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100) %>%
		mutate(solution_1 = 1)
	
	pu_sega_tmp <- sega_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_sega_tmp$area_h_tmp <- as.numeric(pu_sega_tmp$area_h_tmp)
	pu_sega <- pu_sega_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100) %>%
		mutate(solution_1 = 1)
	

	
# =============================================================================
#  Set up initial locked in and locked out areas.
# =============================================================================
	
	# ----------------------
	#  Locked in
	locked_in_tmp <- rbind(pu_sega, pu_deram)
	
	# ----------------------
	#  Locked out
	locked_out <- as(crock_kina_grid, 'Spatial')


	
# =============================================================================
#  Run prioritzation with no BLM constraint
# =============================================================================	
	


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges within SFI area
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1_sfi <- problem(x = pu_sfi_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert1_sfi <- solve(p_vert1_sfi)
	s_vert1_sfi_sf <- st_as_sf(s_vert1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_vert1_sfi_area_h <- as.integer(sum(st_area(s_vert1_sfi_sf)) * 0.0001)
	s_vert1_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1_sfi <- problem(x = pu_sfi_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1_sfi <- solve(p_fly1_sfi)
	s_fly1_sfi_sf <- st_as_sf(s_fly1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_sfi_area_h <- as.integer(sum(st_area(s_fly1_sfi_sf)) * 0.0001)
	s_fly1_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for plants
	p_plant1_sfi <- problem(x = pu_sfi_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1_sfi <- solve(p_plant1_sfi)
	s_plant1_sfi_sf <- st_as_sf(s_plant1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_sfi_area_h <- as.integer(sum(st_area(s_plant1_sfi_sf)) * 0.0001)
	s_plant1_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1_sfi <- problem(x = pu_sfi_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.3) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1_sfi <- solve(p_cond1_sfi)
	s_cond1_sfi_sf <- st_as_sf(s_cond1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_sfi_area_h <- as.integer(sum(st_area(s_cond1_sfi_sf)) * 0.0001)
	s_cond1_sfi_area_h
	

	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1_sfi <- problem(x = pu_sfi_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.25) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1_sfi <- solve(p_acd1_sfi)
	s_acd1_sfi_sf <- st_as_sf(s_acd1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_sfi_area_h <- as.integer(sum(st_area(s_acd1_sfi_sf)) * 0.0001)
	s_acd1_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1_sfi <- problem(x = pu_sfi_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.45) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1_sfi <- solve(p_corr1_sfi)
	s_corr1_sfi_sf <- st_as_sf(s_corr1_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_sfi_area_h <- as.integer(sum(st_area(s_corr1_sfi_sf)) * 0.0001)
	s_corr1_sfi_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_sp_sfi <- as(s_vert1_sfi_sf, 'Spatial')
	s_fly_sp_sfi <- as(s_fly1_sfi_sf, 'Spatial')
	s_plant_sp_sfi <- as(s_plant1_sfi_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_feat_r_sfi <- rasterize(s_vert_sp_sfi, r_template, field = s_vert_sp_sfi$solution_1)
	prior_fly_feat_r_sfi <- rasterize(s_fly_sp_sfi, r_template, field = s_fly_sp_sfi$solution_1)
	prior_plant_feat_r_sfi <- rasterize(s_plant_sp_sfi, r_template, field = s_plant_sp_sfi$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_sfi_u <- st_union(s_cond1_sfi_sf)
	s_acd1_sfi_u <- st_union(s_acd1_sfi_sf)
	s_corr1_sfi_u <- st_union(s_corr1_sfi_sf)

	# ----------------------
	#  Convert to sp objects
	s_cond_sp_sfi <- as(s_cond1_sfi_u, 'Spatial')
	s_acd_sp_sfi <- as(s_acd1_sfi_u, 'Spatial')
	s_corr_sp_sfi <- as(s_corr1_sfi_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_feat_r_sfi <- raster::mask(elev_conn_feat_r, s_cond_sp_sfi)
	prior_acd_feat_r_sfi <- raster::mask(acd_feat_r, s_acd_sp_sfi)
	prior_corr_feat_r_sfi <- raster::mask(corr_feat_r, s_corr_sp_sfi)
	
	

# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together in SFI area.
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_in_sfi <- stack(prior_vert_feat_r_sfi, prior_fly_feat_r_sfi, prior_plant_feat_r_sfi, 
		prior_cond_feat_r_sfi, prior_corr_feat_r_sfi, prior_acd_feat_r_sfi)

	# ----------------------
	#  Set up problem
	equal_targets <- 0.45
	p1_sfi <- problem(x = pu_sfi_in, features = all_feat_in_sfi, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(equal_targets) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1_sfi <- solve(p1_sfi)
	s1_sfi_sf <- st_as_sf(s1_sfi) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_sfi_area_h <- as.integer(sum(st_area(s1_sfi_sf)) * 0.0001)
	s1_sfi_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_sfi1 <- as(s1_sfi_sf, 'Spatial')

	# ----------------------
	#  Set up problem for round 3.
	p2_sfi <- problem(x = pu_sfi_in, features = all_feat_in_sfi, cost_column = "area_h") %>%
		add_max_utility_objective(50000) %>%
		add_locked_in_constraints(locked_in_sfi1) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s2_sfi <- solve(p2_sfi)
	s2_sfi_sf <- st_as_sf(s2_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s2_sfi_area_h <- as.integer(sum(st_area(s2_sfi_sf)) * 0.0001)
	s2_sfi_area_h
	
	
	
# =============================================================================
#  STEP 3: Set up and solve problem for species ranges within Idris area
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1_idris <- problem(x = pu_idris_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.4) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert1_idris <- solve(p_vert1_idris)
	s_vert1_idris_sf <- st_as_sf(s_vert1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_vert1_idris_area_h <- as.integer(sum(st_area(s_vert1_idris_sf)) * 0.0001)
	s_vert1_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1_idris <- problem(x = pu_idris_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.45) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1_idris <- solve(p_fly1_idris)
	s_fly1_idris_sf <- st_as_sf(s_fly1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_idris_area_h <- as.integer(sum(st_area(s_fly1_idris_sf)) * 0.0001)
	s_fly1_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for plants
	p_plant1_idris <- problem(x = pu_idris_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.4) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1_idris <- solve(p_plant1_idris)
	s_plant1_idris_sf <- st_as_sf(s_plant1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_idris_area_h <- as.integer(sum(st_area(s_plant1_idris_sf)) * 0.0001)
	s_plant1_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1_idris <- problem(x = pu_idris_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.65) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1_idris <- solve(p_cond1_idris)
	s_cond1_idris_sf <- st_as_sf(s_cond1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_idris_area_h <- as.integer(sum(st_area(s_cond1_idris_sf)) * 0.0001)
	s_cond1_idris_area_h
	

	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1_idris <- problem(x = pu_idris_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.4) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1_idris <- solve(p_acd1_idris)
	s_acd1_idris_sf <- st_as_sf(s_acd1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_idris_area_h <- as.integer(sum(st_area(s_acd1_idris_sf)) * 0.0001)
	s_acd1_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1_idris <- problem(x = pu_idris_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_relative_targets(0.95) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1_idris <- solve(p_corr1_idris)
	s_corr1_idris_sf <- st_as_sf(s_corr1_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_idris_area_h <- as.integer(sum(st_area(s_corr1_idris_sf)) * 0.0001)
	s_corr1_idris_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_sp_idris <- as(s_vert1_idris_sf, 'Spatial')
	s_fly_sp_idris <- as(s_fly1_idris_sf, 'Spatial')
	s_plant_sp_idris <- as(s_plant1_idris_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_feat_r_idris <- rasterize(s_vert_sp_idris, r_template, field = s_vert_sp_idris$solution_1)
	prior_fly_feat_r_idris <- rasterize(s_fly_sp_idris, r_template, field = s_fly_sp_idris$solution_1)
	prior_plant_feat_r_idris <- rasterize(s_plant_sp_idris, r_template, field = s_plant_sp_idris$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_idris_u <- st_union(s_cond1_idris_sf)
	s_acd1_idris_u <- st_union(s_acd1_idris_sf)
	s_corr1_idris_u <- st_union(s_corr1_idris_sf)
	
	# ----------------------
	#  Convert to sp objects
	s_cond_sp_idris <- as(s_cond1_idris_u, 'Spatial')
	s_acd_sp_idris <- as(s_acd1_idris_u, 'Spatial')
	s_corr_sp_idris <- as(s_corr1_idris_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_feat_r_idris <- raster::mask(elev_conn_feat_r, s_cond_sp_idris)
	prior_acd_feat_r_idris <- raster::mask(acd_feat_r, s_acd_sp_idris)
	prior_corr_feat_r_idris <- raster::mask(corr_feat_r, s_corr_sp_idris)
	


# =============================================================================
#  STEP 4: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together in Idris area.
# =============================================================================	

	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in_idris <- stack(prior_vert_feat_r_idris, prior_fly_feat_r_idris, prior_plant_feat_r_idris, 
		prior_cond_feat_r_idris, prior_corr_feat_r_idris, prior_acd_feat_r_idris)

	# ----------------------
	#  Set up problem
	equal_targets <- 0.65
	p1_idris <- problem(x = pu_idris_in, features = all_feat_in_idris, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(equal_targets) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1_idris <- solve(p1_idris)
	s1_idris_sf <- st_as_sf(s1_idris) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_idris_area_h <- as.integer(sum(st_area(s1_idris_sf)) * 0.0001)
	s1_idris_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_idris1 <- as(s1_idris_sf, 'Spatial')

	# ----------------------
	#  Set up problem for round 3.
	p2_idris <- problem(x = pu_idris_in, features = all_feat_in_idris, cost_column = "area_h") %>%
		add_max_utility_objective(100000) %>%
		add_locked_in_constraints(locked_in_idris1) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s2_idris <- solve(p2_idris)
	s2_idris_sf <- st_as_sf(s2_idris) %>%
		dplyr::filter(solution_1 == 1)
	s2_idris_area_h <- as.integer(sum(st_area(s2_idris_sf)) * 0.0001)
	s2_idris_area_h
	

	
# =============================================================================
#  STEP 5: Combine outputs from SFI and Idris prioritizations and Deramakot as new locked in area 
#   for prioritization across Sabah 
# =============================================================================	
	
	#  Combine prioritized areas of SFI and Idris
	locked_in_sfi_idris <-  rbind(s2_sfi_sf, s2_idris_sf)
	locked_in_all_sf <- rbind(locked_in_sfi_idris, locked_in_tmp)
	locked_in_all <- as(locked_in_all_sf, 'Spatial')
	
	
	
# =============================================================================
#  STEP 6: Set up and solve problem for species ranges for whole of Sabah using already 
#   prioritized area in SFI/Idris, locked in Dermakot, and locked out space between 
#   Crocker Range - Kinabalu Park
# =============================================================================

	# ----------------------
	#  Set up problem 2 for vertebrates
	p_vert2 <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_vert2 <- solve(p_vert2)
	s_vert2_sf <- st_as_sf(s_vert2) %>%
		dplyr::filter(solution_1 == 1)
	s_vert2_area_h <- as.integer(sum(st_area(s_vert2_sf)) * 0.0001)
	s_vert2_area_h
	
	
	# ----------------------
	#  Set up problem 2 for butterflies
	p_fly2 <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>% # try increase to 0.25
		add_feature_weights(fly_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_fly2 <- solve(p_fly2)
	s_fly2_sf <- st_as_sf(s_fly2) %>%
		dplyr::filter(solution_1 == 1)
	s_fly2_area_h <- as.integer(sum(st_area(s_fly2_sf)) * 0.0001)
	s_fly2_area_h
	
	
	# ----------------------
	#  Set up problem 2 for plants
	p_plant2 <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_plant2 <- solve(p_plant2)
	s_plant2_sf <- st_as_sf(s_plant2) %>%
		dplyr::filter(solution_1 == 1)
	s_plant2_area_h <- as.integer(sum(st_area(s_plant2_sf)) * 0.0001)
	s_plant2_area_h
	
	
	# ----------------------
	#  Set up problem 2 for Condatis connectivity
	p_cond2 <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.275) %>%
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_cond2 <- solve(p_cond2)
	s_cond2_sf <- st_as_sf(s_cond2) %>%
		dplyr::filter(solution_1 == 1)
	s_cond2_area_h <- as.integer(sum(st_area(s_cond2_sf)) * 0.0001)
	s_cond2_area_h
	
	
	# ----------------------
	#  Set up problem 2 for ACD
	p_acd2 <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.15) %>%
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_acd2 <- solve(p_acd2)
	s_acd2_sf <- st_as_sf(s_acd2) %>%
		dplyr::filter(solution_1 == 1)
	s_acd2_area_h <- as.integer(sum(st_area(s_acd2_sf)) * 0.0001)
	s_acd2_area_h
	
	
	# ----------------------
	#  Set up problem 2 for corridors
	p_corr2 <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.05) %>%
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_corr2 <- solve(p_corr2)
	s_corr2_sf <- st_as_sf(s_corr2) %>%
		dplyr::filter(solution_1 == 1)
	s_corr2_area_h <- as.integer(sum(st_area(s_corr2_sf)) * 0.0001)
	s_corr2_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_sp <- as(s_vert2_sf, 'Spatial')
	s_fly_sp <- as(s_fly2_sf, 'Spatial')
	s_plant_sp <- as(s_plant2_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond2_u <- st_union(s_cond2_sf)
	s_acd2_u <- st_union(s_acd2_sf)
	s_corr2_u <- st_union(s_corr2_sf)

	# ----------------------
	#  Convert to sp objects
	s_cond_sp <- as(s_cond2_u, 'Spatial')
	s_acd_sp <- as(s_acd2_u, 'Spatial')
	s_corr_sp <- as(s_corr2_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_feat_r <- raster::mask(elev_conn_feat_r, s_cond_sp)
	prior_acd_feat_r <- raster::mask(acd_feat_r, s_acd_sp)
	prior_corr_feat_r <- raster::mask(corr_feat_r, s_corr_sp)
	
	
	
# =============================================================================
#  STEP 7: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, prior_cond_feat_r, 
		prior_corr_feat_r, prior_acd_feat_r)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.85
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	#  Solve round 2.
	s2 <- solve(p2)
	s2_sf <- st_as_sf(s2) %>%
		dplyr::filter(solution_1 == 1)
	s2_area_h <- as.integer(sum(st_area(s2_sf)) * 0.0001)
	s2_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new <- as(s2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(410000) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	s3_sf <- st_as_sf(s3) %>%
		dplyr::filter(solution_1 == 1)
	s3_area_h <- as.integer(sum(st_area(s3_sf)) * 0.0001)
	s3_area_h
	
	
# =============================================================================
#  STEP 8: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	
	
	# ----------------------
	#  Save output.
	scen4_tmp <- as(s3_sf, "Spatial")
	scen4_sf <- st_as_sf(scen4_tmp) %>%	
		mutate(ras_val = 1)
	scen4_out <- s3
	save(scen4_sf, file = "scenario_outputs/scen4/scen4_sf.Rdata")
	save(scen4_out, file = "scenario_outputs/scen4/scen4_out.Rdata")
	
	
	
	
	
	
	
	
	
	
	
# =============================================================================
#  Run prioritzation with BLM constraint
# =============================================================================	
	


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges within SFI area
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1_blm_sfi <- problem(x = pu_sfi_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert1_blm_sfi <- solve(p_vert1_blm_sfi)
	s_vert1_blm_sfi_sf <- st_as_sf(s_vert1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_vert1_blm_sfi_area_h <- as.integer(sum(st_area(s_vert1_blm_sfi_sf)) * 0.0001)
	s_vert1_blm_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1_blm_sfi <- problem(x = pu_sfi_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1_blm_sfi <- solve(p_fly1_blm_sfi)
	s_fly1_blm_sfi_sf <- st_as_sf(s_fly1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_blm_sfi_area_h <- as.integer(sum(st_area(s_fly1_blm_sfi_sf)) * 0.0001)
	s_fly1_blm_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for plants
	p_plant1_blm_sfi <- problem(x = pu_sfi_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1_blm_sfi <- solve(p_plant1_blm_sfi)
	s_plant1_blm_sfi_sf <- st_as_sf(s_plant1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_blm_sfi_area_h <- as.integer(sum(st_area(s_plant1_blm_sfi_sf)) * 0.0001)
	s_plant1_blm_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1_blm_sfi <- problem(x = pu_sfi_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.3) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1_blm_sfi <- solve(p_cond1_blm_sfi)
	s_cond1_blm_sfi_sf <- st_as_sf(s_cond1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_blm_sfi_area_h <- as.integer(sum(st_area(s_cond1_blm_sfi_sf)) * 0.0001)
	s_cond1_blm_sfi_area_h
	

	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1_blm_sfi <- problem(x = pu_sfi_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.25) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1_blm_sfi <- solve(p_acd1_blm_sfi)
	s_acd1_blm_sfi_sf <- st_as_sf(s_acd1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_blm_sfi_area_h <- as.integer(sum(st_area(s_acd1_blm_sfi_sf)) * 0.0001)
	s_acd1_blm_sfi_area_h
	
	
	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1_blm_sfi <- problem(x = pu_sfi_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(0.45) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1_blm_sfi <- solve(p_corr1_blm_sfi)
	s_corr1_blm_sfi_sf <- st_as_sf(s_corr1_blm_sfi) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_blm_sfi_area_h <- as.integer(sum(st_area(s_corr1_blm_sfi_sf)) * 0.0001)
	s_corr1_blm_sfi_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_blm_sp_sfi <- as(s_vert1_blm_sfi_sf, 'Spatial')
	s_fly_blm_sp_sfi <- as(s_fly1_blm_sfi_sf, 'Spatial')
	s_plant_blm_sp_sfi <- as(s_plant1_blm_sfi_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_blm_feat_r_sfi <- rasterize(s_vert_blm_sp_sfi, r_template, field = s_vert_blm_sp_sfi$solution_1)
	prior_fly_blm_feat_r_sfi <- rasterize(s_fly_blm_sp_sfi, r_template, field = s_fly_blm_sp_sfi$solution_1)
	prior_plant_blm_feat_r_sfi <- rasterize(s_plant_blm_sp_sfi, r_template, field = s_plant_blm_sp_sfi$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_blm_sfi_u <- st_union(s_cond1_blm_sfi_sf)
	s_acd1_blm_sfi_u <- st_union(s_acd1_blm_sfi_sf)
	s_corr1_blm_sfi_u <- st_union(s_corr1_blm_sfi_sf)
		
	# s_corr1_u_tmp1 <- st_cast(st_union(s_corr1_sf), 'POLYGON')
		# s_corr1_u_tmp2 <- st_geometry(s_corr1_u_tmp1)
		# s_corr1_u_tmp3 <- as(s_corr1_u_tmp2, 'Spatial')
		# s_corr1_u_tmp4 <- st_as_sf(s_corr1_u_tmp3)
		# s_corr1_u_tmp5 <- s_corr1_u_tmp4[-8,]
	# s_corr1_u <- st_cast(s_corr1_u_tmp5, 'MULTIPOLYGON')
	
	# ----------------------
	#  Convert to sp objects
	s_cond_blm_sp_sfi <- as(s_cond1_blm_sfi_u, 'Spatial')
	s_acd_blm_sp_sfi <- as(s_acd1_blm_sfi_u, 'Spatial')
	s_corr_blm_sp_sfi <- as(s_corr1_blm_sfi_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_blm_feat_r_sfi <- raster::mask(elev_conn_feat_r, s_cond_blm_sp_sfi)
	prior_acd_blm_feat_r_sfi <- raster::mask(acd_feat_r, s_acd_blm_sp_sfi)
	prior_corr_blm_feat_r_sfi <- raster::mask(corr_feat_r, s_corr_blm_sp_sfi)
	
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together in SFI area.
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 1.
	all_feat_blm_in_sfi <- stack(prior_vert_blm_feat_r_sfi, prior_fly_blm_feat_r_sfi, prior_plant_blm_feat_r_sfi, 
		prior_cond_blm_feat_r_sfi, prior_corr_blm_feat_r_sfi, prior_acd_blm_feat_r_sfi)

	# ----------------------
	#  Set up problem
	equal_targets <- 0.45
	p1_sfi_blm <- problem(x = pu_sfi_in, features = all_feat_blm_in_sfi, cost_column = "area_h") %>%
		add_max_features_objective(50000) %>%
		add_relative_targets(equal_targets) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1_sfi_blm <- solve(p1_sfi_blm)
	s1_sfi_blm_sf <- st_as_sf(s1_sfi_blm) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_sfi_blm_area_h <- as.integer(sum(st_area(s1_sfi_blm_sf)) * 0.0001)
	s1_sfi_blm_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_sfi1_blm <- as(s1_sfi_blm_sf, 'Spatial')

	# ----------------------
	#  Set up problem for round 3.
	p2_sfi_blm <- problem(x = pu_sfi_in, features = all_feat_blm_in_sfi, cost_column = "area_h") %>%
		add_max_utility_objective(50000) %>%
		add_locked_in_constraints(locked_in_sfi1_blm) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s2_sfi_blm <- solve(p2_sfi_blm)
	s2_sfi_blm_sf <- st_as_sf(s2_sfi_blm) %>%
		dplyr::filter(solution_1 == 1)
	s2_sfi_blm_area_h <- as.integer(sum(st_area(s2_sfi_blm_sf)) * 0.0001)
	s2_sfi_blm_area_h
	
	
	
# =============================================================================
#  STEP 3: Set up and solve problem for species ranges within Idris area
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1_blm_idris <- problem(x = pu_idris_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert1_blm_idris <- solve(p_vert1_blm_idris)
	s_vert1_blm_idris_sf <- st_as_sf(s_vert1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_vert1_blm_idris_area_h <- as.integer(sum(st_area(s_vert1_blm_idris_sf)) * 0.0001)
	s_vert1_blm_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1_blm_idris <- problem(x = pu_idris_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1_blm_idris <- solve(p_fly1_blm_idris)
	s_fly1_blm_idris_sf <- st_as_sf(s_fly1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_blm_idris_area_h <- as.integer(sum(st_area(s_fly1_blm_idris_sf)) * 0.0001)
	s_fly1_blm_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for plants
	p_plant1_blm_idris <- problem(x = pu_idris_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1_blm_idris <- solve(p_plant1_blm_idris)
	s_plant1_blm_idris_sf <- st_as_sf(s_plant1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_blm_idris_area_h <- as.integer(sum(st_area(s_plant1_blm_idris_sf)) * 0.0001)
	s_plant1_blm_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1_blm_idris <- problem(x = pu_idris_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.6) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1_blm_idris <- solve(p_cond1_blm_idris)
	s_cond1_blm_idris_sf <- st_as_sf(s_cond1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_blm_idris_area_h <- as.integer(sum(st_area(s_cond1_blm_idris_sf)) * 0.0001)
	s_cond1_blm_idris_area_h
	

	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1_blm_idris <- problem(x = pu_idris_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.25) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1_blm_idris <- solve(p_acd1_blm_idris)
	s_acd1_blm_idris_sf <- st_as_sf(s_acd1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_blm_idris_area_h <- as.integer(sum(st_area(s_acd1_blm_idris_sf)) * 0.0001)
	s_acd1_blm_idris_area_h
	
	
	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1_blm_idris <- problem(x = pu_idris_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(0.6) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1_blm_idris <- solve(p_corr1_blm_idris)
	s_corr1_blm_idris_sf <- st_as_sf(s_corr1_blm_idris) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_blm_idris_area_h <- as.integer(sum(st_area(s_corr1_blm_idris_sf)) * 0.0001)
	s_corr1_blm_idris_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_blm_sp_idris <- as(s_vert1_blm_idris_sf, 'Spatial')
	s_fly_blm_sp_idris <- as(s_fly1_blm_idris_sf, 'Spatial')
	s_plant_blm_sp_idris <- as(s_plant1_blm_idris_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_blm_feat_r_idris <- rasterize(s_vert_blm_sp_idris, r_template, field = s_vert_blm_sp_idris$solution_1)
	prior_fly_blm_feat_r_idris <- rasterize(s_fly_blm_sp_idris, r_template, field = s_fly_blm_sp_idris$solution_1)
	prior_plant_blm_feat_r_idris <- rasterize(s_plant_blm_sp_idris, r_template, field = s_plant_blm_sp_idris$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_blm_idris_u <- st_union(s_cond1_blm_idris_sf)
	s_acd1_blm_idris_u <- st_union(s_acd1_blm_idris_sf)
	s_corr1_blm_idris_u <- st_union(s_corr1_blm_idris_sf)
	
	# ----------------------
	#  Convert to sp objects
	s_cond_blm_sp_idris <- as(s_cond1_blm_idris_u, 'Spatial')
	s_acd_blm_sp_idris <- as(s_acd1_blm_idris_u, 'Spatial')
	s_corr_blm_sp_idris <- as(s_corr1_blm_idris_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_blm_feat_r_idris <- raster::mask(elev_conn_feat_r, s_cond_blm_sp_idris)
	prior_acd_blm_feat_r_idris <- raster::mask(acd_feat_r, s_acd_blm_sp_idris)
	prior_corr_blm_feat_r_idris <- raster::mask(corr_feat_r, s_corr_blm_sp_idris)
	
	
	
# =============================================================================
#  STEP 4: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together in Idris area.
# =============================================================================	

	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_blm_in_idris <- stack(prior_vert_blm_feat_r_idris, prior_fly_blm_feat_r_idris, prior_plant_blm_feat_r_idris, 
		prior_cond_blm_feat_r_idris, prior_corr_blm_feat_r_idris, prior_acd_blm_feat_r_idris)

	# ----------------------
	#  Set up problem
	equal_targets <- 0.65
	p1_idris_blm <- problem(x = pu_idris_in, features = all_feat_blm_in_idris, cost_column = "area_h") %>%
		add_max_features_objective(100000) %>%
		add_relative_targets(equal_targets) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1_idris_blm <- solve(p1_idris_blm)
	s1_idris_blm_sf <- st_as_sf(s1_idris_blm) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_idris_blm_area_h <- as.integer(sum(st_area(s1_idris_blm_sf)) * 0.0001)
	s1_idris_blm_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_idris1_blm <- as(s1_idris_blm_sf, 'Spatial')

	# ----------------------
	#  Set up problem for round 3.
	p2_idris_blm <- problem(x = pu_idris_in, features = all_feat_blm_in_idris, cost_column = "area_h") %>%
		add_max_utility_objective(100000) %>%
		add_locked_in_constraints(locked_in_idris1_blm) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s2_idris_blm <- solve(p2_idris_blm)
	s2_idris_blm_sf <- st_as_sf(s2_idris_blm) %>%
		dplyr::filter(solution_1 == 1)
	s2_idris_blm_area_h <- as.integer(sum(st_area(s2_idris_blm_sf)) * 0.0001)
	s2_idris_blm_area_h
	

	
# =============================================================================
#  STEP 5: Combine outputs from SFI and Idris prioritizations and Deramakot as new locked in area 
#   for prioritization across Sabah 
# =============================================================================	
	
	#  Combine prioritized areas of SFI and Idris
	locked_in_sfi_idris_blm <-  rbind(s2_sfi_blm_sf, s2_idris_blm_sf)
	locked_in_all_blm_sf <- rbind(locked_in_sfi_idris_blm, locked_in_tmp)
	locked_in_all_blm <- as(locked_in_all_blm_sf, 'Spatial')
	
	
	
# =============================================================================
#  STEP 6: Set up and solve problem for species ranges for whole of Sabah using already 
#   prioritized area in SFI/Idris, locked in Dermakot, and locked out space between 
#   Crocker Range - Kinabalu Park
# =============================================================================

	# ----------------------
	#  Set up problem 2 for vertebrates
	p_vert2_blm <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_vert2_blm <- solve(p_vert2_blm)
	s_vert2_blm_sf <- st_as_sf(s_vert2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_vert2_blm_area_h <- as.integer(sum(st_area(s_vert2_blm_sf)) * 0.0001)
	s_vert2_blm_area_h
	
	
	# ----------------------
	#  Set up problem 2 for butterflies
	p_fly2_blm <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>% # try increase to 0.25
		add_feature_weights(fly_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_fly2_blm <- solve(p_fly2_blm)
	s_fly2_blm_sf <- st_as_sf(s_fly2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_fly2_blm_area_h <- as.integer(sum(st_area(s_fly2_blm_sf)) * 0.0001)
	s_fly2_blm_area_h
	
	
	# ----------------------
	#  Set up problem 2 for plants
	p_plant2_blm <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_plant2_blm <- solve(p_plant2_blm)
	s_plant2_blm_sf <- st_as_sf(s_plant2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_plant2_blm_area_h <- as.integer(sum(st_area(s_plant2_blm_sf)) * 0.0001)
	s_plant2_blm_area_h
	
	
	# ----------------------
	#  Set up problem 2 for Condatis connectivity
	p_cond2_blm <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.275) %>%
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_cond2_blm <- solve(p_cond2_blm)
	s_cond2_blm_sf <- st_as_sf(s_cond2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_cond2_blm_area_h <- as.integer(sum(st_area(s_cond2_blm_sf)) * 0.0001)
	s_cond2_blm_area_h
	
	
	# ----------------------
	#  Set up problem 2 for ACD
	p_acd2_blm <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.2) %>%
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_acd2_blm <- solve(p_acd2_blm)
	s_acd2_blm_sf <- st_as_sf(s_acd2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_acd2_blm_area_h <- as.integer(sum(st_area(s_acd2_blm_sf)) * 0.0001)
	s_acd2_blm_area_h
	
	
	# ----------------------
	#  Set up problem 2 for corridors
	p_corr2_blm <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(0.65) %>%
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 2.
	s_corr2_blm <- solve(p_corr2_blm)
	s_corr2_blm_sf <- st_as_sf(s_corr2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s_corr2_blm_area_h <- as.integer(sum(st_area(s_corr2_blm_sf)) * 0.0001)
	s_corr2_blm_area_h
	
	
	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_blm_sp <- as(s_vert2_blm_sf, 'Spatial')
	s_fly_blm_sp <- as(s_fly2_blm_sf, 'Spatial')
	s_plant_blm_sp <- as(s_plant2_blm_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_blm_feat_r <- rasterize(s_vert_blm_sp, r_template, field = s_vert_blm_sp$solution_1)
	prior_fly_blm_feat_r <- rasterize(s_fly_blm_sp, r_template, field = s_fly_blm_sp$solution_1)
	prior_plant_blm_feat_r <- rasterize(s_plant_blm_sp, r_template, field = s_plant_blm_sp$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond2_blm_u <- st_union(s_cond2_blm_sf)
	s_acd2_blm_u <- st_union(s_acd2_blm_sf)
	s_corr2_blm_u <- st_union(s_corr2_blm_sf)

	# ----------------------
	#  Convert to sp objects
	s_cond_blm_sp <- as(s_cond2_blm_u, 'Spatial')
	s_acd_blm_sp <- as(s_acd2_blm_u, 'Spatial')
	s_corr_blm_sp <- as(s_corr2_blm_u, 'Spatial')

	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_blm_feat_r <- raster::mask(elev_conn_feat_r, s_cond_blm_sp)
	prior_acd_blm_feat_r <- raster::mask(acd_feat_r, s_acd_blm_sp)
	prior_corr_blm_feat_r <- raster::mask(corr_feat_r, s_corr_blm_sp)
	
	
	
# =============================================================================
#  STEP 7: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in_blm <- stack(prior_vert_blm_feat_r, prior_fly_blm_feat_r, prior_plant_blm_feat_r,
		prior_cond_blm_feat_r, prior_corr_blm_feat_r, prior_acd_blm_feat_r)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.8
	p2_blm <- problem(x = pu_in, features = all_feat_in_blm, cost_column = "area_h") %>%
		add_max_features_objective(410000) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in_all_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	#  Solve round 2.
	s2_blm <- solve(p2_blm)
	s2_blm_sf <- st_as_sf(s2_blm) %>%
		dplyr::filter(solution_1 == 1)
	s2_blm_area_h <- as.integer(sum(st_area(s2_blm_sf)) * 0.0001)
	s2_blm_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_blm <- as(s2_blm_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem for round 3.
	p3_blm <- problem(x = pu_in, features = all_feat_in_blm, cost_column = "area_h") %>%
		add_max_utility_objective(410000) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s3_blm <- solve(p3_blm)
	s3_blm_sf <- st_as_sf(s3_blm) %>%
		dplyr::filter(solution_1 == 1)
	s3_blm_area_h <- as.integer(sum(st_area(s3_blm_sf)) * 0.0001)
	s3_blm_area_h
	
	

# =============================================================================
#  STEP 8: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	
	
	# ----------------------
	#  Save output.
	scen4_blm_tmp <- as(s3_blm_sf, "Spatial")
	scen4_blm_sf <- st_as_sf(scen4_blm_tmp) %>%	
		mutate(ras_val = 1)
	scen4_blm_out <- s3_blm
	
	save(scen4_blm_out, file = "scenario_outputs/scen4/scen4_blm_out.Rdata")
	save(scen4_blm_sf, file = "scenario_outputs/scen4/scen4_blm_sf.Rdata")
	
