###############################################################################
#  Run optimization with 'Prioritizr' package for taxa
#  March 1, 2018; updated May 5, 2018.
#  Script to generate optimize selected land units for protection of a variety of input 
#   conservation features.
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
# install.packages("C:/gurobi751/win64/R/gurobi_7.5-1.zip", repos = NULL)
library(gurobi)
#  devtools::install_github("prioritizr/priortizr") # can use official CRAN version
library(prioritizr)
library(ggplot2)



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Boundaries
	load(file = "C:/Users/saraw/Desktop/5_5_18/sabah_boder.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/sabah_tpa.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/acd_agg_sf_for.Rdata")
	
	# ----------------------
	#  Load species ranges (1 layer per species).
	vert_feat_in_single <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_all.grd")
	fly_feat_in_single <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_all.grd")	
	plant_feat_in_single <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plant_all.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_rep_weight.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_rep_weight.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plant_rep_weight.Rdata")
	
	# ----------------------
	#  Set up problem for connectivity and carbon input layers
	elev_conn_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/elev_conn_feat_in.grd")
	corr_feat_r <- raster("C:/Users/saraw/Desktop/tmp/corr_out_r_small.grd")
	acd_feat_r <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/acd_feat_in.grd")
	
	
	
# =============================================================================
#  Set up planning units WITHOUT existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/sfi_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/idris_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/deram_grid.Rdata")

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
	
	sfi_idris_grid <- rbind(sfi_grid, idris_grid)
	pu_sf_tmp <- sfi_idris_grid %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_sf_tmp$area_h_tmp <- as.numeric(pu_sf_tmp$area_h_tmp)
	pu_sf <- pu_sf_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100)
	pu_sfi_idris_in <- as(pu_sf, "Spatial")
	
	locked_out <- as(deram_grid, 'Spatial')
	
	
	
# # =============================================================================
# #  Set up planning units WITH existing TPAs
# # =============================================================================
	
	# # ----------------------
	# #  Load planning units with existing TPAs
	# load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid_w_pa.Rdata")
	# load(file = "C:/Users/saraw/Desktop/5_5_18/tpa_grid.Rdata")
	# load(file = "C:/Users/saraw/Desktop/5_5_18/sfi_grid_w_pa.Rdata")
	# load(file = "C:/Users/saraw/Desktop/5_5_18/idris_grid_w_pa.Rdata")
	# load(file = "C:/Users/saraw/Desktop/5_5_18/deram_grid_w_pa.Rdata")
	
	# const_cost_h <- 500
	# pu_sf_tmp <- sa_grid_w_pa %>%
		# mutate(area_tmp = st_area(.) * 0.0001) %>%
		# separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	# pu_sf_tmp$area_h_tmp <- as.numeric(pu_sf_tmp$area_h_tmp)
	# pu_sf <- pu_sf_tmp %>%
		# mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		# mutate(const_cost = const_cost_h) %>%
		# dplyr::select(-area_h_tmp) %>%
		# dplyr::filter(area_h > 100)
	# pu_in <- as(pu_sf, "Spatial")
	
	# locked_in <- (tpa_grid, 'Spatial')	
	# tpa_area_h <- as.integer(sum(st_area(tpa_grid)) * 0.0001)
	
	# sfi_idris_grid_w_pa <- rbind(sfi_grid_w_pa, idris_grid_w_pa)
	# pu_sf_tmp <- sfi_idris_grid_w_pa %>%
		# mutate(area_tmp = st_area(.) * 0.0001) %>%
		# separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	# pu_sf_tmp$area_h_tmp <- as.numeric(pu_sf_tmp$area_h_tmp)
	# pu_sf <- pu_sf_tmp %>%
		# mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		# mutate(const_cost = const_cost_h) %>%
		# dplyr::select(-area_h_tmp) %>%
		# dplyr::filter(area_h > 100)
	# pu_sfi_idris_in <- as(pu_sf, "Spatial")
	
	# locked_out <- as(deram_grid_w_pa, 'Spatial')
	
	
# =============================================================================
#  Set up raster template.
# =============================================================================
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	temp <- fly_feat_in_single[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	

# =============================================================================
#  Initiate problem and solve using single repetition of species ranges and weighting.
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1 <- problem(x = pu_sfi_idris_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert1 <- solve(p_vert1)
	s_vert1_sf <- st_as_sf(s_vert1) %>%
		dplyr::filter(solution_1 == 1)
	s_vert1_area_h <- as.integer(sum(st_area(s_vert1_sf)) * 0.0001)
	s_vert1_area_h
	

	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1 <- problem(x = pu_sfi_idris_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1 <- solve(p_fly1)
	s_fly1_sf <- st_as_sf(s_fly1) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_area_h <- as.integer(sum(st_area(s_fly1_sf)) * 0.0001)
	s_fly1_area_h
	
	
	
	# ----------------------
	#  Set up problem 1 for plants
	p_plant1 <- problem(x = pu_sfi_idris_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1 <- solve(p_plant1)
	s_plant1_sf <- st_as_sf(s_plant1) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_area_h <- as.integer(sum(st_area(s_plant1_sf)) * 0.0001)
	s_plant1_area_h
	
	

	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1 <- problem(x = pu_sfi_idris_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.6) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1 <- solve(p_cond1)
	s_cond1_sf <- st_as_sf(s_cond1) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_area_h <- as.integer(sum(st_area(s_cond1_sf)) * 0.0001)
	s_cond1_area_h
	


	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1 <- problem(x = pu_sfi_idris_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.45) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1 <- solve(p_acd1)
	s_acd1_sf <- st_as_sf(s_acd1) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_area_h <- as.integer(sum(st_area(s_acd1_sf)) * 0.0001)
	s_acd1_area_h
	
	

	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1 <- problem(x = pu_sfi_idris_in, features = corr_feat_r, cost_column = "area_h") %>%
		#add_max_utility_objective(200000) %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.95) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1 <- solve(p_corr1)
	s_corr1_sf <- st_as_sf(s_corr1) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_area_h <- as.integer(sum(st_area(s_corr1_sf)) * 0.0001)
	s_corr1_area_h
	
	
	
	
	s_vert_sp_si <- as(s_vert1_sf, 'Spatial')
	s_fly_sp_si <- as(s_fly1_sf, 'Spatial')
	s_plant_sp_si <- as(s_plant1_sf, 'Spatial')

	prior_vert_feat_r_si <- rasterize(s_vert_sp_si, r_template, field = s_vert_sp_si$solution_1)
	prior_fly_feat_r_si <- rasterize(s_fly_sp_si, r_template, field = s_fly_sp_si$solution_1)
	prior_plant_feat_r_si <- rasterize(s_plant_sp_si, r_template, field = s_plant_sp_si$solution_1)
	
	
	
	s_cond1_u <- st_union(s_cond1_sf)
	s_acd1_u <- st_union(s_acd1_sf)
	s_corr1_u_tmp1 <- st_cast(st_union(s_corr1_sf), 'POLYGON')
		s_corr1_u_tmp2 <- st_geometry(s_corr1_u_tmp1)
		s_corr1_u_tmp3 <- as(s_corr1_u_tmp2, 'Spatial')
		s_corr1_u_tmp4 <- st_as_sf(s_corr1_u_tmp3)
		s_corr1_u_tmp5 <- s_corr1_u_tmp4[-8,]
	s_corr1_u <- st_cast(s_corr1_u_tmp5, 'MULTIPOLYGON')
	
	s_cond_sp_si <- as(s_cond1_u, 'Spatial')
	s_acd_sp_si <- as(s_acd1_u, 'Spatial')
	s_corr_sp_si <- as(s_corr1_u, 'Spatial')
	
	prior_cond_feat_r_si <- raster::mask(elev_conn_feat_r, s_cond_sp_si)
	prior_acd_feat_r_si <- raster::mask(acd_feat_r, s_acd_sp_si)
	prior_corr_feat_r_si <- raster::mask(corr_feat_r, s_corr_sp_si)
	
	

	
	
	
	
	
	
	
	# ----------------------
	#  Set up problem 2 for vertebrates
	p_vert2 <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
		add_max_features_objective(350000) %>%
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
		add_max_features_objective(350000) %>%
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
		add_max_features_objective(350000) %>%
		add_relative_targets(0.2) %>%
		add_locked_in_constraints(locked_in_all) %>%
		#add_locked_out_constraints(locked_out) %>%
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
		add_max_features_objective(350000) %>%
		add_relative_targets(0.18) %>%
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
		add_max_features_objective(350000) %>%
		add_relative_targets(0.65) %>%
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
	
	
	
	s_vert_sp <- as(s_vert2_sf, 'Spatial')
	s_fly_sp <- as(s_fly2_sf, 'Spatial')
	s_plant_sp <- as(s_plant2_sf, 'Spatial')

	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	
	
	s_cond2_u <- st_union(s_cond2_sf)
	s_acd2_u <- st_union(s_acd2_sf)
	s_corr2_u <- st_union(s_corr2_sf)

	s_cond_sp <- as(s_cond2_u, 'Spatial')
	s_acd_sp <- as(s_acd2_u, 'Spatial')
	s_corr_sp <- as(s_corr2_u, 'Spatial')
	
	prior_cond_feat_r <- raster::mask(elev_conn_feat_r, s_cond_sp)
	prior_acd_feat_r <- raster::mask(acd_feat_r, s_acd_sp)
	prior_corr_feat_r <- raster::mask(corr_feat_r, s_corr_sp)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_vert2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for vertebrates
	p_vert3 <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_vert3 <- solve(p_vert3)
	s_vert3_sf <- st_as_sf(s_vert3) %>%
		dplyr::filter(solution_1 == 1)
	s_vert3_area_h <- as.integer(sum(st_area(s_vert3_sf)) * 0.0001)
	s_vert3_area_h
	s_vert_sf_upd <- s_vert3_sf
	
	
	
	
	# ----------------------
	#  Set up problem 1 for butterflies
	p_fly1 <- problem(x = pu_sfi_idris_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_fly1 <- solve(p_fly1)
	s_fly1_sf <- st_as_sf(s_fly1) %>%
		dplyr::filter(solution_1 == 1)
	s_fly1_area_h <- as.integer(sum(st_area(s_fly1_sf)) * 0.0001)
	s_fly1_area_h
	
		# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_fly1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 2 for butterflies
	p_fly2 <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_locked_in_constraints(locked_in) %>%
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
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_fly2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for butterflies
	p_fly3 <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_fly3 <- solve(p_fly3)
	s_fly3_sf <- st_as_sf(s_fly3) %>%
		dplyr::filter(solution_1 == 1)
	s_fly3_area_h <- as.integer(sum(st_area(s_fly3_sf)) * 0.0001)
	s_fly3_area_h
	s_fly_sf_upd <- s_fly3_sf
	

	# ----------------------
	#  Set up problem 1 for plants
	p_plant1 <- problem(x = pu_sfi_idris_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.35) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant1 <- solve(p_plant1)
	s_plant1_sf <- st_as_sf(s_plant1) %>%
		dplyr::filter(solution_1 == 1)
	s_plant1_area_h <- as.integer(sum(st_area(s_plant1_sf)) * 0.0001)
	s_plant1_area_h
	
	r_plant1 <- feat_rep(p_plant1, s_plant1[1, "solution_1"])
	r_plant1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_plant1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 2 for butterflies
	p_plant2 <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.2) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_locked_in_constraints(locked_in) %>%
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
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_plant2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for butterflies
	p_plant3 <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_feature_weights(plant_rep_weight) %>% 
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_plant3 <- solve(p_plant3)
	s_plant3_sf <- st_as_sf(s_plant3) %>%
		dplyr::filter(solution_1 == 1)
	s_plant3_area_h <- as.integer(sum(st_area(s_plant3_sf)) * 0.0001)
	s_plant3_area_h
	s_plant_sf_upd <- s_plant3_sf
	
	
	
	# ----------------------
	#  Save sf outputs
	save(s_vert_sf_upd, file = "C:/Users/saraw/Desktop/5_5_18/s_vert_sf_upd.Rdata")
	save(s_fly_sf_upd, file = "C:/Users/saraw/Desktop/5_5_18/s_fly_sf_upd.Rdata")
	save(s_plant_sf_upd, file = "C:/Users/saraw/Desktop/5_5_18/s_plant_sf_upd.Rdata")
	
	

# =============================================================================
#  Convert to rasters to serve as feature input layers in next round of prioritization.
# =============================================================================
	
	# ----------------------
	#  Convert outputs sp objects
	s_vert_sp_si <- as(s_vert1_sf, 'Spatial')
	s_fly_sp_si <- as(s_fly1_sf, 'Spatial')
	s_plant_sp_si <- as(s_plant1_sf, 'Spatial')
	
	s_vert_sp <- as(s_vert_sf_upd, 'Spatial')
	s_fly_sp <- as(s_fly_sf_upd, 'Spatial')
	s_plant_sp <- as(s_plant_sf_upd, 'Spatial')
	
	# ----------------------
	#  Rasterize
	prior_vert_feat_r_si <- rasterize(s_vert_sp_si, r_template, field = s_vert_sp_si$solution_1)
	prior_fly_feat_r_si <- rasterize(s_fly_sp_si, r_template, field = s_fly_sp_si$solution_1)
	prior_plant_feat_r_si <- rasterize(s_plant_sp_si, r_template, field = s_plant_sp_si$solution_1)
	
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	# ----------------------
	#  Save rasters
	writeRaster(prior_vert_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r_upd.grd")
	writeRaster(prior_fly_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r_upd.grd")
	writeRaster(prior_plant_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r_upd.grd")
	
	writeRaster(prior_vert_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r_upd_si.grd")
	writeRaster(prior_fly_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r_upd_si.grd")
	writeRaster(prior_plant_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r_upd_si.grd")
	
	
	
	
	
	
	# ----------------------
	#  Set up problem for Condatis elevational connectivity output
	elev_conn_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/elev_conn_feat_in.grd")
	
	# ----------------------
	#  Set up problem 1 for condatis
	p_cond1 <- problem(x = pu_sfi_idris_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.6) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_cond1 <- solve(p_cond1)
	s_cond1_sf <- st_as_sf(s_cond1) %>%
		dplyr::filter(solution_1 == 1)
	s_cond1_area_h <- as.integer(sum(st_area(s_cond1_sf)) * 0.0001)
	s_cond1_area_h
	
	r_cond1 <- feat_rep(p_cond1, s_cond1[1, "solution_1"])
	r_cond1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_cond1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 2 for Condatis
	p_cond2 <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.3) %>%
		add_locked_in_constraints(locked_in) %>%
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
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_cond2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for Condatis
	p_cond3 <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_cover_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_cond3 <- solve(p_cond3)
	s_cond3_sf <- st_as_sf(s_cond3) %>%
		dplyr::filter(solution_1 == 1)
	s_cond3_area_h <- as.integer(sum(st_area(s_cond3_sf)) * 0.0001)
	s_cond3_area_h
	
	
	
	# ----------------------
	#  Set up problem for ACD output
	acd <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/acd_feat_in.grd")
	
	# ----------------------
	#  Set up problem 1 for ACD
	p_acd1 <- problem(x = pu_sfi_idris_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.45) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_acd1 <- solve(p_acd1)
	s_acd1_sf <- st_as_sf(s_acd1) %>%
		dplyr::filter(solution_1 == 1)
	s_acd1_area_h <- as.integer(sum(st_area(s_acd1_sf)) * 0.0001)
	s_acd1_area_h
	
	r_acd1 <- feat_rep(p_acd1, s_acd1[1, "solution_1"])
	r_acd1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_acd1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 2 for ACD
	p_acd2 <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.4) %>%
		add_locked_in_constraints(locked_in) %>%
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
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_acd2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for ACD
	p_acd3 <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_cover_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_acd3 <- solve(p_acd3)
	s_acd3_sf <- st_as_sf(s_acd3) %>%
		dplyr::filter(solution_1 == 1)
	s_acd3_area_h <- as.integer(sum(st_area(s_acd3_sf)) * 0.0001)
	s_acd3_area_h
	
	
	
	
	corr_feat_r <- raster("C:/Users/saraw/Desktop/tmp/corr_out_r_small.grd")
	# ----------------------
	#  Set up problem 1 for corridor connectivity
	p_corr1 <- problem(x = pu_sfi_idris_in, features = corr_feat_r, cost_column = "area_h") %>%
		#add_max_utility_objective(200000) %>%
		add_max_features_objective(200000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.8) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1 <- solve(p_corr1)
	s_corr1_sf <- st_as_sf(s_corr1) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_area_h <- as.integer(sum(st_area(s_corr1_sf)) * 0.0001)
	s_corr1_area_h
	
	r_corr1 <- feat_rep(p_corr1, s_corr1[1, "solution_1"])
	r_corr1

	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_corr1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 2 for corridor connectivity
	p_corr2 <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.8) %>%
		add_locked_in_constraints(locked_in) %>%
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
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s_corr2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem 3 for corridor connectivity
	p_corr3 <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		add_max_cover_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s_corr3 <- solve(p_corr3)
	s_corr3_sf <- st_as_sf(s_corr3) %>%
		dplyr::filter(solution_1 == 1)
	s_corr3_area_h <- as.integer(sum(st_area(s_corr3_sf)) * 0.0001)
	s_corr3_area_h
	
	
	
	
	# ----------------------
	#  Convert outputs sp objects
	s_acd_sp_si <- as(s_acd1_sf, 'Spatial')
	s_cond_sp_si <- as(s_cond1_sf, 'Spatial')
	s_corr_sp_si <- as(s_corr1_sf, 'Spatial')
	
	s_acd_sp <- as(s_acd3_sf, 'Spatial')
	s_cond_sp <- as(s_cond3_sf, 'Spatial')
	s_corr_sp <- as(s_corr3_sf, 'Spatial')
	
	# ----------------------
	#  Rasterize
	prior_acd_feat_r_si <- rasterize(s_acd_sp_si, r_template, field = s_acd_sp_si$solution_1)
	prior_cond_feat_r_si <- rasterize(s_cond_sp_si, r_template, field = s_cond_sp_si$solution_1)
	prior_corr_feat_r_si <- rasterize(s_corr_sp_si, r_template, field = s_corr_sp_si$solution_1)
	
	prior_acd_feat_r <- rasterize(s_vert_sp, r_template, field = s_acd_sp$solution_1)
	prior_cond_feat_r <- rasterize(s_fly_sp, r_template, field = s_cond_sp$solution_1)
	prior_corr_feat_r <- rasterize(s_plant_sp, r_template, field = s_corr_sp$solution_1)
	
	# ----------------------
	#  Save rasters
	writeRaster(prior_vert_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r_upd.grd")
	writeRaster(prior_fly_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r_upd.grd")
	writeRaster(prior_plant_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r_upd.grd")
	
	writeRaster(prior_vert_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r_upd_si.grd")
	writeRaster(prior_fly_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r_upd_si.grd")
	writeRaster(prior_plant_feat_r_si, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r_upd_si.grd")
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# ----------------------
	#  Convert outputs sp objects
	s_cond_sp <- as(s_cond3_sf, 'Spatial')
	s_acd_sp <- as(s_acd3_sf, 'Spatial')
	s_corr_sp <- as(s_corr3_sf, 'Spatial')
	
	# ----------------------
	#  Rasterize
	prior_cond_feat_r <- rasterize(s_cond_sp, r_template, field = s_cond_sp$solution_1)
	prior_acd_feat_r <- rasterize(s_acd_sp, r_template, field = s_acd_sp$solution_1)
	prior_corr_feat_r <- rasterize(s_corr_sp, r_template, field = s_corr_sp$solution_1)
	
	# ----------------------
	#  Save rasters
	writeRaster(prior_vert_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r.grd")
	writeRaster(prior_fly_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r.grd")
	writeRaster(prior_plant_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r.grd")
	
	