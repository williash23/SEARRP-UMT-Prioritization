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
#  devtools::install_github("prioritizr/priortizr") #### Now unecessary as there is an official
#   CRAN version. 
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
	plant_feat_in_single <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plants_all.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_rep_weight.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_rep_weight.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plant_rep_weight.Rdata")
	
	
	
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
#  Initiate problem and solve first iteration using single repetition of species ranges and weighting.
# =============================================================================
	
	# ----------------------
	#  Set up problem for vertebrates
	p_vert <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(350000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_boundary_penalties(penalty = 0.00000001) %>%
		#add_locked_in_constraints(locked_in) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_vert <- solve(p_vert)
	s_vert_sf <- st_as_sf(s_vert) %>%
		dplyr::filter(solution_1 == 1)
	vert_area_h <- as.integer(sum(st_area(s_vert_sf)) * 0.0001)
	vert_area_h
	
	r_vert <- feat_rep(p_vert, s_vert[, "solution_1"])
	r_vert

	
	# ----------------------
	#  Set up problem for butterflies
	p_fly <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(350000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
		add_feature_weights(fly_rep_weight) %>% 
		add_boundary_penalties(penalty = 0.00000001) %>%
		#add_locked_in_constraints(locked_in) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve 
	s_fly <- solve(p_fly)
	s_fly_sf <- st_as_sf(s_fly) %>%
		dplyr::filter(solution_1 == 1)
	fly_area_h <- as.integer(sum(st_area(s_fly_sf)) * 0.0001)
	fly_area_h
	
	r_fly <- feat_rep(p_fly, s_fly[, "solution_1"])
	r_fly

	
	# ----------------------
	#  Set up problem for plants
	p_plant <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(350000 + tpa_area_h) %>%
		add_relative_targets(0.1) %>%
		add_feature_weights(vert_plant_weight) %>% 
		add_boundary_penalties(penalty = 0.00000001) %>%
		#add_locked_in_constraints(locked_in) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_plant <- solve(p_plant)
	s_plant_sf <- st_as_sf(s_plant) %>%
		dplyr::filter(solution_1 == 1)
	plant_area_h <- as.integer(sum(st_area(s_plant_sf)) * 0.0001)
	plant_area_h
	
	r_plant <- feat_rep(p_plant, s_plant[, "solution_1"])
	r_plant
	
	
	# ----------------------
	#  Save sf outputs
	save(s_vert_sf, file = "C:/Users/saraw/Desktop/5_5_18/s_vert_sf.Rdata")
	save(s_fly_sf, file = "C:/Users/saraw/Desktop/5_5_18/s_fly_sf.Rdata")
	save(s_plant_sf, file = "C:/Users/saraw/Desktop/5_5_18/s_plant_sf.Rdata")
	
	

# =============================================================================
#  Convert to rasters to serve as feature input layers in next round of prioritization.
# =============================================================================
	
	# ----------------------
	#  Convert outputs sp objects
	s_vert_sp <- as(s_vert_sf, 'Spatial')
	s_fly_sp <- as(s_fly_sf, 'Spatial')
	s_plant_sp <- as(s_plant_sf, 'Spatial')
	
	# ----------------------
	#  Rasterize
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	# ----------------------
	#  Save rasters
	writeRaster(prior_vert_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r.grd")
	writeRaster(prior_fly_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r.grd")
	writeRaster(prior_plant_feat_r, "C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r.grd")
	