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

	
	
##############################################################################
############################## SCENARIO 1 #####################################
##############################################################################


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges, connectivity and carbon 
#   for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1 <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
	p_fly1 <- problem(x = pu_in, features = fly_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.2) %>%
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
	p_plant1 <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
	p_cond1 <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.325) %>%
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
	p_acd1 <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.2) %>%
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
	p_corr1 <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		#add_max_utility_objective(350000) %>%
		add_max_features_objective(350000) %>%
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
	

	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_sp <- as(s_vert1_sf, 'Spatial')
	s_fly_sp <- as(s_fly1_sf, 'Spatial')
	s_plant_sp <- as(s_plant1_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_u <- st_union(s_cond1_sf)
	s_acd1_u <- st_union(s_acd1_sf)
	s_corr1_u <- st_union(s_corr1_sf)
	# s_corr1_u_tmp1 <- st_cast(st_union(s_corr1_sf), 'POLYGON')
		# s_corr1_u_tmp2 <- st_geometry(s_corr1_u_tmp1)
		# s_corr1_u_tmp3 <- as(s_corr1_u_tmp2, 'Spatial')
		# s_corr1_u_tmp4 <- st_as_sf(s_corr1_u_tmp3)
		# s_corr1_u_tmp5 <- s_corr1_u_tmp4[-8,]
	# s_corr1_u <- st_cast(s_corr1_u_tmp5, 'MULTIPOLYGON')
	
	# ----------------------
	#  Convert to sp objects
	s_cond_sp <- as(s_cond1_u, 'Spatial')
	s_acd_sp <- as(s_acd1_u, 'Spatial')
	s_corr_sp <- as(s_corr1_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_feat_r <- raster::mask(elev_conn_feat_r, s_cond_sp)
	prior_acd_feat_r <- raster::mask(acd_feat_r, s_acd_sp)
	prior_corr_feat_r <- raster::mask(corr_feat_r, s_corr_sp)
	
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, prior_cond_feat_r, prior_corr_feat_r, prior_acd_feat_r)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.4
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(equal_targets) %>%
		#add_locked_in_constraints(locked_in_all) %>%
		#add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.00000002, 0.5) %>%
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
		add_max_utility_objective(350000) %>%
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
	
	# ----------------------
	#  Save solution as rasters.
	# Update scenario 1: 
	upd1_tmp <- as(s3_sf, "Spatial")
	upd1_sf <- st_as_sf(upd1_tmp) %>%	
		mutate(ras_val = 1)
	upd1_sp <- as(upd1_sf, "Spatial")
	upd1_r <- rasterize(upd1_sp, r_template, field = upd1_sp$ras_val)
	#writeRaster(upd1_r, file = "C:/Users/saraw/Desktop/5_9_18/upd_map1.tif")
	
	# ----------------------
	#  Map 
	map1 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		#geom_sf(data = s_vert1_sf, fill = "turquoise4", colour = "turquoise4", alpha = 0.5) +
		#geom_sf(data = s_fly1_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.5) +
		#geom_sf(data = s_plant1_sf, fill = "springgreen4", colour = "springgreen4", alpha = 0.5) +
		#geom_sf(data = s2_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		#geom_sf(data = s1_sf, fill = "darkorange2", colour = "darkorange2", alpha = 0.8) +
		#geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map1

	
	
	# ----------------------
	#  Map individual inputs
	map1a <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = upd1_sf,  fill = "darkorange2",  colour = "darkorange2", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		geom_sf(data = s_vert1_sf, fill = "turquoise4", colour = "turquoise4", alpha = 0.5) +
		geom_sf(data = s_fly1_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.5) +
		geom_sf(data = s_plant1_sf, fill = "springgreen4", colour = "springgreen4", alpha = 0.5) +
		#geom_sf(data = upd1_sfi_sf, fill = "mediumpurple2", colour = "mediumpurple2", alpha = 0.7) +
		#geom_sf(data = upd1_idris_sf, fill = "darkslategray3", colour = "darkslategray3", alpha = 0.7) +
		geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map1a
	
	
	
	
	
	
	
##############################################################################
############################## SCENARIO 2 #####################################
##############################################################################


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges, connectivity and carbon for whole of Sabah
# =============================================================================

	# ----------------------
	#  Set up problem 1 for vertebrates
	p_vert1 <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
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
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
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
	p_plant1 <- problem(x = pu_in, features = plant_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
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
	p_cond1 <- problem(x = pu_in, features = elev_conn_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
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
	p_acd1 <- problem(x = pu_in, features = acd_feat_r, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
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
	p_corr1 <- problem(x = pu_in, features = corr_feat_r, cost_column = "area_h") %>%
		#add_max_utility_objective(350000) %>%
		add_max_features_objective(350000) %>%
		#add_max_features_objective(200000 + tpa_area_h) %>%
		add_relative_targets(0.4) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve
	s_corr1 <- solve(p_corr1)
	s_corr1_sf <- st_as_sf(s_corr1) %>%
		dplyr::filter(solution_1 == 1)
	s_corr1_area_h <- as.integer(sum(st_area(s_corr1_sf)) * 0.0001)
	s_corr1_area_h
	

	# ----------------------
	#  Convert prioritized species ranges into sp objects
	s_vert_sp <- as(s_vert1_sf, 'Spatial')
	s_fly_sp <- as(s_fly1_sf, 'Spatial')
	s_plant_sp <- as(s_plant1_sf, 'Spatial')

	# ----------------------
	#  Convert prioritized species ranges into rasters to serve as input for next step.
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	# ----------------------
	#  Union prioritized connectivity and carbon into single sf object polygons
	s_cond1_u <- st_union(s_cond1_sf)
	s_acd1_u <- st_union(s_acd1_sf)
	s_corr1_u_tmp1 <- st_cast(st_union(s_corr1_sf), 'POLYGON')
		s_corr1_u_tmp2 <- st_geometry(s_corr1_u_tmp1)
		s_corr1_u_tmp3 <- as(s_corr1_u_tmp2, 'Spatial')
		s_corr1_u_tmp4 <- st_as_sf(s_corr1_u_tmp3)
		s_corr1_u_tmp5 <- s_corr1_u_tmp4[-8,]
	s_corr1_u <- st_cast(s_corr1_u_tmp5, 'MULTIPOLYGON')
	
	# ----------------------
	#  Convert to sp objects
	s_cond_sp <- as(s_cond1_u, 'Spatial')
	s_acd_sp <- as(s_acd1_u, 'Spatial')
	s_corr_sp <- as(s_corr1_u, 'Spatial')
	
	# ----------------------
	#  Mask prioritized area of original rasters so that we maintain continuous values.
	prior_cond_feat_r <- raster::mask(elev_conn_feat_r, s_cond_sp)
	prior_acd_feat_r <- raster::mask(acd_feat_r, s_acd_sp)
	prior_corr_feat_r <- raster::mask(corr_feat_r, s_corr_sp)
	
	
	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, prior_cond_feat_r, prior_corr_feat_r, prior_acd_feat_r)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.8
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in_all) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.00000002, 0.5) %>%
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
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.00000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	s3_sf <- st_as_sf(s3) %>%
		dplyr::filter(solution_1 == 1)
	s3_area_h <- as.integer(sum(st_area(s3_sf)) * 0.0001)
	s3_area_h
	
	
	# ----------------------
	#  Save solution as rasters.
	# Update scenario 1: 
	upd2_tmp <- as(s3_sf, "Spatial")
	upd2_sf <- st_as_sf(upd2_tmp) %>%	
		mutate(ras_val = 1)
	upd2_sp <- as(upd2_sf, "Spatial")
	upd2_r <- rasterize(upd2_sp, r_template, field = upd2_sp$ras_val)
	writeRaster(upd2_r, file = "C:/Users/saraw/Desktop/5_9_18/upd_map2.tif")
	
	
	# ----------------------
	#  Map 
	map2 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		#geom_sf(data = s_vert1_sf, fill = "turquoise4", colour = "turquoise4", alpha = 0.5) +
		#geom_sf(data = s_fly1_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.5) +
		#geom_sf(data = s_plant1_sf, fill = "springgreen4", colour = "springgreen4", alpha = 0.5) +
		geom_sf(data = s3_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		geom_sf(data = s1_sf, fill = "darkorange2", colour = "darkorange2", alpha = 0.8) +
		geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map2
	
	