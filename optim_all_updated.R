###############################################################################
#  Run optimization with 'Prioritizr' package for vertebrate species ranges.
#  May 5, 2018.
#  Script to generate optimize selected land units for protection of all input
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
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/sabah_boder.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/sabah_tpa.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/acd_agg_sf_for.Rdata")
	
	# ----------------------
	# Prioritized taxa inputs.
	prior_vert_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r.grd")
	prior_fly_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r.grd")	
	prior_plant_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r.grd")
	
	# ----------------------
	#  Connectivity inputs
	#corr_feat_r <- raster("C:/Users/saraw/Desktop/tmp/corr_out_r_small.grd")
	#elev_conn_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/elev_conn_feat_in.grd")
	prior_cond_r <- raster("
	prior_corr_r <- raster("
	
	# ----------------------
	#  ACD input
	#acd_feat_r <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/acd_feat_in.grd")
	prior_acd_r <- raster("
	
	# ----------------------
	#  All inputs in single raster stack.	
	#all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, elev_conn_feat_r, corr_feat_r, acd_feat_r)
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, prior_cond_feat_r, prior_corr_feat_r, prior_acd_feat_r)
	all_feat_in_si <- stack(prior_vert_feat_r_si, prior_fly_feat_r_si, prior_plant_feat_r_si, prior_cond_feat_r_si, prior_corr_feat_r_si, prior_acd_feat_r_si)

	
# =============================================================================
#  Set up planning units without existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	load(file = "C:/Users/saraw/Desktop/5_5_18/sfi_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/idris_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid.Rdata")
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
# #  Set up planning units with existing TPAs
# # =============================================================================
	
	# # ----------------------
	# #  Planning unit grid with existing TPAs
	# load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid_w_pa.Rdata")
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
	
	
	# sfi_idris_grid <- rbind(sfi_grid_w_pa, idris_grid_w_pa)
	# pu_sf_tmp <- sfi_idris_grid %>%
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
	#  Raster template for outputs
	r_mat <- matrix(0, nrow(all_feat_in[[1]]), ncol(all_feat_in[[1]]))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(all_feat_in[[1]])
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 
	

	
# =============================================================================
#  Scenario 1 - 200,000 from SFI-Idris combined, no constraints
# =============================================================================

	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.6
	p1 <- problem(x = pu_sfi_idris_in, features = all_feat_in_si, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		add_relative_targets(equal_targets) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1 <- solve(p1)
	s1_sf <- st_as_sf(s1) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_area_h <- as.integer(sum(st_area(s1_sf)) * 0.0001)
	s1_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in <- as(s1_sf, 'Spatial')
	
	
	p2 <- problem(x = pu_sfi_idris_in, features = all_feat_in_si, cost_column = "area_h") %>%
		add_max_cover_objective(200000) %>%
		add_binary_decisions() %>%
		add_locked_in_constraints(locked_in) %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s2 <- solve(p2)
	s2_sf <- st_as_sf(s2) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s2_area_h <- as.integer(sum(st_area(s2_sf)) * 0.0001)
	s2_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_all <- as(s2_sf, 'Spatial')
	
	
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, prior_cond_feat_r, prior_corr_feat_r, prior_acd_feat_r)

	
	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.8
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
	# Update 1: 
	upd1_tmp <- as(s3_sf, "Spatial")
	upd1_sf <- st_as_sf(upd1_tmp) %>%	
		mutate(ras_val = 1)
	upd1_sp <- as(upd1_sf, "Spatial")
	upd1_r <- rasterize(upd1_sp, r_template, field = upd1_sp$ras_val)
	writeRaster(upd1_r, file = "C:/Users/saraw/Desktop/5_9_18/upd_map1.tif")
	
	

# =============================================================================
#  Scenario 2 - 200,000 from SFI-Idris combined, neighbor constraints
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.3
	p1 <- problem(x = pu_sfi_idris_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		add_relative_targets(equal_targets) %>%
		add_neighbor_constraints(3) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1 <- solve(p1)
	s1_sf <- st_as_sf(s1) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_area_h <- as.integer(sum(st_area(s1_sf)) * 0.0001)
	s1_area_h
	
	# ----------------------
	#  Set up locked in areas for round 2.
	locked_in <- as(s1_sf, 'Spatial')

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.3
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_neighbor_constraints(3) %>%
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
	locked_in <- as(s2_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem for round 3.
	p3 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_neighbor_constraints(3) %>%
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
	# Update 2: 
	upd2_tmp <- as(s3_sf, "Spatial")
	upd2_sf <- st_as_sf(upd2_tmp) %>%	
		mutate(ras_val = 1)
	upd2_sp <- as(upd1_sf, "Spatial")
	upd2_r <- rasterize(upd2_sp, r_template, field = upd2_sp$ras_val)
	writeRaster(upd2_r, file = "C:/Users/saraw/Desktop/5_5_18/upd_map2_w_si_all_feat.tif")
	
	
	
# =============================================================================
#  Scenario 3 - 200,000 from SFI-Idris combined, boundary length modifier
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.3
	p1 <- problem(x = pu_sfi_idris_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(200000) %>%
		add_relative_targets(equal_targets) %>%
		add_boundary_penalties(0.0000002, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1 <- solve(p1)
	s1_sf <- st_as_sf(s1) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1_area_h <- as.integer(sum(st_area(s1_sf)) * 0.0001)
	s1_area_h
	
	# ----------------------
	#  Set up locked in areas for round 2.
	locked_in <- as(s1_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.3
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in) %>%
		add_locked_out_constraints(locked_out) %>%
		add_boundary_penalties(0.00000001, 0.5) %>%
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
	#  Does not require round 3.

	# ----------------------
	#  Save solution as rasters.
	# Update 3: 
	upd3_tmp <- as(s2_sf, "Spatial")
	upd3_sf <- st_as_sf(upd3_tmp) %>%	
		mutate(ras_val = 1)
	upd3_sp <- as(upd3_sf, "Spatial")
	upd3_r <- rasterize(upd3_sp, r_template, field = upd3_sp$ras_val)
	writeRaster(upd3_r, file = "C:/Users/saraw/Desktop/5_5_18/upd_map3_w_si_all_feat.tif")

	
	
# =============================================================================
#  Plots
# =============================================================================	

	# ----------------------
	#  Plot updated map 1
	upd1_sfi_sf <- st_intersection(upd1_sf, sfi_RT)
	upd1_idris_sf <- st_intersection(upd1_sf, idris_RT)
	
	map1 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = upd1_sf,  fill = "darkorange2",  colour = "darkorange2", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		geom_sf(data = upd1_sfi_sf, fill = "mediumpurple2", colour = "mediumpurple2", alpha = 0.7) +
		geom_sf(data = upd1_idris_sf, fill = "darkslategray3", colour = "darkslategray3", alpha = 0.7) +
		geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map1

	# ----------------------
	#  Plot updated map 2
	upd2_sfi_sf <- st_intersection(upd2_sf, sfi_RT)
	upd2_idris_sf <- st_intersection(upd2_sf, idris_RT)
	
	map2 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = upd2_sf,  fill = "darkorange2",  colour = "darkorange2", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		geom_sf(data = upd2_sfi_sf, fill = "mediumpurple2", colour = "mediumpurple2", alpha = 0.7) +
		geom_sf(data = upd2_idris_sf, fill = "darkslategray3", colour = "darkslategray3", alpha = 0.7) +
		geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map2
	
	# ----------------------
	#  Plot updated map 3
	upd3_sfi_sf <- st_intersection(upd3_sf, sfi_RT)
	upd3_idris_sf <- st_intersection(upd3_sf, idris_RT)
	
	map3 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = upd3_sf,  fill = "darkorange2",  colour = "darkorange2", alpha = 0.7) +
		geom_sf(data = sfi_RT, fill = "mediumpurple4", colour = "mediumpurple4") +
		geom_sf(data = idris_RT, fill = "darkslategray4", colour = "darkslategray4") +
		geom_sf(data = upd3_sfi_sf, fill = "mediumpurple2", colour = "mediumpurple2", alpha = 0.7) +
		geom_sf(data = upd3_idris_sf, fill = "darkslategray3", colour = "darkslategray3", alpha = 0.7) +
		geom_sf(data = deram,  fill = "tomato4",  colour = "tomato4") +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map3
	
	
	
# =============================================================================	
###############################################################################



