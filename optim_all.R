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
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/mainland_sabah_boder.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/sabah_tpa.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/5_5_18/acd_agg_sf_for.Rdata")
	
	# ----------------------
	# Prioritized taxa inputs.
	prior_vert_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_vert_feat_r.grd")
	prior_fly_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_fly_feat_r.grd")	
	prior_plant_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/prior_plant_feat_r.grd")
	
	# ----------------------
	#  Connectivity inputs
	corr_feat_r <- raster("C:/Users/saraw/Desktop/tmp/corr_out_r_small.grd")
	elev_conn_feat_r <- raster("C:/Users/saraw/Desktop/5_5_18/elev_conn_feat_in.grd")
	
	# ----------------------
	#  ACD input
	acd_feat_r <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/acd_feat_in.grd")
	
	# ----------------------
	#  All inputs in single raster stack.	
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, elev_conn_feat_r, corr_feat_r, acd_feat_r)

	
	
# =============================================================================
#  Set up planning units without existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grid
	load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid.Rdata")
		
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
	
	

# # =============================================================================
# #  Set up planning units with existing TPAs
# # =============================================================================
	
	# # ----------------------
	# #  Planning unit grid with existing TPAs
	# load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid_w_pa.Rdata")

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
#  Scenario 1
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.4
	p1 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
	#  Set up problem for round 2.
	locked_in <- as(s1_sf, 'Spatial') 

	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
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
	#  Save solution as rasters.
	# Prelim 1: No constraints
	prelim1_tmp <- as(s2_sf, "Spatial")
	prelim1_sf <- st_as_sf(prelim1_tmp) %>%	
		mutate(ras_val = 1)
	prelim1_sp <- as(prelim1_sf, "Spatial")
	prelim1_r <- rasterize(prelim1_sp, r_template, field = prelim1_sp$ras_val)
	writeRaster(prelim1_r, file = "C:/Users/saraw/Desktop/5_5_18/prelim_map1.tif")
	
	
	
# =============================================================================
#  Scenario 2
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.4
	p1 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(equal_targets) %>%
		#add_boundary_penalties(0.0000002, 0.5) %>%
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
	#  Set up problem for round 2.
	locked_in <- as(s1_sf, 'Spatial') 
	
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
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
	#  Save solution as rasters.
	# Prelim 2: Min 2000 ha clusters
	prelim2_tmp <- as(s2_sf, "Spatial")
	prelim2_sf <- st_as_sf(prelim2_tmp) %>%	
		mutate(ras_val = 1)
	prelim2_sp <- as(prelim2_sf, "Spatial")
	prelim2_r <- rasterize(prelim2_sp, r_template, field = prelim2_sp$ras_val)
	writeRaster(prelim2_r, file = "C:/Users/saraw/Desktop/5_5_18/prelim_map2.tif")

	
	
# =============================================================================
#  Scenario 3
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	equal_targets <- 0.4
	p1 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
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
	#  Does not require round 2.
	
	# ----------------------
	#  Save solution as rasters.
	# Prelim 3: Boundary constraint
	prelim3_tmp <- as(s1_sf, "Spatial")
	prelim3_sf <- st_as_sf(prelim3_tmp) %>%	
		mutate(ras_val = 1)
	prelim3_sp <- as(prelim3_sf, "Spatial")
	prelim3_r <- rasterize(prelim3_sp, r_template, field = prelim3_sp$ras_val)
	writeRaster(prelim3_r, file = "C:/Users/saraw/Desktop/5_5_18/prelim_map3.tif")



# =============================================================================
#  Scenario 1b
# =============================================================================
	
	# ----------------------
	#  Set up problem for round 1.
	unequal_targets <- c(0.5, 0.5, 0.5, 0.25, 0.25, 0.25)
	p1b <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(unequal_targets) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	# Solve round 1.
	s1b <- solve(p1b)
	s1b_sf <- st_as_sf(s1b) %>%
		dplyr::filter(solution_1 == 1) %>%
		mutate(locked_in = 1)
	s1b_area_h <- as.integer(sum(st_area(s1b_sf)) * 0.0001)
	s1b_area_h
	
	# ----------------------
	#  Set up problem for round 2.
	locked_in <- as(s1b_sf, 'Spatial') 

	p2b <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(350000) %>%
		add_locked_in_constraints(locked_in) %>%
		add_binary_decisions() %>%
		add_gurobi_solver(time_limit = 30)
		
	# ----------------------
	#  Solve round 2.
	s2b <- solve(p2b)
	s2b_sf <- st_as_sf(s2b) %>%
		dplyr::filter(solution_1 == 1)
	s2b_area_h <- as.integer(sum(st_area(s2b_sf)) * 0.0001)
	s2b_area_h
	
	# ----------------------
	#  Save solution as rasters.
	# Prelim 1: No constraints
	prelim1b_tmp <- as(s2b_sf, "Spatial")
	prelim1b_sf <- st_as_sf(prelim1b_tmp) %>%	
		mutate(ras_val = 1)
	prelim1b_sp <- as(prelim1b_sf, "Spatial")
	prelim1b_r <- rasterize(prelim1b_sp, r_template, field = prelim1b_sp$ras_val)
	writeRaster(prelim1b_r, file = "C:/Users/saraw/Desktop/5_5_18/prelim_map1b.tif")
	


	
# =============================================================================
#  Plots
# =============================================================================	

	# ----------------------
	#  Plot preliminary map 1
	 map1 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		#geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = prelim1_sf, fill = "darkorange3", colour = "darkorange3", alpha = 0.7) +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map1
	
	# ----------------------
	#  Plot preliminary map 2
	 map2 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = prelim2_sf, fill = "darkorange3", colour = "darkorange3", alpha = 0.7) +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map2
	
	# ----------------------
	#  Plot preliminary map 3
	 map3 <- ggplot() +
		geom_sf(data = mainland_sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = acd_agg_sf_for, colour = "#7d9983", fill = "#7d9983") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = prelim3_sf, fill = "darkorange3", colour = "darkorange3", alpha = 0.7) +
			coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	map3
	
	
	
# =============================================================================	
###############################################################################
