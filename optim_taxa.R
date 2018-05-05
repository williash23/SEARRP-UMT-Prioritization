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
# #  Set up planning units WITH existing TPAs
# # =============================================================================
	
	# # ----------------------
	# #  Load planning units with existing TPAs
	# load(file = "C:/Users/saraw/Desktop/5_5_18/sa_grid_w_pa.Rdata")
	# load(file = "C:/Users/saraw/Desktop/5_5_18/tpa_grid.Rdata")
	
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# Rep ########
	
	# ----------------------
	#  Load species ranges that have been replicated by weighting.
	vert_feat_in_rep <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/vert_feat_rep.grd")
	fly_feat_in_rep <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/fly_feat_rep.grd")	
	plant_feat_in_rep <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/feature_inputs/plants_feat_rep.grd")
		

	# ----------------------
	#  Set up problem for vertebrates.
	p_vert_rep <- problem(x = pu_in, features = vert_feat_in_rep, cost_column = "area_h") %>%
		add_max_features_objective(300000) %>%
		add_relative_targets(0.3) %>%
		add_neighbor_constraints(2) %>%
		#add_boundary_penalties(penalty = 0.00000001) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve vertebrates problem.
	s_vert_rep <- solve(p_vert_rep)
	s_vert_rep_sf <- st_as_sf(s_vert_rep) %>%
		dplyr::filter(solution_1 == 1)
	vert_area_h_rep <- as.integer(sum(st_area(s_vert_rep_sf)) * 0.0001)
	r_vert <- feat_rep(p_vert_rep, s_vert_rep[, "solution_1"])
	
	# Single #########
	# ----------------------
	#  Set up problem for vertebrates.
	p_vert_rep <- problem(x = pu_in, features = vert_feat_in_single, cost_column = "area_h") %>%
		add_max_features_objective(350000) %>%
		add_relative_targets(0.1) %>%
		add_feature_weights(vert_rep_weight) %>% 
		add_boundary_penalties(penalty = 0.00000001) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve vertebrates problem.
	s_vert_rep <- solve(p_vert_rep)
	s_vert_rep_sf <- st_as_sf(s_vert_rep) %>%
		dplyr::filter(solution_1 == 1)
	vert_area_h_rep <- as.integer(sum(st_area(s_vert_rep_sf)) * 0.0001)
	vert_area_h_rep
	
	r_vert <- feat_rep(p_vert_rep, s_vert_rep[, "solution_1"])
	r_vert
	
	
	# ----------------------
	#  Set up problem for butterflies.
	p_fly_rep <- problem(x = pu_in, features = fly_feat_in_rep, cost_column = "area_h") %>%
		add_max_utility_objective(300000) %>%
		#add_max_cover_objective(300000) %>%
		#add_feature_weights(w) %>%
		#add_boundary_penalties(500, 0.5) %>%
		add_neighbor_constraints(2) %>%
		#add_contiguity_constraints() %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve butterflies problem.
	s_fly_rep <- solve(p_fly_rep)
	s_fly_rep_sf <- st_as_sf(s_fly_rep) %>%
		dplyr::filter(solution_1 == 1)
	fly_area_h_rep <- as.integer(sum(st_area(s_fly_rep_sf)) * 0.0001)
	

	
	
	
	# ----------------------
	#  Set up problem for plants.
	p_plant_rep <- problem(x = pu_in, features = plant_feat_in_rep, cost_column = "area_h") %>%
		add_max_utility_objective(300000) %>%
		#add_max_cover_objective(300000) %>%
		#add_feature_weights(w) %>%
		#add_boundary_penalties(500, 0.5) %>%
		#add_neighbor_constraints(2) %>%
		#add_contiguity_constraints() %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve plants problem.
	s_plant_rep <- solve(p_plant_rep)
	s_plant_rep_sf <- st_as_sf(s_plant_rep) %>%
		dplyr::filter(solution_1 == 1)
	plant_area_h_rep <- as.integer(sum(st_area(s_plant_rep_sf)) * 0.0001)
	
	
	
	
	# ----------------------
	#  Save outputs.
	save(s_vert_rep_sf, file = "C:/Users/saraw/Desktop/tmp/s_vert_rep_sf.Rdata")
	save(s_fly_rep_sf, file = "C:/Users/saraw/Desktop/tmp/s_fly_rep_sf.Rdata")
	save(s_plant_rep_sf, file = "C:/Users/saraw/Desktop/tmp/s_plant_rep_sf.Rdata")
	
	# ----------------------
	#  Convert outputs to raster format for input in next round of prioritization.
	s_vert_sp <- as(s_vert_rep_sf, 'Spatial')
	s_fly_sp <- as(s_fly_rep_sf, 'Spatial')
	s_plant_sp <- as(s_plant_rep_sf, 'Spatial')
	
	prior_vert_feat_r <- rasterize(s_vert_sp, r_template, field = s_vert_sp$solution_1)
	prior_fly_feat_r <- rasterize(s_fly_sp, r_template, field = s_fly_sp$solution_1)
	prior_plant_feat_r <- rasterize(s_plant_sp, r_template, field = s_plant_sp$solution_1)
	
	# ----------------------
	#  Save outputs.
	writeRaster(prior_vert_feat_r, "C:/Users/saraw/Desktop/tmp/prior_vert_feat_r.grd")
	writeRaster(prior_fly_feat_r, "C:/Users/saraw/Desktop/tmp/prior_fly_feat_r.grd")
	writeRaster(prior_plant_feat_r, "C:/Users/saraw/Desktop/tmp/prior_plant_feat_r.grd")
	
	

	
	
	
	# ----------------------
	#  Save plot to see planning units selected over time.
	 solution_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = pts_wt_sf_plot_sc, aes(colour = wt_sc), alpha = 0.4) +
		#scale_colour_distiller(type = "seq", palette = "Reds", direction = 1) +		
		geom_sf(data = acd_agg_sf_for, colour = "#1B792F", fill = "#1B792F", alpha = 0.3) +
		#geom_sf(data = sa_grid_sf, colour = "black", fill = "transparent") +
		geom_sf(data = ssk_pa_near_sf_clust_1k, fill = "darkorange2", colour = "grey50", alpha = 0.7) +
		#geom_sf(data = all_sol_sf_lg,  aes(fill = factor(group))) +
		#geom_sf(data = all_sol_sf_clust,  aes(fill = factor(group))) +
		geom_sf(data = s_fly_rep_sf,  fill = "darkred", alpha = 0.7) +
		#scale_fill_distiller(type = "seq", palette = "Reds", direction = 1) +	
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	solution_p

	
	
	
	
	
	
	
	
	