###############################################################################
#  Prioritization - Prepare species ranges including ranges across existing Protected Areas
#   inputs in chunks
#  October 12, 2018
#  Prioritize each taxon individually - in percentage chunks, so that species range inputs
#   becomes coverage across all of Sabah and value is order in which that cell was
#   prioritized. 
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
library(gurobi)
library(prioritizr)
library(ggplot2)


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
#  Set up planning units WITHOUT existing TPAs
# =============================================================================
	
	# ----------------------
	#  Planning unit grids
	load(file = "planning_unit_grids/sa_grid.Rdata")

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

	
	
# =============================================================================
#  Run prioritzation with no BLM constraint
# =============================================================================	
	

# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================
	
	
		
		start <- Sys.time()
		start

		chunks <- seq(0.1, 2, by = 0.1)
	
		rev_vert_rep_weight <- vert_rep_weight_w_pa_range
		rev_vert_rep_weight[rev_vert_rep_weight == 4] <- 10
		rev_vert_rep_weight[rev_vert_rep_weight == 2] <- 4
		rev_vert_rep_weight[rev_vert_rep_weight == 10] <- 2
	
		v_rel_tar <- chunks[1]*(1/rev_vert_rep_weight)
		
		# ----------------------
		#  Vertebrates
		p_vert_1 <- problem(x = pu_in, features = vert_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(v_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		# ----------------------
		#  Solve round 2.
		s_vert_1 <- solve(p_vert_1)
		s_vert_1_sf  <- st_as_sf(s_vert_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_vert_1_area_h <- as.integer(sum(st_area(s_vert_1_sf)) * 0.0001)
		s_vert_1_area_h
		locked_in <- as(s_vert_1_sf, "Spatial")
		
		
		for(i in 2:length(chunks)){
			
			v_rel_tar <- chunks[i]*(1/rev_vert_rep_weight)
			
			# ----------------------
			#  Set up problem 2 for vertebrates
			p_vert_tmp <- problem(x = pu_in, features = vert_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(v_rel_tar) %>%				
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_vert_tmp <- solve(p_vert_tmp)
			s_vert_tmp_sf <- st_as_sf(s_vert_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
				
			nam <- paste("s_vert", i, "sf", sep = "_")
			assign(nam, s_vert_tmp_sf)
			
			# ----------------------
			#  Set up new locked in area.
			locked_in <- as(s_vert_tmp_sf, "Spatial")
			
			}
			
		new_vert_1_sf <- s_vert_1_sf
		
		for(i in 2:length(chunks)){
			
			tmp1 <- get(paste("s_vert", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_vert", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_vert", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
			
		all_vert_sf_w_pa_range <- rbind(new_vert_1_sf, new_vert_2_sf, new_vert_3_sf, 
			new_vert_4_sf, new_vert_5_sf, new_vert_6_sf, new_vert_7_sf, 
			new_vert_8_sf, new_vert_9_sf, new_vert_10_sf, new_vert_11_sf, 
			new_vert_12_sf, new_vert_13_sf, new_vert_14_sf, new_vert_15_sf, 
			new_vert_16_sf, new_vert_17_sf, new_vert_18_sf, new_vert_19_sf, 
			new_vert_20_sf) %>%
			mutate(ras_val = 1 - (order/length(chunks)))
		all_vert_sp_w_pa_range <- as(all_vert_sf_w_pa_range, "Spatial")
		
		save(all_vert_sf_w_pa_range, file = "all_vert_sf_w_pa_range.Rdata")
		prior_vert_feat_r_w_pa_range <- rasterize(all_vert_sp_w_pa_range, r_template, 
			field = all_vert_sp_w_pa_range$ras_val)
		writeRaster(prior_vert_feat_r_w_pa_range, file = "chunks/prior_vert_feat_r_w_pa_range.grd")
	
		end <- Sys.time()
		el1 <- end - start
		el1
		
		
		start <- Sys.time()
		start
		
		chunks <- seq(0.1, 1, by = 0.1)
		
		rev_fly_rep_weight <- fly_rep_weight_w_pa_range
		rev_fly_rep_weight[rev_fly_rep_weight == 2] <- 10
		rev_fly_rep_weight[rev_fly_rep_weight == 1] <- 2
		rev_fly_rep_weight[rev_fly_rep_weight == 10] <- 1
	
		f_rel_tar <- chunks[1]*(1/rev_fly_rep_weight)
		
		# ----------------------
		#  Butterflies
		p_fly_1 <- problem(x = pu_in, features = fly_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(f_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver() 

		# ----------------------
		#  Solve round 2.
		s_fly_1 <- solve(p_fly_1)
		
		s_fly_1_sf  <- st_as_sf(s_fly_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_fly_1_area_h <- as.integer(sum(st_area(s_fly_1_sf)) * 0.0001)
		s_fly_1_area_h
		locked_in <- as(s_fly_1_sf, "Spatial")
		
		for(i in 2:length(chunks)){
		
			f_rel_tar <- chunks[i]*(1/rev_fly_rep_weight)
		
			# ----------------------
			#  Set up problem 2 for butterflies
			p_fly_tmp <- problem(x = pu_in, features = fly_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(f_rel_tar) %>% 
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_fly_tmp <- solve(p_fly_tmp)
			s_fly_tmp_sf <- st_as_sf(s_fly_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
			
			nam <- paste("s_fly", i, "sf", sep = "_")
			assign(nam, s_fly_tmp_sf)		
			
			# ----------------------
			#  Set up new locked in area.
			locked_in <- as(s_fly_tmp_sf, "Spatial")
			}
		
		new_fly_1_sf <- s_fly_1_sf
		
		for(i in 2:length(chunks)){
			
			tmp1 <- get(paste("s_fly", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_fly", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_fly", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
		
		all_fly_sf_w_pa_range <- rbind(new_fly_1_sf, new_fly_2_sf, new_fly_3_sf, 
			new_fly_4_sf, new_fly_5_sf, new_fly_6_sf, new_fly_7_sf, 
			new_fly_8_sf, new_fly_9_sf, new_fly_10_sf) %>%
			mutate(ras_val = 1 - (order/length(chunks)))
		all_fly_sp_w_pa_range <- as(all_fly_sf_w_pa_range, "Spatial")
		
		save(all_fly_sf_w_pa_range, file = "all_fly_sf_w_pa_range.Rdata")
		prior_fly_feat_r_w_pa_range <- rasterize(all_fly_sp_w_pa_range, r_template, 
			field = all_fly_sp_w_pa_range$ras_val)
		writeRaster(prior_fly_feat_r_w_pa_range, file = "chunks/prior_fly_feat_r_w_pa_range.grd")
		
		end <- Sys.time()
		el2 <- end - start
		el2
		
		
		
		start <- Sys.time()
		start
		
		chunks <- seq(0.1, 1, by = 0.1)
		
		rev_plant_rep_weight <- plant_rep_weight_w_pa_range
		rev_plant_rep_weight[rev_plant_rep_weight == 2] <- 10
		rev_plant_rep_weight[rev_plant_rep_weight == 1] <- 2
		rev_plant_rep_weight[rev_plant_rep_weight == 10] <- 1
	
		p_rel_tar <- chunks[1]*(1/rev_plant_rep_weight)
		
		# ----------------------
		#  Plants
		p_plant_1 <- problem(x = pu_in, features = plant_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(p_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		# ----------------------
		#  Solve round 2.
		s_plant_1 <- solve(p_plant_1)
		s_plant_1_sf  <- st_as_sf(s_plant_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_plant_1_area_h <- as.integer(sum(st_area(s_plant_1_sf)) * 0.0001)
		s_plant_1_area_h
		
		locked_in <- as(s_plant_1_sf, "Spatial")
		
		for(i in 2:length(chunks)){
		
			p_rel_tar <- chunks[i]*(1/rev_plant_rep_weight)
			
			# ----------------------
			#  Set up problem 2 for butterflies
			p_plant_tmp <- problem(x = pu_in, features = plant_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(p_rel_tar) %>%
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_plant_tmp <- solve(p_plant_tmp)
			s_plant_tmp_sf <- st_as_sf(s_plant_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
			
			nam <- paste("s_plant", i, "sf", sep = "_")
			assign(nam, s_plant_tmp_sf)		
			
			# ----------------------
			#  Set up new locked in area.
			locked_in <- as(s_plant_tmp_sf, "Spatial")
			}
	
		new_plant_1_sf <- s_plant_1_sf
			
		for(i in 2:length(chunks)){
			
			tmp1 <- get(paste("s_plant", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_plant", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_plant", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
		
		all_plant_sf_w_pa_range <- rbind(new_plant_1_sf, new_plant_2_sf, new_plant_3_sf, 
			new_plant_4_sf, new_plant_5_sf, new_plant_6_sf, new_plant_7_sf, 
			new_plant_8_sf, new_plant_9_sf, new_plant_10_sf) %>%
			mutate(ras_val = 1 - (order/length(chunks)))
		all_plant_sp_w_pa_range <- as(all_plant_sf_w_pa_range, "Spatial")
		
		save(all_plant_sf_w_pa_range, file = "all_plant_sf_w_pa_range.Rdata")
		prior_plant_feat_r_w_pa_range <- rasterize(all_plant_sp_w_pa_range, r_template, 
			field = all_plant_sp_w_pa_range$ras_val)
		writeRaster(prior_plant_feat_r_w_pa_range, file = "chunks/prior_plant_feat_r_w_pa_range.grd")
	
		end <- Sys.time()
		el3 <- end - start
		el3
	
	
	
		
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(prior_vert_feat_r_w_pa_range, prior_fly_feat_r_w_pa_range, 
		prior_plant_feat_r_w_pa_range, elev_conn_feat_in_w_pa_range, 
		corr_conn_feat_in_w_pa_range, acd_feat_in_w_pa_range)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.3
	p2 <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(423500) %>%
		add_relative_targets(equal_targets) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
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
		add_max_utility_objective(423500) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3 <- solve(p3)
	s3_sf <- st_as_sf(s3) %>%
		dplyr::filter(solution_1 == 1)
	s3_area_h <- as.integer(sum(st_area(s3_sf)) * 0.0001)
	s3_area_h
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen5_tmp <- as(s3_sf, "Spatial")
	scen5_chunks_sf <- st_as_sf(scen5_tmp) %>%	
		mutate(ras_val = 1)
	scen5_chunks_out <- s3
	save(scen5_chunks_sf, file = "scenario_outputs/scen5/scen5_chunks_sf.Rdata")
	save(scen5_chunks_out, file = "scenario_outputs/scen5/scen5_chunks_out.Rdata")
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
# =============================================================================
#  Run prioritzation with BLM constraint
# =============================================================================	
	


# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================

		start <- Sys.time()
		start
		# ----------------------
		#  Set up problem 2 for vertebrates
		p_vert_blm_1 <- problem(x = pu_in, features = vert_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(0.1) %>%
			add_feature_weights(vert_rep_weight) %>% 
			add_boundary_penalties(0.0000001, 0.5) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		# ----------------------
		#  Solve round 2.
		s_vert_blm_1 <- solve(p_vert_blm_1)
		s_vert_blm_1_sf  <- st_as_sf(s_vert_blm_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_vert_blm_1_area_h <- as.integer(sum(st_area(s_vert_blm_1_sf)) * 0.0001)
		s_vert_blm_1_area_h
		locked_in_blm <- as(s_vert_blm_1_sf, "Spatial")
		
		for(i in 2:10){
		
			# ----------------------
			#  Set up problem 2 for vertebrates
			p_vert_blm_tmp <- problem(x = pu_in, features = vert_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(i/10) %>%
				add_feature_weights(vert_rep_weight) %>% 
				add_boundary_penalties(0.0000001, 0.5) %>%
				add_locked_in_constraints(locked_in_blm) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_vert_blm_tmp <- solve(p_vert_blm_tmp)
			s_vert_blm_tmp_sf <- st_as_sf(s_vert_blm_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
				
			nam <- paste("s_vert_blm", i, "sf", sep = "_")
			assign(nam, s_vert_blm_tmp_sf)
			
			# ----------------------
			#  Set up new locked in area.
			locked_in_blm <- as(s_vert_tmp_blm_sf, "Spatial")
			}
			
		new_vert_blm_1_sf <- s_vert_blm_1_sf
		
		for(i in 2:10){
			
			tmp1 <- get(paste("s_vert_blm", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_vert_blm", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_vert_blm", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
			
		all_vert_blm_sf <- rbind(new_vert_blm_1_sf, new_vert_blm_2_sf, new_vert_blm_3_sf, 
			new_vert_blm_4_sf, new_vert_blm_5_sf, new_vert_blm_6_sf, new_vert_blm_7_sf, 
			new_vert_blm_8_sf, new_vert_blm_9_sf, new_vert_blm_10_sf) %>%
			mutate(ras_val = 1 - (order/10))
		all_vert_blm_sp <- as(all_vert_blm_sf, "Spatial")
		
		prior_vert_blm_feat_r <- rasterize(all_vert_blm_sp, r_template, field = all_vert_blm_sp$ras_val)
		writeRaster(prior_vert_blm_feat_r, file = "prior_vert_blm_feat_r.tif")
	
		end <- Sys.time()
		el1 <- end - start
		el1
		
		
		start <- Sys.time()
		start
		
		# ----------------------
		#  Set up problem 2 for butterflies
		p_fly_blm_1 <- problem(x = pu_in, features = fly_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(0.1) %>%
			add_feature_weights(fly_rep_weight) %>% 
			add_boundary_penalties(0.0000001, 0.5) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		# ----------------------
		#  Solve round 2.
		s_fly_blm_1 <- solve(p_fly_blm_1)
		s_fly_blm_1_sf  <- st_as_sf(s_fly_blm_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_fly_blm_1_area_h <- as.integer(sum(st_area(s_fly_blm_1_sf)) * 0.0001)
		s_fly_blm_1_area_h
		locked_in_blm <- as(s_fly_blm_1_sf, "Spatial")
		
		for(i in 2:10){
		
			# ----------------------
			#  Set up problem 2 for butterflies
			p_fly_blm_tmp <- problem(x = pu_in, features = fly_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(i/10) %>%
				add_feature_weights(fly_rep_weight) %>% 
				add_boundary_penalties(0.0000001, 0.5) %>%
				add_locked_in_constraints(locked_in_blm) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_fly_blm_tmp <- solve(p_fly_blm_tmp)
			s_fly_blm_tmp_sf <- st_as_sf(s_fly_blm_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
				
			nam <- paste("s_fly_blm", i, "sf", sep = "_")
			assign(nam, s_fly_blm_tmp_sf)
			
			# ----------------------
			#  Set up new locked in area.
			locked_in_blm <- as(s_fly_tmp_blm_sf, "Spatial")
			}
			
		new_fly_blm_1_sf <- s_fly_blm_1_sf
		
		for(i in 2:10){
			
			tmp1 <- get(paste("s_fly_blm", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_fly_blm", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_fly_blm", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
			
		all_fly_blm_sf <- rbind(new_fly_blm_1_sf, new_fly_blm_2_sf, new_fly_blm_3_sf, 
			new_fly_blm_4_sf, new_fly_blm_5_sf, new_fly_blm_6_sf, new_fly_blm_7_sf, 
			new_fly_blm_8_sf, new_fly_blm_9_sf, new_fly_blm_10_sf) %>%
			mutate(ras_val = 1 - (order/10))
		all_fly_blm_sp <- as(all_fly_blm_sf, "Spatial")
		
		prior_fly_blm_feat_r <- rasterize(all_fly_blm_sp, r_template, field = all_fly_blm_sp$ras_val)
		writeRaster(prior_fly_blm_feat_r, file = "prior_fly_blm_feat_r.tif")
	
		
		end <- Sys.time()
		el2 <- end - start
	
		
		start <- Sys.time()
		start
		
		# ----------------------
		#  Set up problem 2 for plants
		p_plant_blm_1 <- problem(x = pu_in, features = plant_feat_in, cost_column = "area_h") %>%
			add_min_set_objective() %>%
			add_relative_targets(0.1) %>%
			add_feature_weights(plant_rep_weight) %>% 
			add_boundary_penalties(0.0000001, 0.5) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		# ----------------------
		#  Solve round 2.
		s_plant_blm_1 <- solve(p_plant_blm_1)
		s_plant_blm_1_sf  <- st_as_sf(s_plant_blm_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			dplyr::mutate(order = 1)
			
		s_plant_blm_1_area_h <- as.integer(sum(st_area(s_plant_blm_1_sf)) * 0.0001)
		s_plant_blm_1_area_h
		locked_in_blm <- as(s_plant_blm_1_sf, "Spatial")
		
		for(i in 2:10){
		
			# ----------------------
			#  Set up problem 2 for plants
			p_plant_blm_tmp <- problem(x = pu_in, features = plant_feat_in, cost_column = "area_h") %>%
				add_min_set_objective() %>%
				add_relative_targets(i/10) %>%
				add_feature_weights(plant_rep_weight) %>% 
				add_boundary_penalties(0.0000001, 0.5) %>%
				add_locked_in_constraints(locked_in_blm) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			# ----------------------
			#  Solve round 2.
			s_plant_blm_tmp <- solve(p_plant_blm_tmp)
			s_plant_blm_tmp_sf <- st_as_sf(s_plant_blm_tmp) %>%
				dplyr::filter(solution_1 == 1) %>%
				dplyr::mutate(order = i)
				
			nam <- paste("s_plant_blm", i, "sf", sep = "_")
			assign(nam, s_plant_blm_tmp_sf)
			
			# ----------------------
			#  Set up new locked in area.
			locked_in_blm <- as(s_plant_tmp_blm_sf, "Spatial")
			}
			
		new_plant_blm_1_sf <- s_plant_blm_1_sf
		
		for(i in 2:10){
			
			tmp1 <- get(paste("s_plant_blm", i-1, "sf", sep = "_"))
			tmp2 <- get(paste("s_plant_blm", i, "sf", sep = "_"))
			tmp3 <- st_intersection(tmp1, tmp2)
			tmp4 <- st_erase(tmp2, tmp3)
			
			nam <- paste("new_plant_blm", i, "sf", sep = "_")
			assign(nam, tmp4)
			}
			
		all_plant_blm_sf <- rbind(new_plant_blm_1_sf, new_plant_blm_2_sf, new_plant_blm_3_sf, 
			new_plant_blm_4_sf, new_plant_blm_5_sf, new_plant_blm_6_sf, new_plant_blm_7_sf, 
			new_plant_blm_8_sf, new_plant_blm_9_sf, new_plant_blm_10_sf) %>%
			mutate(ras_val = 1 - (order/10))
		all_plant_blm_sp <- as(all_plant_blm_sf, "Spatial")
		
		prior_plant_blm_feat_r <- rasterize(all_plant_blm_sp, r_template, field = all_plant_blm_sp$ras_val)
		writeRaster(prior_plant_blm_feat_r, file = "prior_plant_blm_feat_r.tif")
	
		end <- Sys.time()
		el3 <- end - start
		el3
	
	


	
# =============================================================================
#  STEP 2: Set up prioritized inputs from above as new inputs for prioritization of 
#   all inputs together for whole of Sabah
# =============================================================================	
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_blm_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, 
		elev_conn_feat_in, corr_conn_feat_in, acd_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.225
	p2_blm <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(423500) %>%
		add_relative_targets(equal_targets) %>%
		add_boundary_penalties(0.0000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
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
	p3_blm <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(423500) %>%
		add_locked_in_constraints(locked_in_new_blm) %>%
		add_boundary_penalties(0.00000001, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3_blm <- solve(p3_blm)
	s3_blm_sf <- st_as_sf(s3_blm) %>%
		dplyr::filter(solution_1 == 1)
	s3_blm_area_h <- as.integer(sum(st_area(s3_blm_sf)) * 0.0001)
	s3_blm_area_h
	
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	
	
	# ----------------------
	#  Save output.
	scen5_blm_tmp <- as(s3_blm_sf, "Spatial")
	scen5_chunks_blm_sf <- st_as_sf(scen5_blm_tmp) %>%	
		mutate(ras_val = 1)
	scen5_chunks_blm_out <- s3_blm
	save(scen5_chunks_blm_sf, file = "scenario_outputs/scen5/scen5_chunks_blm_sf.Rdata")
	save(scen5_chunks_blm_out, file = "scenario_outputs/scen5/scen5_chunks_blm_out.Rdata")
	
	
	
	
	
	
	
	
	
# =============================================================================
#  Run prioritzation with no BLM constraint with PA locked in
# =============================================================================	
	
	
	
	# ----------------------
	#  Planning unit grids
	load(file = "planning_unit_grids/sa_grid_w_pa.Rdata")
	load(file = "planning_unit_grids/tpa_grid.Rdata")
	
	const_cost_h <- 500
	pu_sf_tmp <- sa_grid_w_pa %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, sep = " ", c("area_h_tmp"), drop = TRUE) 
	pu_sf_tmp$area_h_tmp <- as.numeric(pu_sf_tmp$area_h_tmp)
	pu_sf <- pu_sf_tmp %>%
		mutate(area_h = ifelse(area_h_tmp < 1, 1, area_h_tmp)) %>%
		mutate(const_cost = const_cost_h) %>%
		dplyr::select(-area_h_tmp) %>%
		dplyr::filter(area_h > 100)
	pu_in <- as(pu_sf, "Spatial")

	locked_in <- as(tpa_grid, "Spatial")

	add_area <- as.integer(sum(st_area(tpa_grid)) * 0.0001)
	
	# ----------------------
	#  Input feature stack from inputs generated in step 3.
	all_feat_in <- stack(prior_vert_feat_r, prior_fly_feat_r, prior_plant_feat_r, 
		elev_conn_feat_in, corr_conn_feat_in, acd_feat_in)

	# ----------------------
	#  Set up problem for round 2.
	equal_targets <- 0.3
	p2_pa <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_features_objective(423500 + add_area) %>%
		add_relative_targets(equal_targets) %>%
		add_locked_in_constraints(locked_in) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve round 2.
	s2_pa <- solve(p2_pa)
	s2_pa_sf <- st_as_sf(s2_pa) %>%
		dplyr::filter(solution_1 == 1)
	s2_pa_area_h <- as.integer(sum(st_area(s2_pa_sf)) * 0.0001)
	s2_pa_area_h
	
	# ----------------------
	#  Set up locked in and locked out areas for round 2.
	locked_in_new_pa <- as(s2_pa_sf, 'Spatial')
	
	# ----------------------
	#  Set up problem for round 3.
	p3_pa <- problem(x = pu_in, features = all_feat_in, cost_column = "area_h") %>%
		add_max_utility_objective(423500 + add_area) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve round 3.
	s3_pa <- solve(p3_pa)
	s3_pa_sf <- st_as_sf(s3_pa) %>%
		dplyr::filter(solution_1 == 1)
	s3_pa_area_h <- as.integer(sum(st_area(s3_pa_sf)) * 0.0001)
	s3_pa_area_h
	
	
# =============================================================================
#  STEP 3: Save overall prioritization outputs as sf object, raster and maps.
# =============================================================================	

	# ----------------------
	#  Save output.
	scen5_tmp <- as(s3_sf, "Spatial")
	scen5_chunks_sf <- st_as_sf(scen5_tmp) %>%	
		mutate(ras_val = 1)
	scen5_chunks_out <- s3
	save(scen5_chunks_sf, file = "scenario_outputs/scen5/scen5_chunks_sf.Rdata")
	save(scen5_chunks_out, file = "scenario_outputs/scen5/scen5_chunks_out.Rdata")
	
	
	
	
	# ----------------------
	#  Map 1
	map1 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = scen5_chunks_sf, fill = "dodgerblue4", colour = "dodgerblue4", alpha = 0.7) +
		geom_sf(data = scen5_chunks_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.7) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("No constraints") +
		theme_bw()
	map1
	