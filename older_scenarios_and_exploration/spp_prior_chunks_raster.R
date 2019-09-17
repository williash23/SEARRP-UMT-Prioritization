###############################################################################
#  Prioritization - Prepare species ranges inputs in chunks
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
range01 <- function(x){(x-r_min)/(r_max-r_min)}


# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Load species ranges (1 layer per species).
	vert_feat_in <- stack("feature_inputs/vert_all.grd")
	fly_feat_in <- stack("feature_inputs/fly_all.grd")	
	plant_feat_in <- stack("feature_inputs/plant_all.grd")
	
	# ----------------------
	#  Load species weighting for single layer species stacks.
	load("feature_inputs/vert_rep_weight.Rdata")
	load("feature_inputs/fly_rep_weight.Rdata")
	load("feature_inputs/plant_rep_weight.Rdata")
	
	# ----------------------
	#  Forest formations inputs - 	raster 13, 18 and 19 (in order) in the stack have no 
	#   area outside of existing PAs so remove for this analysis.
	forest_form_feat_in <- stack("feature_inputs/forest_form_feat_in.grd")
	forest_form_feat_in <- forest_form_feat_in[[-13]]
	forest_form_feat_in <- forest_form_feat_in[[-18]]
	forest_form_feat_in <- forest_form_feat_in[[-19]]
	forest_form_feat_in[is.na(forest_form_feat_in)] <- 0
	
	# ----------------------
	#  Create raster following template of study area 
	temp <- fly_feat_in[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	

# =============================================================================
#  Set up planning units WITHOUT existing TPAs
# =============================================================================

	# ----------------------
	#  Planning unit grids
	pu_r_tmp <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 
	pu_r <- raster::mask(pu_r_tmp, tpa_r, inverse = TRUE)
	
	area_ha <- cellStats(pu_r, sum)
	
	
# =============================================================================
#  Run prioritzation with no BLM constraint
# =============================================================================	
	

# =============================================================================
#  STEP 1: Set up and solve problem for species ranges for whole of Sabah
# =============================================================================
	
		
		# ----------------------
		#  Vertebrates
		chunks <- seq(1, 15, by = 1)
	
		adj_vert_rep_weight <- vert_rep_weight
		adj_vert_rep_weight[adj_vert_rep_weight == 4] <- 0.1
		adj_vert_rep_weight[adj_vert_rep_weight == 3] <- 0.075
		adj_vert_rep_weight[adj_vert_rep_weight == 2] <- 0.05
	
		v_rel_tar <- chunks[1]*adj_vert_rep_weight
		
		p_vert_1 <- problem(x = pu_r, features = vert_feat_in) %>%
			add_max_features_objective(area_ha*(1/10)) %>%
			add_relative_targets(v_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		s_vert_1 <- solve(p_vert_1)
		s_vert_1[s_vert_1 == 1] <- 1
		locked_in <- s_vert_1
		
		for(i in 2:length(chunks)){
			
			v_rel_tar <- chunks[i]*adj_vert_rep_weight
			v_rel_tar[v_rel_tar > 1] <- 1
			
			p_vert_tmp <- problem(x = pu_r, features = vert_feat_in) %>%
				add_max_features_objective(area_ha*(i/10)) %>%
				add_relative_targets(v_rel_tar) %>%				
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()

			s_vert_tmp <- solve(p_vert_tmp)
			s_vert_tmp[s_vert_tmp == 1] <- i
			locked_in <- s_vert_tmp
			locked_in[locked_in > 1] <- 1
			
			nam <- paste("s_vert", i, sep = "_")
			assign(nam, s_vert_tmp)		
			}
		
		
		i = i +1
		area_left <- area_ha - (cellStats(locked_in, sum) * 100)
		
		p_vert_tmp <- problem(x = pu_r, features = vert_feat_in) %>%
			add_max_utility_objective(area_ha) %>%
			add_locked_in_constraints(locked_in) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		s_vert_tmp <- solve(p_vert_tmp)
		s_vert_tmp[s_vert_tmp == 1] <- i
		locked_in[locked_in > 1] <- 1
		locked_in <- s_vert_tmp
	
		nam <- paste("s_vert", i, sep = "_")
		assign(nam, s_vert_tmp)		
		
		vert_stack <- stack(s_vert_1, s_vert_2, s_vert_3, s_vert_4, s_vert_5, s_vert_6, s_vert_7, s_vert_8,
			 s_vert_9, s_vert_10, s_vert_11)
		vert_stack[vert_stack == 0] <- NA
		vert_stack_min <- stackApply(vert_stack, indices = rep(1, nlayers(vert_stack)), fun = min)
		prior_vert_feat_r_tmp <- (nlayers(vert_stack) + 1) - vert_stack_min
		r_min <- cellStats(prior_vert_feat_r_tmp, min)
		r_max <- cellStats(prior_vert_feat_r_tmp, max)
		prior_vert_feat_r <- raster::calc(prior_vert_feat_r_tmp, fun = range01)
		
		writeRaster(prior_vert_feat_r, file = "chunks/prior_vert_feat_r.grd")
	
		
		
		# ----------------------
		#  Butterflies
		chunks <- seq(0.1, 1, by = 0.1)
		
		rev_fly_rep_weight <- fly_rep_weight
		rev_fly_rep_weight[rev_fly_rep_weight == 2] <- 10
		rev_fly_rep_weight[rev_fly_rep_weight == 1] <- 2
		rev_fly_rep_weight[rev_fly_rep_weight == 10] <- 1
	
		f_rel_tar <- chunks[1]*(1/rev_fly_rep_weight)
		
		p_fly_1 <- problem(x = pu_r, features = fly_feat_in) %>%
			add_max_features_objective(area_ha*(1/10)) %>%
			add_relative_targets(f_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver() 

		s_fly_1 <- solve(p_fly_1)
		s_fly_1[s_fly_1 == 1] <- 1
		locked_in <- s_fly_1
		
		for(i in 2:length(chunks)){
		
			f_rel_tar <- chunks[i]*(1/rev_fly_rep_weight)

			p_fly_tmp <- problem(x = pu_r, features = fly_feat_in) %>%
				add_max_features_objective(area_ha*(i/10)) %>%
				add_relative_targets(f_rel_tar) %>% 
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()

			s_fly_tmp <- solve(p_fly_tmp)
			s_fly_tmp[s_fly_tmp == 1] <- i
			locked_in <- s_fly_tmp
			locked_in[locked_in > 1] <- 1
			
			nam <- paste("s_fly", i, sep = "_")
			assign(nam, s_fly_tmp)		
			}
		
		i = i +1
		area_left <- area_ha - (cellStats(locked_in, sum) * 100)
		
		p_fly_tmp <- problem(x = pu_r, features = fly_feat_in) %>%
			add_max_utility_objective(area_ha) %>%
			add_locked_in_constraints(locked_in) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		s_fly_tmp <- solve(p_fly_tmp)
		s_fly_tmp[s_fly_tmp == 1] <- i
		locked_in[locked_in > 1] <- 1
		locked_in <- s_fly_tmp
	
		nam <- paste("s_fly", i, sep = "_")
		assign(nam, s_fly_tmp)		
		
		fly_stack <- stack(s_fly_1, s_fly_2, s_fly_3, s_fly_4, s_fly_5, s_fly_6, s_fly_7, s_fly_8,
			 s_fly_9, s_fly_10, s_fly_11)
		fly_stack[fly_stack == 0] <- NA
		fly_stack_min <- stackApply(fly_stack, indices = rep(1, nlayers(fly_stack)), fun = min)
		prior_fly_feat_r_tmp <- (nlayers(fly_stack) + 1) - fly_stack_min
		r_min <- cellStats(prior_fly_feat_r_tmp, min)
		r_max <- cellStats(prior_fly_feat_r_tmp, max)
		prior_fly_feat_r <- raster::calc(prior_fly_feat_r_tmp, fun = range01)
		
		writeRaster(prior_fly_feat_r, file = "chunks/prior_fly_feat_r.grd")
		
		
		
		# ----------------------
		#  Plants
		chunks <- seq(0.1, 2, by = 0.1)
		
		rev_plant_rep_weight <- plant_rep_weight
		rev_plant_rep_weight[rev_plant_rep_weight == 2] <- 10
		rev_plant_rep_weight[rev_plant_rep_weight == 1] <- 2
		rev_plant_rep_weight[rev_plant_rep_weight == 10] <- 1
	
		p_rel_tar <- chunks[1]*(1/rev_plant_rep_weight)
		
		p_plant_1 <- problem(x = pu_r, features = plant_feat_in) %>%
			add_max_features_objective(area_ha*(1/10)) %>%
			add_relative_targets(p_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		s_plant_1 <- solve(p_plant_1)
		s_plant_1[s_plant_1 == 1] <- 1
		locked_in <- s_plant_1
		
		for(i in 2:7){
		
			p_rel_tar <- chunks[i]*(1/rev_plant_rep_weight)
			p_rel_tar[p_rel_tar > 1] <- 1
			
			p_plant_tmp <- problem(x = pu_r, features = plant_feat_in) %>%
				add_max_features_objective(area_ha*(i/10)) %>%
				add_relative_targets(p_rel_tar) %>%
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			s_plant_tmp <- solve(p_plant_tmp)
			s_plant_tmp[s_plant_tmp == 1] <- i
			locked_in <- s_plant_tmp
			locked_in[locked_in > 1] <- 1
			
			nam <- paste("s_plant", i, sep = "_")
			assign(nam, s_plant_tmp)		
			}
	
	
		for(i in 7:length(chunks)){
		
			p_rel_tar <- chunks[i]*(1/rev_plant_rep_weight)
			p_rel_tar[p_rel_tar > 1] <- 1
			
			p_plant_tmp <- problem(x = pu_r, features = plant_feat_in) %>%
				add_max_cover_objective(area_ha*(i/10)) %>%
				#add_feature_weights(p_rel_tar) %>%
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			s_plant_tmp <- solve(p_plant_tmp)
			s_plant_tmp[s_plant_tmp == 1] <- i
			locked_in <- s_plant_tmp
			locked_in[locked_in > 1] <- 1
			
			nam <- paste("s_plant", i, sep = "_")
			assign(nam, s_plant_tmp)		
			}
			
		i = i +1
		area_left <- area_ha - (cellStats(locked_in, sum) * 100)
		
		p_plant_tmp <- problem(x = pu_r, features = plant_feat_in) %>%
			add_max_utility_objective(area_ha) %>%
			add_locked_in_constraints(locked_in) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		s_plant_tmp <- solve(p_plant_tmp)
		s_plant_tmp[s_plant_tmp == 1] <- i
		locked_in[locked_in > 1] <- 1
		locked_in <- s_plant_tmp
	
		nam <- paste("s_plant", i, sep = "_")
		assign(nam, s_plant_tmp)		
		
		plant_stack <- stack(s_plant_1, s_plant_2, s_plant_3, s_plant_4, s_plant_5, s_plant_6, 
			s_plant_7, s_plant_8, s_plant_9, s_plant_10, s_plant_11)
		plant_stack[plant_stack == 0] <- NA
		plant_stack_min <- stackApply(plant_stack, indices = rep(1, nlayers(plant_stack)), fun = min)
		prior_plant_feat_r_tmp <- (nlayers(plant_stack) + 1) - plant_stack_min
		r_min <- cellStats(prior_plant_feat_r_tmp, min)
		r_max <- cellStats(prior_plant_feat_r_tmp, max)
		prior_plant_feat_r <- raster::calc(prior_plant_feat_r_tmp, fun = range01)
		
		writeRaster(prior_plant_feat_r, file = "chunks/prior_plant_feat_r.grd")
	
	
	
	
		# ----------------------
		#  Forest formations
		chunks <- seq(0.1, 1, by = 0.1)
		
		f_rel_tar <- chunks[1]
		
		p_form_1 <- problem(x = pu_r, features = forest_form_feat_in) %>%
			add_max_features_objective(area_ha*(1/10)) %>%
			add_relative_targets(f_rel_tar) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
	
		s_form_1 <- solve(p_form_1)
		s_form_1[s_form_1 == 1] <- 1
		locked_in <- s_form_1
		
		for(i in 2:length(chunks)){
		
			f_rel_tar <- chunks[i]

			p_form_tmp <- problem(x = pu_r, features = forest_form_feat_in) %>%
				add_max_features_objective(area_ha*(i/10)) %>%
				add_relative_targets(f_rel_tar) %>%
				add_locked_in_constraints(locked_in) %>%
				add_binary_decisions() %>%
				add_gurobi_solver()
		
			s_form_tmp <- solve(p_form_tmp)
			s_form_tmp[s_form_tmp == 1] <- i
			locked_in <- s_form_tmp
			locked_in[locked_in > 1] <- 1
			
			nam <- paste("s_form", i, sep = "_")
			assign(nam, s_form_tmp)		
			}
	
		i = i +1
		area_left <- area_ha - (cellStats(locked_in, sum) * 100)
		
		p_form_tmp <- problem(x = pu_r, features = forest_form_feat_in) %>%
			add_max_utility_objective(area_ha) %>%
			add_locked_in_constraints(locked_in) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		s_form_tmp <- solve(p_form_tmp)
		s_form_tmp[s_form_tmp == 1] <- i
		locked_in[locked_in > 1] <- 1
		locked_in <- s_form_tmp
	
		nam <- paste("s_form", i, sep = "_")
		assign(nam, s_form_tmp)		
		
		form_stack <- stack(s_form_1, s_form_2, s_form_3, s_form_4, s_form_5, s_form_6, 
			s_form_7, s_form_8, s_form_9, s_form_10, s_form_11)
		form_stack[form_stack == 0] <- NA
		form_stack_min <- stackApply(form_stack, indices = rep(1, nlayers(form_stack)), fun = min)
		prior_form_feat_r_tmp <- (nlayers(form_stack) + 1) - form_stack_min
		r_min <- cellStats(prior_form_feat_r_tmp, min)
		r_max <- cellStats(prior_form_feat_r_tmp, max)
		prior_form_feat_r <- raster::calc(prior_form_feat_r_tmp, fun = range01)

		writeRaster(prior_form_feat_r, file = "chunks/prior_form_feat_r.grd")
	