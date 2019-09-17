
# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(smoothr)
library(units)
library(fasterize)
library(dplyr)
library(SiMRiv)
library(units)
library(doParallel)
library(foreach)



setwd("C:/Users/saraw/Documents/Prioritization/Metapop_Capacity/")

load("sabah_tpa.Rdata")
load("main_sabah_sf.Rdata")

load("conn_prior_patches.Rdata")
conn_prior_r <- raster("conn_prior_r.tif")



# =============================================================================
#  Load base layer input data for prioritization with connectivity layers
# =============================================================================
	
	# ----------------------
	# Raster study area land cover (values from 0 to 9). 
	sa_cov <- conn_prior_r
	
	# ----------------------
	# Land cover patches - Connectivity prioritization output
	patch_sf <- conn_prior_patches %>%
		dplyr::mutate(id = 1:n()) 
	
		

# =============================================================================
#  Resistance matrix.
# =============================================================================
	
	# ----------------------
	#  Generate resistance matrix by reclassify forest cover into a resistance map using values 
	#   specified for each land cover type.
	
	# Reclassification:
	# Non-forest & Plantation = 0 -- > resist = 0.95
	# PA and prioritized area = 1 -- > resist = 0.0001
	# Intact forest = 1 -- > resist = 0.0001
	# Logged forest/regrowth = 2 -- > resist = 0.01

	rc_val_for <- c(0, 0.95,
		1, 0.0001,
		2, 0.01)
		
	rc_mat_for <- matrix(rc_val_for, 
		ncol = 2, 
		byrow = TRUE)
	resist <- raster::reclassify(sa_cov, rc_mat_for)	
	
	
	
# =============================================================================
#  Movement simulation parameters
# =============================================================================
	
	# ----------------------
	#  Specify dispersal max and median distances (m).
	disp_dist_max <- 1000
	disp_dist_med <- 200
	
	# ----------------------
	#  Number of simulations (i.e., dispersers), and steps per dispserser path
	nsim <- 1000
	nstep <- 1000

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	turn_conc <- c(0.75)
	
	# ----------------------
	#  Perceptual range of individual in meters.
	percep_range <- c(res(sa_cov)[1]*4) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- 10

	# ----------------------
	#  Weibull distribution of dispersal success based on max distance.
	disp_dist <- seq(0, disp_dist_max, by = 1) 
	disp_fail_probs <- pweibull(disp_dist, 1, disp_dist_max)
	disp_success_probs <- 1 - disp_fail_probs
	disp_persist <- as.data.frame(cbind(disp_dist, disp_success_probs))
	
	

# =============================================================================
#  Simulate movement from a focal patch 
# =============================================================================
	
	# ----------------------
	#  Initiate cluster for parallel processing.
	cl <- makeCluster(detectCores() - 4)
	registerDoParallel(cl)
	
	start_tm <- Sys.time()
	start_tm
	
	# ----------------------
	#  Generate and save end locations from nsim dispersers from each habitat patch.
	out_locs_conn_prior <- foreach(j = 1:nrow(patch_sf),  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {

		# ----------------------	
		#  Select focal patch.
		focal_patch <- patch_sf %>%
			dplyr::filter(id == j)
		
			if(focal_patch$layer == 0){
		
				end_locs_sf <- NA
				} else {
			
				# ----------------------
				#  Run movement simulations for each individual and store output in a list.
				sim_out_list <- list()
				
				for(w in 1:nsim){
				
					focal_patch_sp <- as(focal_patch, 'Spatial')
					bound <- st_boundary(focal_patch)
					bound_sp <- as(bound, "Spatial")
				
						# ----------------------
						#  Randomly select starting location for individuals, which is constrained to boundary 
						#   of the focal patch.
						init_loc <-  spsample(bound_sp, n = 1, "random")
						init_loc_df <- as.data.frame(init_loc) 
						init_loc_m <- as.matrix(init_loc_df[,1:2])

						# ----------------------
						#  Generate the CRW using specified parameters.
						spp_crw <- SiMRiv::species(state(turn_conc, perceptualRange("cir", percep_range), step_l, "CorrelatedRW"))
						
						# ----------------------
						#  Simulate a single movement path from CRW model
						sim <- SiMRiv::simulate(spp_crw, nstep, resist = resist, coords = c(init_loc_m[1], init_loc_m[2]))

						# ----------------------			
						#  Obtain dispersal covariate values and combine correct successful dispersal probability
						#   for each step's distance.
						sim_move_df <- as.data.frame(sim[,1:2])
						colnames(sim_move_df) =  c("x", "y")
						coordinates(sim_move_df) = ~ x + y
						crs(sim_move_df) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs")
						sim_move_sf_tmp <- st_as_sf(sim_move_df)
						
						first_loc <-  sim_move_sf_tmp %>%
							slice(1)
											
						sim_move_sf_tmp2 <- sim_move_sf_tmp %>%
							mutate(sim_path = w) %>%
							mutate(step_num = 1:n()) %>%
							mutate(disp_dist = round(as.numeric(st_distance(., first_loc)))) %>%
							left_join(disp_persist, by = 'disp_dist') 
						sim_move_sf_tmp2[is.na(sim_move_sf_tmp2)] <- 0
						
						sim_move_df_tmp <- sim_move_sf_tmp2 %>%
							rowwise() %>%
							mutate(cont_path = rbinom(1, 1, disp_success_probs)) %>%
							as.data.frame() 
				
						end_step_tmp <- which(sim_move_df_tmp$cont_path == 0)[1]
						if(is.na(end_step_tmp)){
							end_step <- nstep} else{
								end_step <- end_step_tmp
								}
						
						# ----------------------
						#  Select only the last step for dispersers ending location
						sim_move_df <- sim_move_df_tmp[1:end_step,]
						sim_move_sf <- st_as_sf(sim_move_df, sf_column_name = "geometry") %>%
							dplyr::select(sim_path, step_num) %>%
							dplyr::filter(step_num == max(step_num))
						
						# ----------------------			
						#  Save single simulation ending location 
						sim_out_list[[w]] <- sim_move_sf
						}
					
				# ----------------------			
				#  Bind the ending locations for all simulations.
				end_locs_sf <- do.call(rbind, sim_out_list)
				}
			
			# ----------------------			
			#  Specify return of end_locs_sf object to be combined within parallel processing 
			#   procedure via rbind.
			return(end_locs_sf)
			}
	
	save(out_locs_conn_prior, file = "out_locs_conn_prior_butterflies.Rdata")		

	end_tm_mov <- Sys.time()
	end_tm_mov
	
	
	
# =============================================================================
#  Calculate probability of ending disperal in all other patches and place within matrix
# =============================================================================
	
	# ----------------------
	# Factor to adjust patch area by resistance value of landcover type
	# Reclassification:
	# Non-forest & Plantation = 0 -- > resist = 0.95
	# PA and prioritized area = 1 -- > resist = 0.0001
	# Intact forest = 1 -- > resist = 0.0001
	# Logged forest/regrowth = 2 -- > resist = 0.01

	low_resist <- min(unique(resist)) # 0.0001 (PA/prioritized area)
	hab_corr_fac <- c((low_resist/.95), # correction for non-forest
		(low_resist/0.0001), # correction for PA/prioritized area
		(low_resist/0.01)) # correction for logged forest
	# Patch types: 
	# 0 = Non-forest
	# 1 = PA/Prioritized area/Intact forest
	# 2 = Logged forest
	
	# ----------------------
	#  Determine number of dispsers from each focal patch that ended up in every other patch from 
	#   output "out_locs" from above.
	out_probs_conn_prior <- foreach(k = 1:nrow(patch_sf), .combine = cbind, .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
	
		mat_val_vec <- rep(NA, nrow(patch_sf))
			
		focal_patch <- patch_sf %>%
				dplyr::filter(id == k)
		
			if(focal_patch$layer == 0){
			
				mat_val_vec <- rep(0, nrow(patch_sf))
				} else {
			
					focal_end_locs <- out_locs_conn_prior[[k]] %>%
						st_transform(st_crs(patch_sf))
					focal_patch_area <- as.numeric(st_area(focal_patch) * 0.0001)
					focal_patch_cov_type <- focal_patch$layer
					focal_patch_corr_area <- focal_patch_area * hab_corr_fac[focal_patch_cov_type + 1]
					
					mat_val_vec <- rep(NA, nrow(patch_sf))
			
					for(p in 1:nrow(patch_sf)){
						
						end_patch <- patch_sf %>%
							dplyr::filter(id == p)
					
						end_patch_cov_type <- end_patch$layer
						disp_sucs <- st_intersection(st_buffer(focal_end_locs, 0.001), st_buffer(end_patch, 0.001))
						prob_disp_sucs <- nrow(disp_sucs)/nsim
						
						end_patch_area <- as.numeric(st_area(end_patch)) * 0.0001
						end_patch_corr_area <- end_patch_area * hab_corr_fac[end_patch_cov_type + 1]
						end_patch_corr_area_sqrt <- sqrt(end_patch_corr_area)
						
						mat_val <- prob_disp_sucs * focal_patch_corr_area * end_patch_corr_area_sqrt
						mat_val_vec[p] <- mat_val
						}
					}
			
		
		# ----------------------			
		#  Specify return of prob_vec object to be combined within parallel processing 
		#   procedure via cbind.
		return(mat_val_vec)
		}

	# ----------------------			
	#  Convert out_probs to df and rename columns
	prob_df <- as.data.frame(out_probs_conn_prior)
	features <- c(sprintf("focal_%02d", seq(1, nrow(prob_df))))
	names(prob_df)[1:nrow(prob_df)] <- features
	
	# ----------------------			
	#  Stop parallel processing
	stopCluster(cl)
	
	
	
# =============================================================================
#  Reshape into easy to use matrix with proper value enforcements.
# =============================================================================
	
	# ----------------------			
	#  Enforce that symmetrical cells in matrix are average of two cell
	#  combinations
	prob_mat_ave <- matrix(NA, nrow = nrow(prob_df), ncol = ncol(prob_df))
	
	for(m in 1:ncol(prob_mat_ave)) {
		for(n in 1:nrow(prob_mat_ave)) {
			val1 <- prob_df[m,n]
			val2 <- prob_df[n,m]
			ave_val <- (val1 + val2)/2
			prob_mat_ave[m,n] <- ave_val
			prob_mat_ave[n,m] <- ave_val
			}
		}
		
	# ----------------------	
	#  Calculate patch area (corrected for habitat type) of each patch and place in 
	#   diaganol of matrix.
	area_df <- patch_sf %>%
		dplyr::mutate(area = as.numeric(st_area(.)) * 0.0001) %>%
		dplyr::mutate(log10area = log10(as.numeric(st_area(.)) * 0.0001)) %>%
		dplyr::mutate(hab_corr_area = area * hab_corr_fac[layer + 1])
	corr_area_vec <- area_df$hab_corr_area
	prob_mat_ave_area <- prob_mat_ave
	
	for(q in 1:nrow(prob_mat_ave_area)) {
		prob_mat_ave_area[q,q] <- corr_area_vec[q]
		}
	
	
	
# =============================================================================
#  Save final output.
# =============================================================================
	
	# ----------------------			
	#  Save.
	prob_mat_conn_prior <- prob_mat_ave_area
	save(prob_mat_conn_prior, file  = "prob_mat_conn_prior_butterflies.Rdata")
	
	end_tm <- Sys.time()
	
	
	rm(prob_mat_ave_area, prob_df, prob_mat_ave, patch_sf, sa_cov)
	
# =============================================================================	
###############################################################################

