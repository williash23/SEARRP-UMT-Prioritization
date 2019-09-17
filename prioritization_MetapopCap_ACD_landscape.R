
# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(smoothr)
library(units)
library(rmapshaper)
 library(fasterize)
 library(dplyr)
 library(SiMRiv)
library(units)
library(smoothr)
library(doParallel)
library(foreach)



 
setwd("C:/Users/saraw/Documents/")
 
acd <- raster("Spatial_data/raw_spat_data/CAO_ACD_30m_unmasked.tif")

no_conn_sol <- st_read("Prioritization/scenario_outputs/12_1_18/no_connectivity_inputs_prioritization_w_pa.shp") %>%
	mutate(value = 100)

 sol <- st_read("Prioritization/scenario_outputs/12_1_18/all_inputs_prioritization_w_pa.shp") %>%
	mutate(value = 100)
	
 load("Prioritization/study_area_boundaries_and_pas/sabah_tpa.Rdata")
 load("Prioritization/study_area_boundaries_and_pas/main_sabah_sf.Rdata")
 
	
acd_agg <- raster::aggregate(acd, fact = 30)

	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(acd_agg), ncol(acd_agg))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(acd_agg)
	projection(r_template) <- CRS("+proj=omerc +lat_0=4 +lonc=115 +alpha=53.31582047222222 +k=0.99984 +x_0=590476.87 +y_0=442857.65 +gamma=53.31582047222222 +ellps=evrstSS +towgs84=-679,669,-48,0,0,0,0 +units=m +no_defs")

	
nc_sol_r <- fasterize(no_conn_sol, r_template, field = "value")
sol_r <- fasterize(sol, r_template, field = "value")
 
 
 
# =============================================================================
#  Reclassify to get elevational bands within intact forest and logged forest/regrowth.
# =============================================================================
	
	# ----------------------
	#  Reclassify based on ACD grouping (higher ACD leads to lower resistance
	rc_val <- c(0, 50.999999, 7,
		51, 100.999999, 6,
		101, 150.99999, 5,
		151, 200.999999, 4,
		201, 250.999999, 3,
		251, 300.9999999, 2,
		301, 350, 1)
		
	rc_mat <- matrix(rc_val, 
		ncol = 3, 
		byrow = TRUE)
	acd_resist <- raster::reclassify(acd_agg, rc_mat, include.lowest = TRUE)
	
	
nc_sol_resist_tmp <- raster::stack(nc_sol_r, acd_resist)
nc_sol_resist <- sum(nc_sol_resist_tmp, na.rm = TRUE)
nc_sol_resist[nc_sol_resist == 0]  <- NA
nc_sol_resist[nc_sol_resist > 7] <- 0

sol_resist_tmp <- raster::stack(sol_r, acd_resist)
sol_resist <- sum(sol_resist_tmp, na.rm = TRUE)
sol_resist[sol_resist == 0]  <- NA
sol_resist[sol_resist > 7] <- 0


	
# =============================================================================
#  Convert to polygons and simplfy to generate habitat patches.
# =============================================================================
	
	# ----------------------
	#  Convert to polygons.
	nc_sol_p_tmp <- spex::polygonize(nc_sol_resist) 
	cov_types <- sort(unique(nc_sol_p_tmp$layer))
	
	area_thresh <- set_units(2, km^2)
	
	for(i in cov_types){
		sf_tmp <- nc_sol_p_tmp %>%
			dplyr::filter(layer == i)
		if(nrow(sf_tmp) > 1){
			sf_u <- st_union(sf_tmp) %>%
				fill_holes(area_thresh) 
			sf_poly <- st_cast(sf_u, "POLYGON") 
			sf_poly_sp <- as(sf_poly, "Spatial")
			sf_out <- st_as_sf(sf_poly_sp) %>%
				dplyr::mutate(value = i)
			}else{
				sf_out <- sf_tmp %>%
					dplyr::mutate(value = i) %>%
					dplyr::select(value)
				}
		nam <- paste("landcov_nc", i, sep = "_")
		assign(nam, sf_out)
		}
	
	nc_sol_p <- rbind(landcov_nc_0, landcov_nc_3, landcov_nc_4, landcov_nc_5, 
		landcov_nc_6, landcov_nc_7) %>%
		st_intersection(st_transform(main_sabah_sf, st_crs(nc_sol_p_tmp)))
	nc_sol_p <- st_as_sf(as(nc_sol_p, "Spatial"))
	# 2531

	nc_sol_p_r <- fasterize(nc_sol_p, r_template, field = "value")
	nc_sol_p_r[nc_sol_p_r == 0] <- 100
	# 100 = PA
	# 1 = 301-350 ACD  ### Not present in raster
	# 2 = 251-300 ACD ### Not present in raster
	# 3 = 201-250 ACD
	# 4 = 151-200 ACD
	# 5 = 101-150 ACD
	# 6 = 51-100 ACD
	# 7 = 0-50 ACD
	
	

	sol_p_tmp <- spex::polygonize(sol_resist) 
	cov_types <- sort(unique(sol_p_tmp$layer))
	
	area_thresh <- set_units(2, km^2)
	
	for(i in cov_types){
		sf_tmp <- sol_p_tmp %>%
			dplyr::filter(layer == i)
		if(nrow(sf_tmp) > 1){
			sf_u <- st_union(sf_tmp) %>%
				fill_holes(area_thresh) 
			sf_poly <- st_cast(sf_u, "POLYGON") 
			sf_poly_sp <- as(sf_poly, "Spatial")
			sf_out <- st_as_sf(sf_poly_sp) %>%
				dplyr::mutate(value = i)
			}else{
				sf_out <- sf_tmp %>%
					dplyr::mutate(value = i) %>%
					dplyr::select(value)
				}
		nam <- paste("landcov", i, sep = "_")
		assign(nam, sf_out)
		}
	
	sol_p <- rbind(landcov_0, landcov_2, landcov_3, landcov_4, landcov_5, landcov_6, landcov_7) %>%
		st_intersection(st_transform(main_sabah_sf, st_crs(nc_sol_p_tmp)))
	sol_p <- st_as_sf(as(sol_p, "Spatial"))
	# 2497 patches
	
	sol_p_r <- fasterize(sol_p, r_template, field = "value")
	# 0 = PA
	# 1 = 301-350 ACD  ### Not present in raster
	# 2 = 251-300 ACD 
	# 3 = 201-250 ACD
	# 4 = 151-200 ACD
	# 5 = 101-150 ACD
	# 6 = 51-100 ACD
	# 7 = 0-50 ACD

	
	
	rc_val <- c(0, 10,
		2, 20,
		3, 30,
		4, 40,
		5, 50,
		6, 60,
		7, 70)
	rc_mat <- matrix(rc_val, 
		ncol =2, 
		byrow = TRUE)
	sol_p_r_rc <- raster::reclassify(sol_p_r, rc_mat, include.lowest = TRUE)
	# 10 = PA
	# 20 = 251-300 ACD 
	# 30 = 201-250 ACD
	# 40 = 151-200 ACD
	# 50 = 101-150 ACD
	# 60 = 51-100 ACD
	# 70 = 0-50 ACD


	both_sol_r <- stack(nc_sol_p_r, sol_p_r_rc)
	both_sol_r_sum <- sum(both_sol_r)
	
	sol_sf <- spex::polygonize(both_sol_r_sum) 
	cov_types <- sort(unique(sol_sf$layer))
	
	area_thresh <- set_units(5, km^2)
	
	for(i in cov_types){
		sf_tmp <- sol_sf %>%
			dplyr::filter(layer == i)
		if(nrow(sf_tmp) > 1){
			sf_u <- st_union(sf_tmp) %>%
				fill_holes(area_thresh) 
			sf_poly <- st_cast(sf_u, "POLYGON") 
			sf_poly_sp <- as(sf_poly, "Spatial")
			sf_out <- st_as_sf(sf_poly_sp) %>%
				dplyr::mutate(value = i)
			}else{
				sf_out <- sf_tmp %>%
					dplyr::mutate(value = i) %>%
					dplyr::select(value)
				}
		nam <- paste("landcov_both", i, sep = "_")
		assign(nam, sf_out)
		}
		
	both_sol_all_p <- rbind(landcov_both_13, landcov_both_14, landcov_both_15, landcov_both_16, landcov_both_17,
		landcov_both_33, landcov_both_43, landcov_both_44, landcov_both_45, landcov_both_55, landcov_both_56,
		landcov_both_66, landcov_both_67, landcov_both_76, landcov_both_77, landcov_both_110, landcov_both_120,
		landcov_both_130, landcov_both_140, landcov_both_150, landcov_both_160, landcov_both_170)
	
	both_sol_all_p <- st_as_sf(as(both_sol_all_p, "Spatial"))
	# 3023
	
	#13 = Conn Sol PA + No Conn Sol 3
	#14 = Conn Sol PA + No Conn Sol 4
	#15 = Conn Sol PA + No Conn Sol 5
	#16 = Conn Sol PA + No Conn Sol 6
	#17 = Conn Sol PA + No Conn Sol 7
	#33 = Conn Sol 3 + No Conn Sol 3
	#43 = Conn Sol 4 + No Conn Sol 3
	#44 = Conn Sol 4+ No Conn Sol 4
	#45 = Conn Sol 4 + No Conn Sol 5
	#55 = Conn Sol 5 + No Conn Sol 5
	#56 = Conn Sol 5 + No Conn Sol 6
	#66 = Conn Sol 6 + No Conn Sol 6
	#67 = Conn Sol 6 + No Conn Sol 7
	#76 = Conn Sol 7 + No Conn Sol 6
	#77 = Conn Sol 7 + No Conn Sol 7
	#110 =  Conn Sol PA + No Conn Sol PA
	#120 = Conn Sol 2 + No Conn Sol PA 
	#130 = Conn Sol 3 + No Conn Sol PA 
	#140 = Conn Sol 4 + No Conn Sol PA 
	#150 = Conn Sol 5 + No Conn Sol PA 
	#160 = Conn Sol 6 + No Conn Sol PA 
	#170 = Conn Sol 72 + No Conn Sol PA 
	
	
	
	# Conn Prior Landscape Patches:
	conn_prior_patches <- both_sol_all_p 
	conn_prior_patches$value[conn_prior_patches$value == 13] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 14] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 15] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 16] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 17] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 33] <- 3
	conn_prior_patches$value[conn_prior_patches$value == 43] <- 4
	conn_prior_patches$value[conn_prior_patches$value == 44] <- 4
	conn_prior_patches$value[conn_prior_patches$value == 45] <- 4
	conn_prior_patches$value[conn_prior_patches$value == 55] <- 5
	conn_prior_patches$value[conn_prior_patches$value == 56] <- 5
	conn_prior_patches$value[conn_prior_patches$value == 66] <- 6
	conn_prior_patches$value[conn_prior_patches$value == 67] <- 6
	conn_prior_patches$value[conn_prior_patches$value == 76] <- 7
	conn_prior_patches$value[conn_prior_patches$value == 77] <- 7
	conn_prior_patches$value[conn_prior_patches$value == 110] <- 0
	conn_prior_patches$value[conn_prior_patches$value == 120] <- 2
	conn_prior_patches$value[conn_prior_patches$value == 130] <- 3
	conn_prior_patches$value[conn_prior_patches$value == 140] <- 4
	conn_prior_patches$value[conn_prior_patches$value == 150] <- 5
	conn_prior_patches$value[conn_prior_patches$value == 160] <- 6
	conn_prior_patches$value[conn_prior_patches$value == 170] <- 7
	conn_prior_r <- fasterize(conn_prior_patches, r_template, field = "value")
	
	
	# No Conn Prior Landscape Patches:
	no_conn_prior_patches <- both_sol_all_p 
	no_conn_prior_patches$value[no_conn_prior_patches$value == 13] <- 3
	no_conn_prior_patches$value[no_conn_prior_patches$value == 14] <- 4
	no_conn_prior_patches$value[no_conn_prior_patches$value == 15] <- 5
	no_conn_prior_patches$value[no_conn_prior_patches$value == 16] <- 6
	no_conn_prior_patches$value[no_conn_prior_patches$value == 17] <- 7
	no_conn_prior_patches$value[no_conn_prior_patches$value == 33] <- 3
	no_conn_prior_patches$value[no_conn_prior_patches$value == 43] <- 3
	no_conn_prior_patches$value[no_conn_prior_patches$value == 44] <- 4
	no_conn_prior_patches$value[no_conn_prior_patches$value == 45] <- 5
	no_conn_prior_patches$value[no_conn_prior_patches$value == 55] <- 5
	no_conn_prior_patches$value[no_conn_prior_patches$value == 56] <- 6
	no_conn_prior_patches$value[no_conn_prior_patches$value == 66] <- 6
	no_conn_prior_patches$value[no_conn_prior_patches$value == 67] <- 7
	no_conn_prior_patches$value[no_conn_prior_patches$value == 76] <- 6
	no_conn_prior_patches$value[no_conn_prior_patches$value == 77] <- 7
	no_conn_prior_patches$value[no_conn_prior_patches$value == 110] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 120] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 130] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 140] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 150] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 160] <- 0
	no_conn_prior_patches$value[no_conn_prior_patches$value == 170] <- 0
	
	no_conn_prior_r <- fasterize(no_conn_prior_patches, r_template, field = "value")
	
	


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
	rc_val_for <- c(0, 0.01,
		1, 0.14,
		2, 0.28,
		3, 0.42, 
		4, 0.56,
		5, 0.70,
		6, 0.84,
		7, 0.98)
	rc_mat_for <- matrix(rc_val_for, 
		ncol = 2, 
		byrow = TRUE)
	resist <- raster::reclassify(sa_cov, rc_mat_for)
	
	
	
# =============================================================================
#  Simulate movement from a focal patch.
# =============================================================================
	
	# ----------------------
	#  Specify dispersal max and median distances (m).
	disp_dist_max <- 152206 ## Rusa unicolor
	disp_dist_med <- 20127 ## Rusa unicolor
	
	# ----------------------
	#  Number of simulations (i.e., dispersers), and steps per dispserser path
	nsim <- 1000
	nstep <- 1000

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	turn_conc <- c(0.98)
	
	# ----------------------
	#  Perceptual range of individual in meters.
	percep_range <- c(res(sa_cov)[1]*4) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- c(res(sa_cov)[1]*0.5)

	# ----------------------
	#  Weibull distribution of dispersal success based on max distance.
	disp_dist <- seq(0, disp_dist_max, by = 1) 
	disp_fail_probs <- pweibull(disp_dist, 1, disp_dist_max)
	disp_success_probs <- 1 - disp_fail_probs
	disp_persist <- as.data.frame(cbind(disp_dist, disp_success_probs))
	
	
	
	
	
	# ----------------------
	#  Initiate cluster for parallel processing.
	cl <- makeCluster(detectCores()-3)
	registerDoParallel(cl)
	
	# ----------------------
	#  Generate and save end locations from nsim dispersers from each habitat patch.
	out_locs_conn_prior <- foreach(j = 1:nrow(patch_sf),  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {

		# ----------------------	
		#  Select focal patch.
		focal_patch <- patch_sf %>%
			dplyr::filter(id == j)
		
			if(focal_patch$value == 0){
		
				end_locs_sf <-  NA} else {
			
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
						crs(sim_move_df) <- CRS("+proj=omerc +lat_0=4 +lonc=115 +alpha=53.31582047222222 +k=0.99984 +x_0=590476.87 +y_0=442857.65 +gamma=53.31582047222222 +ellps=evrstSS +towgs84=-679,669,-48,0,0,0,0 +units=m +no_defs")
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
	
	#save(out_locs_conn_prior, file = "Prioritization/Metapop_Capacity/out_locs_conn_prior.Rdata")		

	
	
# =============================================================================
#  Calculate probability of ending disperal in all other patches and place within matrix
# =============================================================================
	
	# ----------------------
	# Factor to adjust patch area by resistance value of landcover type
	low_resist <- min(rc_mat_for[,2])
	hab_corr_fac <- low_resist/rc_mat_for[,2]
	
	# ----------------------
	#  Determine number of dispsers from each focal patch that ended up in every other patch from 
	#   output "out_locs" from above.
	out_probs_conn_prior <- foreach(k = 1:nrow(patch_sf), .combine = cbind, .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
	
		mat_val_vec <- rep(NA, nrow(patch_sf))
			
		focal_patch <- patch_sf %>%
				dplyr::filter(id == k)
			
		focal_end_locs <- out_locs_conn_prior[[k]] %>%
			st_transform(st_crs(sol_p))
		focal_patch_area <- as.numeric(st_area(focal_patch) * 0.0001)
		focal_patch_cov_type <- focal_patch$value
		focal_patch_corr_area <- focal_patch_area * hab_corr_fac[focal_patch_cov_type + 1]
		
		mat_val_vec <- rep(NA, nrow(patch_sf))

		for(p in 1:nrow(patch_sf)){
			
			end_patch <- patch_sf %>%
				dplyr::filter(id == p)
		
			end_patch_cov_type <- end_patch$value
			disp_sucs <- st_intersection(st_buffer(focal_end_locs, 0.001), st_buffer(end_patch, 0.001))
			prob_disp_sucs <- nrow(disp_sucs)/nsim
			
			end_patch_area <- as.numeric(st_area(end_patch)) * 0.0001
			end_patch_corr_area <- end_patch_area * hab_corr_fac[end_patch_cov_type + 1]
			end_patch_corr_area_sqrt <- sqrt(end_patch_corr_area)
			
			mat_val <- prob_disp_sucs * focal_patch_corr_area * end_patch_corr_area_sqrt
			mat_val_vec[p] <- mat_val
			}
		
		
		# ----------------------			
		#  Specify return of prob_vec object to be combined within parallel processing 
		#   procedure via cbind.
		return(mat_val_vec)
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
		dplyr::mutate(hab_corr_area = area * hab_corr_fac[value + 1])
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
	save(prob_mat_conn_prior, file  = "prob_mat_con_prior.Rdata")
	
	rm(prob_mat_ave_area, prob_df, prob_mat_ave, patch_sf, sa_cov)
	
# =============================================================================	
###############################################################################










# =============================================================================
#  Load base layer input data for prioritization with connectivity layers
# =============================================================================
	
	# ----------------------
	# Raster study area land cover (values from 0 to 9). 
	sa_cov <- no_conn_prior_r
	
	# ----------------------
	# Land cover patches - Connectivity prioritization output
	patch_sf <- no_conn_prior_patches %>%
		dplyr::mutate(id = 1:n()) 
	
		

# =============================================================================
#  Resistance matrix.
# =============================================================================
	
	# ----------------------
	#  Generate resistance matrix by reclassify forest cover into a resistance map using values 
	#   specified for each land cover type.
	rc_val_for <- c(0, 0.01,
		1, 0.14,
		2, 0.28,
		3, 0.42, 
		4, 0.56,
		5, 0.70,
		6, 0.84,
		7, 0.98)
	rc_mat_for <- matrix(rc_val_for, 
		ncol = 2, 
		byrow = TRUE)
	resist <- raster::reclassify(sa_cov, rc_mat_for)
	
	
	
# =============================================================================
#  Simulate movement from a focal patch.
# =============================================================================
	
	# ----------------------
	#  Specify dispersal max and median distances (m).
	disp_dist_max <- 152206 ## Rusa unicolor
	disp_dist_med <- 20127 ## Rusa unicolor
	
	# ----------------------
	#  Number of simulations (i.e., dispersers), and steps per dispserser path
	nsim <- 1000
	nstep <- 1000

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	turn_conc <- c(0.98)
	
	# ----------------------
	#  Perceptual range of individual in meters.
	percep_range <- c(res(sa_cov)[1]*4) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- c(res(sa_cov)[1]*0.5)

	# ----------------------
	#  Weibull distribution of dispersal success based on max distance.
	disp_dist <- seq(0, disp_dist_max, by = 1) 
	disp_fail_probs <- pweibull(disp_dist, 1, disp_dist_max)
	disp_success_probs <- 1 - disp_fail_probs
	disp_persist <- as.data.frame(cbind(disp_dist, disp_success_probs))
	
	
	
	
	
	# ----------------------
	#  Initiate cluster for parallel processing.
	cl <- makeCluster(detectCores() - 2)
	registerDoParallel(cl)
	
	# ----------------------
	#  Generate and save end locations from nsim dispersers from each habitat patch.
	out_locs_no_conn_prior <- foreach(j = 1:nrow(patch_sf),  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {

		# ----------------------	
		#  Select focal patch.
		focal_patch <- patch_sf %>%
			dplyr::filter(id == j)
		
			if(focal_patch$value == 0){
		
				end_locs_sf <-  NA} else {
			
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
						crs(sim_move_df) <- CRS("+proj=omerc +lat_0=4 +lonc=115 +alpha=53.31582047222222 +k=0.99984 +x_0=590476.87 +y_0=442857.65 +gamma=53.31582047222222 +ellps=evrstSS +towgs84=-679,669,-48,0,0,0,0 +units=m +no_defs")
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
	
	#save(out_locs_no_conn_prior, file = "Prioritization/Metapop_Capacity/out_locs_no_conn_prior.Rdata")		

	
	
# =============================================================================
#  Calculate probability of ending disperal in all other patches and place within matrix
# =============================================================================
	
	# ----------------------
	# Factor to adjust patch area by resistance value of landcover type
	low_resist <- min(rc_mat_for[,2])
	hab_corr_fac <- low_resist/rc_mat_for[,2]
	
	# ----------------------
	#  Determine number of dispsers from each focal patch that ended up in every other patch from 
	#   output "out_locs" from above.
	out_probs_no_conn_prior <- foreach(k = 1:nrow(patch_sf), .combine = cbind, .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
	
		mat_val_vec <- rep(NA, nrow(patch_sf))
			
		focal_patch <- patch_sf %>%
				dplyr::filter(id == k)
			
		focal_end_locs <- out_locs_no_conn_prior[[k]] %>%
			st_transform(st_crs(sol_p))
		focal_patch_area <- as.numeric(st_area(focal_patch) * 0.0001)
		focal_patch_cov_type <- focal_patch$value
		focal_patch_corr_area <- focal_patch_area * hab_corr_fac[focal_patch_cov_type + 1]
		
		mat_val_vec <- rep(NA, nrow(patch_sf))

		for(p in 1:nrow(patch_sf)){
			
			end_patch <- patch_sf %>%
				dplyr::filter(id == p)
		
			end_patch_cov_type <- end_patch$value
			disp_sucs <- st_intersection(st_buffer(focal_end_locs, 0.001), st_buffer(end_patch, 0.001))
			prob_disp_sucs <- nrow(disp_sucs)/nsim
			
			end_patch_area <- as.numeric(st_area(end_patch)) * 0.0001
			end_patch_corr_area <- end_patch_area * hab_corr_fac[end_patch_cov_type + 1]
			end_patch_corr_area_sqrt <- sqrt(end_patch_corr_area)
			
			mat_val <- prob_disp_sucs * focal_patch_corr_area * end_patch_corr_area_sqrt
			mat_val_vec[p] <- mat_val
			}
		
		
		# ----------------------			
		#  Specify return of prob_vec object to be combined within parallel processing 
		#   procedure via cbind.
		return(mat_val_vec)
		}

	# ----------------------			
	#  Convert out_probs to df and rename columns
	prob_df <- as.data.frame(out_probs_no_conn_prior)
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
		dplyr::mutate(hab_corr_area = area * hab_corr_fac[value + 1])
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
	prob_mat_no_conn_prior <- prob_mat_ave_area
	save(prob_mat_no_conn_prior, file  = "prob_mat_no_conn_prior.Rdata")
	
	
	rm(prob_mat_ave_area, prob_df, prob_mat_ave, patch_sf, sa_cov)
	
# =============================================================================	
###############################################################################









	
	