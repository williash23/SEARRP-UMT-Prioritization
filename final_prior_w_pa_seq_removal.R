# =============================================================================
#  Set up and solve final prioritization using prioritized inputs from each category
# =============================================================================
	
	# ----------------------
	# Stack all prioritized areas for each category of input with single input removed.
	
	# all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_form_blm, 
		# s_cond_blm, s_acd_blm)
	# all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_form_blm, 
		# s_cond_blm, s_corr_blm)
	# all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_form_blm, 
		# s_acd_blm, s_corr_blm)
	# all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_plant_blm, s_cond_blm, 
		# s_acd_blm, s_corr_blm)
	# all_feat_blm_in <- stack(s_vert_blm, s_fly_blm, s_form_blm, s_cond_blm, 
		# s_acd_blm, s_corr_blm)
	# all_feat_blm_in <- stack(s_vert_blm, s_plant_blm, s_form_blm, s_cond_blm, 
		# s_acd_blm, s_corr_blm)
	# all_feat_blm_in <- stack(s_fly_blm, s_plant_blm, s_form_blm, s_cond_blm, 
		# s_acd_blm, s_corr_blm)

	# ----------------------
	#  Set up problem with all prioritized inputs
	p_all_blm <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_boundary_penalties(0.0000000002, 0.5) %>%
		add_locked_in_constraints(locked_in_pas) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
		
	# ----------------------
	#  Solve.
	s_all <- solve(p_all_blm)
	cellStats(s_all, sum)
	#plot(s_all)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	sol <- s_all
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol

	# ----------------------
	#  Process to remove tiny "speckled" selections
	s_all_sp <- rasterToPolygons(s_all, fun = function(x){x>0}, dissolve = TRUE)
	s_all_sf <- st_as_sf(s_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s_all_sf_drop <- s_all_sf %>%
		drop_crumbs(area_thresh)
	s_all_sf_fill <- s_all_sf_drop %>%
		fill_holes(area_thresh)
	s_all_new_sp <- as(s_all_sf_fill, "Spatial")
	s_all_new <- rasterize(s_all_new_sp, pu_r)
	locked_in_stack_tmp <- stack(locked_in_pas, s_all_new)
	locked_in_sum <- raster::calc(locked_in_stack_tmp, sum, na.rm = TRUE)
	locked_in_sum[locked_in_sum > 0] <- 1

	# ----------------------
	#  Set up locked in areas (previous solution) for next round of selection.
	locked_in_new <- locked_in_sum
	#cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Select additional area to fulfill budget with increased BLM
	p2_all <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_features_objective(prior_area_ha) %>%
		add_relative_targets(0.875) %>%
		add_locked_in_constraints(locked_in_new) %>%
		add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve 
	s2_all <- solve(p2_all)
	#cellStats(s2_all, sum)
	#plot(s2_all)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#   Check coverage of each input.
	
	# prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat_blm_in))

	# sol <- s2_all
	# cov_area_ha <- cellStats(sol, sum)
	# sol[sol > 1] <- 1
	# sol[is.na(sol)] <- 0
	# all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	# for(i in 1:nlayers(all_feat_blm_in_m)){
		# range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		# range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		# cov <- range_cov_sol_sum/range_cov_feat_sum
		# prop_cov_all_feat_sol[i] <- cov
		# }
	
	# prop_cov_all_feat_sol
	
	# ----------------------
	#  Process to remove tiny "speckled" selections again.
	s2_all_sp <- rasterToPolygons(s2_all, fun = function(x){x>0}, dissolve = TRUE)
	s2_all_sf <- st_as_sf(s2_all_sp)
	area_thresh <- units::set_units(1000, ha)
	s2_all_sf_drop <- s2_all_sf %>%
		drop_crumbs(area_thresh)
	s2_all_new_sp <- as(s2_all_sf_drop, "Spatial")
	s2_all_new <- rasterize(s2_all_new_sp, pu_r)
	locked_in_stack_tmp2 <- stack(locked_in_pas, s2_all_new)
	locked_in_sum2 <- raster::calc(locked_in_stack_tmp2, sum, na.rm = TRUE)
	locked_in_sum2[locked_in_sum2 > 0] <- 1
		
	# ----------------------
	#  Set up locked in areas (previous solution) for final round of selection
	locked_in_new <- locked_in_sum2
	#cellStats(locked_in_new, sum)
	#plot(locked_in_new)
	#plot(locked_in_pas, add = T)

	# ----------------------
	#  Select additional area to fulfill final remaining land area. 
	p3_all <- problem(x = pu_r, features = all_feat_blm_in) %>%
		add_max_utility_objective(prior_area_ha) %>%
		#add_max_features_objective(prior_area_ha) %>%
		#add_relative_targets(0.9) %>%
		add_locked_in_constraints(locked_in_new) %>%
		#add_boundary_penalties(0.00000000075, 0.5) %>%
		add_binary_decisions() %>%
		add_gurobi_solver()
	
	# ----------------------
	#  Solve.
	s3_all <- solve(p3_all)
	#cellStats(s3_all, sum)
	#plot(s3_all)
	#plot(locked_in_pas, add = T)
	
	# ----------------------
	#  Rename final output.
	final_sol <- s3_all
	
	
	# ----------------------
	#  Check final coverage of each category of input
	prop_cov_prior_input_sol <- numeric(length = nlayers(all_feat_blm_in))

	sol <- final_sol
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_blm_in_m <- raster::mask(all_feat_blm_in, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat_blm_in_m)){
		range_cov_feat_sum <- cellStats(all_feat_blm_in[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_blm_in_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_prior_input_sol[i] <- cov
		}
	
	prop_cov_prior_input_sol
	
	
	
# =============================================================================
#  Sequentially remove each input and resolve.
# =============================================================================
	
	
	# sol_no_corr <- s3_all
	# sol_no_acd <- s3_all
	# sol_no_cond <- s3_all	
	
	# sol_no_form <- s3_all
	# sol_no_plant <- s3_all
	# sol_no_fly <- s3_all
	# sol_no_vert <- s3_all
	
	
	
	all_sol_stack_w_init_blm <- stack(final_sol, sol_no_vert, sol_no_fly, sol_no_plant, sol_no_form, 
		sol_no_cond, sol_no_acd, sol_no_corr)
	# writeRaster(all_sol_stack_w_init_blm, file = "all_sol_stack_w_init_blm.tif")
		
		
		
		
		
		myColors <- c('darkred',  'dodgerblue4')
		myKey <- list(text=list(lab = c("Selected Units", "PAs")),
			  rectangles=list(col = myColors, border = FALSE),
			  space='bottom', columns=2)
			
		p2 <- levelplot(all_sol_stack_w_init_blm, layout=c(4, 2), margin = FALSE, col.regions=c('transparent', 'dodgerblue4'),
			xlab="", ylab="", scales=list(draw=FALSE), colorkey=FALSE, key = myKey,
			names.attr = c("All Inputs", "No Verts", "No Butterflies", "No Plants", "No Forest Forms",
				"No Elev Conn", "No ACD", "No Corridors"), 
			par.settings=list(layout.heights=list(key.bottom=2,key.top=1)))
		p2 + layer(sp.polygons(main_sabah_tpa_sp, lwd=0.8, col = "transparent", fill = "red", alpha = 0.4))
		
