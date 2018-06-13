

	# ----------------------
	#  Scenario 5
	scen5_sp <- as(scen5_sf, "Spatial")
	scen5_blm_sp <- as(scen5_blm_sf, "Spatial")

	fly_prop_cov_scen5 <- rep(NA, nlayers(fly_feat_in_single))
	plant_prop_cov_scen5 <- rep(NA, nlayers(plant_feat_in_single))
	vert_prop_cov_scen5 <- rep(NA, nlayers(vert_feat_in_single))

	for(i in 1:nlayers(fly_feat_in_single)){
		fly_range_cov_r <- raster::mask(fly_feat_in_single[[i]], scen5_sp)
		fly_range_cov_spp_sum <- cellStats(fly_feat_in_single[[i]], sum, na.rm = TRUE)
		fly_range_cov_sol_sum <- cellStats(fly_range_cov_r, sum, na.rm = TRUE)
		fly_cov <- fly_range_cov_sol_sum/fly_range_cov_spp_sum
		fly_prop_cov_scen5[i] <- fly_cov
		}
	
	for(i in 1:nlayers(plant_feat_in_single)){
		plant_range_cov_r <- raster::mask(plant_feat_in_single[[i]], scen5_sp)
		plant_range_cov_spp_sum <- cellStats(plant_feat_in_single[[i]], sum, na.rm = TRUE)
		plant_range_cov_sol_sum <- cellStats(plant_range_cov_r, sum, na.rm = TRUE)
		plant_cov <- plant_range_cov_sol_sum/plant_range_cov_spp_sum
		plant_prop_cov_scen5[i] <- plant_cov
		}
		
	for(i in 1:nlayers(vert_feat_in_single)){
		vert_range_cov_r <- raster::mask(vert_feat_in_single[[i]], scen5_sp)
		vert_range_cov_spp_sum <- cellStats(vert_feat_in_single[[i]], sum, na.rm = TRUE)
		vert_range_cov_sol_sum <- cellStats(vert_range_cov_r, sum, na.rm = TRUE)
		vert_cov <- vert_range_cov_sol_sum/vert_range_cov_spp_sum
		vert_prop_cov_scen5[i] <- vert_cov
		}
	
	elev_conn_cov_r <- raster::mask(elev_conn_feat_r, scen5_sp)
	elev_conn_cov_sum <- cellStats(elev_conn_feat_r, sum, na.rm = TRUE)
	elev_conn_cov_sol_sum <- cellStats(elev_conn_cov_r, sum, na.rm = TRUE)
	elev_conn_cov <- elev_conn_cov_sol_sum/elev_conn_cov_sum
	elev_conn_prop_cov_scen5 <- elev_conn_cov
	
	corr_cov_r <- raster::mask(corr_feat_r, scen5_sp)
	corr_cov_sum <- cellStats(corr_feat_r, sum, na.rm = TRUE)
	corr_cov_sol_sum <- cellStats(corr_cov_r, sum, na.rm = TRUE)
	corr_cov <- corr_cov_sol_sum/corr_cov_sum
	corr_prop_cov_scen5 <- corr_cov
	
	acd_cov_r <- raster::mask(acd_feat_r, scen5_sp)
	acd_cov_sum <- cellStats(acd_feat_r, sum, na.rm = TRUE)
	acd_cov_sol_sum <- cellStats(acd_cov_r, sum, na.rm = TRUE)
	acd_cov <- acd_cov_sol_sum/acd_cov_sum
	acd_prop_cov_scen5 <- acd_cov
	
	all_cov_scen5 <- c(fly_prop_cov_scen5, plant_prop_cov_scen5, vert_prop_cov_scen5,
		elev_conn_prop_cov_scen5, corr_prop_cov_scen5, acd_prop_cov_scen5)
	cats_scen5 <- c(rep("Butterfly Spp", length(fly_prop_cov_scen5)), 
		rep("Plant Spp", length(plant_prop_cov_scen5)),
		rep("Vertebrate Spp", length(vert_prop_cov_scen5)),
		"Elevational Connectivity", "Corridors", "ACD")
	all_cov_scen5_df <- as.data.frame(cats_scen5) %>%
		cbind(all_cov_scen5)
	names(all_cov_scen5_df) <- c("Proportion", "Category")	
	
	
	fly_prop_cov_scen5_blm <- rep(NA, nlayers(fly_feat_in_single))
	plant_prop_cov_scen5_blm <- rep(NA, nlayers(plant_feat_in_single))
	vert_prop_cov_scen5_blm <- rep(NA, nlayers(vert_feat_in_single))
	
	for(i in 1:nlayers(fly_feat_in_single)){
		fly_range_cov_r <- raster::mask(fly_feat_in_single[[i]], scen5_blm_sp)
		fly_range_cov_spp_sum <- cellStats(fly_feat_in_single[[i]], sum, na.rm = TRUE)
		fly_range_cov_sol_sum <- cellStats(fly_range_cov_r, sum, na.rm = TRUE)
		fly_cov <- fly_range_cov_sol_sum/fly_range_cov_spp_sum
		fly_prop_cov_scen5_blm[i] <- fly_cov
		}
	
	for(i in 1:nlayers(plant_feat_in_single)){
		plant_range_cov_r <- raster::mask(plant_feat_in_single[[i]], scen5_blm_sp)
		plant_range_cov_spp_sum <- cellStats(plant_feat_in_single[[i]], sum, na.rm = TRUE)
		plant_range_cov_sol_sum <- cellStats(plant_range_cov_r, sum, na.rm = TRUE)
		plant_cov <- plant_range_cov_sol_sum/plant_range_cov_spp_sum
		plant_prop_cov_scen5_blm[i] <- plant_cov
		}
		
	for(i in 1:nlayers(vert_feat_in_single)){
		vert_range_cov_r <- raster::mask(vert_feat_in_single[[i]], scen5_blm_sp)
		vert_range_cov_spp_sum <- cellStats(vert_feat_in_single[[i]], sum, na.rm = TRUE)
		vert_range_cov_sol_sum <- cellStats(vert_range_cov_r, sum, na.rm = TRUE)
		vert_cov <- vert_range_cov_sol_sum/vert_range_cov_spp_sum
		vert_prop_cov_scen5_blm[i] <- vert_cov
		}
		
		
	elev_conn_cov_r <- raster::mask(elev_conn_feat_r, scen5_blm_sp)
	elev_conn_cov_sum <- cellStats(elev_conn_feat_r, sum, na.rm = TRUE)
	elev_conn_cov_sol_sum <- cellStats(elev_conn_cov_r, sum, na.rm = TRUE)
	elev_conn_cov <- elev_conn_cov_sol_sum/elev_conn_cov_sum
	elev_conn_prop_cov_scen5_blm <- elev_conn_cov
	
	corr_cov_r <- raster::mask(corr_feat_r, scen5_blm_sp)
	corr_cov_sum <- cellStats(corr_feat_r, sum, na.rm = TRUE)
	corr_cov_sol_sum <- cellStats(corr_cov_r, sum, na.rm = TRUE)
	corr_cov <- corr_cov_sol_sum/corr_cov_sum
	corr_prop_cov_scen5_blm <- corr_cov
	
	acd_cov_r <- raster::mask(acd_feat_r, scen5_blm_sp)
	acd_cov_sum <- cellStats(acd_feat_r, sum, na.rm = TRUE)
	acd_cov_sol_sum <- cellStats(acd_cov_r, sum, na.rm = TRUE)
	acd_cov <- acd_cov_sol_sum/acd_cov_sum
	acd_prop_cov_scen5_blm <- acd_cov
		
	all_cov_scen5_blm <- c(fly_prop_cov_scen5_blm, plant_prop_cov_scen5_blm, 
		vert_prop_cov_scen5_blm, elev_conn_prop_cov_scen5_blm, 
		corr_prop_cov_scen5_blm, acd_prop_cov_scen5_blm)
	cats_scen5_blm <- c(rep("Butterfly Spp", length(fly_prop_cov_scen5_blm)), 
		rep("Plant Spp", length(plant_prop_cov_scen5_blm)),
		rep("Vertebrate Spp", length(vert_prop_cov_scen5_blm)),
		"Elevational Connectivity", "Corridors", "ACD")
	all_cov_scen5_blm_df <- as.data.frame(cats_scen5_blm) %>%
		cbind(all_cov_scen5_blm)
	names(all_cov_scen5_blm_df) <- c("Proportion", "Category")

	
	# ----------------------
	#  Coverage boxplots
	scen5_cov_p <- ggplot(all_cov_scen5_df, aes(x = Category, y = Proportion)) +
		geom_boxplot() +
		coord_flip()
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 5 \nWith BLM") +
		theme_bw()
	scen5_cov_p
		
