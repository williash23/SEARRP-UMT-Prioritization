

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
library(smoothr)


setwd("C:/Users/saraw/Documents/Prioritization/")


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
	#  Forest formations inputs - 	raster 13, 18 and 19 (in order) in the stack have no 
	#   area outside of existing PAs so remove for this analysis.
	forest_form_feat_in <- stack("feature_inputs/forest_form_feat_in_w_pa_range.grd")
	forest_form_feat_in <- forest_form_feat_in[[-13]]
	forest_form_feat_in <- forest_form_feat_in[[-18]]
	forest_form_feat_in <- forest_form_feat_in[[-19]]
	forest_form_feat_in[is.na(forest_form_feat_in)] <- 0
	
		
	# ----------------------
	#  Set up problem for connectivity and carbon input layers
	elev_conn_feat_in <- raster("feature_inputs/elev_conn_feat_in_w_pa_range.grd")
	corr_conn_feat_in <- raster("feature_inputs/new_corr_conn_feat_in.grd")
	acd_feat_in <- raster("feature_inputs/acd_feat_in_w_pa_range.grd")
	
	
	# ----------------------
	#  Planning unit grids
	pu_r <- raster("mainland_sabah_planning_units.tif")
	tpa_r <- raster("mainland_sabah_tpa.tif") 
	area_ha <- cellStats(pu_r, sum)
	

	# ----------------------
	#   Mask species ranges by the solution.
	#sol_sf <- st_read("scenario_outputs/12_1_18/all_inputs_prioritization_w_pa.shp")
	#sol_sf <- st_read("scenario_outputs/12_1_18/no_connectivity_inputs_prioritization_w_pa.shp")
	#cov_area_ha <- as.numeric(st_area(sol_sf)) * 0.0001
	#sol <- as(sol_sf, "Spatial")
	#all_feat_m <- raster::mask(all_feat, sol)  #, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	
	
	
	
	sol_tmp <- stack("scenario_outputs/12_1_18/all_sol_stack_181201.grd")
	sol <- sol_tmp[[1]]
	#sol <- raster("scenario_outputs/12_1_18/no_conn_sol.tif")
	cov_area_ha <- cellStats(sol, sum)
	sol_no_pa <- raster::mask(sol, tpa_r, inverse = TRUE)
	sol_area_ha <- cellStats(sol_no_pa, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	

	all_feat <- stack(vert_feat_in, fly_feat_in, plant_feat_in, forest_form_feat_in, 
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)
	all_feat_m <- raster::mask(all_feat, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)
	
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat))

	for(i in 1:nlayers(all_feat)){
		range_cov_feat_sum <- cellStats(all_feat[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	vert_cov <- prop_cov_all_feat_sol[1:81]
	fly_cov <- prop_cov_all_feat_sol[82:158]
	plant_cov <- prop_cov_all_feat_sol[159:555]
	form_cov <- prop_cov_all_feat_sol[556:574]
	elev_conn_cov <- prop_cov_all_feat_sol[575]
	acd_cov <- prop_cov_all_feat_sol[576]
	corr_conn_cov <- prop_cov_all_feat_sol[577]
	
	# Remove rare plants
	plant_cov1 <- prop_cov_all_feat_sol[159:206]
	plant_cov2 <- prop_cov_all_feat_sol[512:529]
	plant_cov <- c(plant_cov1, plant_cov2)
	
	summary(vert_cov)
	summary(fly_cov)
	summary(plant_cov)
	summary(form_cov)
	summary(elev_conn_cov)
	summary(corr_conn_cov)
	summary(acd_cov)
	
	length(vert_cov[vert_cov > .3199999])
	length(fly_cov[fly_cov > .3199999])
	length(plant_cov[plant_cov > .3199999])
	length(form_cov[form_cov > .3199999])
	
	
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, rep(1, length(plant_cov)),
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))

	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		#geom_point(data = cat_cov, aes(x = Category, y = Proportion, colour = Category)) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p
		
		
		
		
		
	### Removing existing TPA network
	
	sol_no_pa <- raster::mask(sol, tpa_r, inverse = TRUE)
	sol_area_ha <- cellStats(sol_no_pa, sum)
	sol <- sol_no_pa
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	

	all_feat <- stack(vert_feat_in, fly_feat_in, plant_feat_in, forest_form_feat_in, 
		elev_conn_feat_in, acd_feat_in, corr_conn_feat_in)
	all_feat <- raster::mask(all_feat, tpa_r, inverse = TRUE) # Remove feature range that is within TPA
	all_feat_m <- raster::mask(all_feat, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)
	
	prop_cov_all_feat_sol <- numeric(length = nlayers(all_feat))

	for(i in 1:nlayers(all_feat)){
		range_cov_feat_sum <- cellStats(all_feat[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_sol[i] <- cov
		}
	
	prop_cov_all_feat_sol
	
	vert_cov <- prop_cov_all_feat_sol[1:81]
	fly_cov <- prop_cov_all_feat_sol[82:158]
	plant_cov <- prop_cov_all_feat_sol[159:555]
	form_cov <- prop_cov_all_feat_sol[556:574]
	elev_conn_cov <- prop_cov_all_feat_sol[575]
	acd_cov <- prop_cov_all_feat_sol[576]
	corr_conn_cov <- prop_cov_all_feat_sol[577]
	
	# Remove rare plants
	plant_cov1 <- prop_cov_all_feat_sol[159:206]
	plant_cov2 <- prop_cov_all_feat_sol[512:529]
	plant_cov <- c(plant_cov1, plant_cov2)
	
	summary(vert_cov)
	summary(fly_cov)
	summary(plant_cov)
	summary(form_cov)
	summary(elev_conn_cov)
	summary(corr_conn_cov)
	summary(acd_cov)
	
	length(vert_cov[vert_cov > .0499999])
	length(fly_cov[fly_cov > .0499999])
	length(plant_cov[plant_cov > .0499999])
	length(form_cov[form_cov > .0499999])
	
	
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, rep(1, length(plant_cov)),
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))

	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		geom_point(data = all_cov_df, aes(x = Category, y = Proportion, colour = Category)) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	
	
	
	
	#cat_cov <- as.data.frame(cbind(c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"), 
		#prop_cov_prior_input_sol))
	#colnames(cat_cov) <- c("Category", "Proportion")
	#cat_cov$Proportion <- as.numeric(as.character(cat_cov$Proportion))
	

	prop_cov_all_feat_tpa <- numeric(length = nlayers(all_feat))

	# ----------------------
	#   Mask species ranges by the solution.
	sol <- locked_in_pas
	cov_area_ha <- cellStats(sol, sum)
	sol[sol > 1] <- 1
	sol[is.na(sol)] <- 0
	all_feat_m <- raster::mask(all_feat, sol, maskvalue = 1, updatevalue = 0, inverse = TRUE)

	for(i in 1:nlayers(all_feat)){
		range_cov_feat_sum <- cellStats(all_feat[[i]], sum, na.rm = TRUE)
		range_cov_sol_sum <- cellStats(all_feat_m[[i]], sum, na.rm = TRUE)
		cov <- range_cov_sol_sum/range_cov_feat_sum
		prop_cov_all_feat_tpa[i] <- cov
		}
	
	add_cov <- prop_cov_all_feat_sol - prop_cov_all_feat_tpa
	
	vert_cov <- add_cov[1:81]
	fly_cov <- add_cov[82:158]
	plant_cov <- add_cov[159:555]
	form_cov <- add_cov[556:574]
	elev_conn_cov <- add_cov[575]
	acd_cov <- add_cov[576]
	corr_conn_cov <- add_cov[577]
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, plant_rep_weight_w_pa_range,
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
		
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))

		
	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		#geom_point(aes(colour = Category, shape = Weight)) +
		#scale_shape_manual(name = "Weight", values = c(1, 2, 8, 19, 15, 5)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p

	
	
	# Remove rare plants
	plant_cov1 <- add_cov[159:206]
	plant_cov2 <- add_cov[512:529]
	plant_cov <- c(plant_cov1, plant_cov2)
	
	
	# ----------------------
	#  Combine output of proportion covered across all input layers
	#   into single DF
	all_cov <- c(vert_cov, fly_cov, plant_cov, form_cov, elev_conn_cov,
		acd_cov, corr_conn_cov)
	cats <- c(rep("VERT", length(vert_cov)), 
		rep("FLY", length(fly_cov)), 
		rep("PLANT", length(plant_cov)),
		rep("FORM", length(form_cov)),
		rep("ELEV", length(elev_conn_cov)),
		rep("ACD", length(acd_cov)),
		rep("CORR", length(corr_conn_cov)))
	weights <- c(vert_rep_weight_w_pa_range, fly_rep_weight_w_pa_range, rep(1, length(plant_cov)),
		rep(1, nlayers(forest_form_feat_in)), 1, 1, 1)
	all_cov_df <- as.data.frame(cats) %>%
		cbind(all_cov, weights)
	names(all_cov_df) <- c("Category", "Proportion", "Weight")	
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 4] <- "CR"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 3] <- "EN"
	all_cov_df$Weight[all_cov_df$Category == "VERT" & all_cov_df$Weight == 2] <- "VU"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "FLY" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 2] <- "END"
	all_cov_df$Weight[all_cov_df$Category == "PLANT" & all_cov_df$Weight == 1] <- "NON-END"
	all_cov_df$Weight[all_cov_df$Weight == 1] <- "NONE"
	all_cov_df <- all_cov_df %>%
		dplyr::filter(!is.na(Proportion)) 
	
		
	all_cov_df$Category <- factor(all_cov_df$Category , levels = c("VERT", 
		"FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"))
	
	
	# ----------------------
	# Box plot - for current scenario
	cov_p <- ggplot(all_cov_df, aes(x = Category, y = Proportion)) +
		#geom_point(aes(colour = Category, shape = Weight)) +
		#scale_shape_manual(name = "Weight", values = c(1, 2, 8, 19, 15, 5)) +
		geom_boxplot(aes(colour = Category)) +
		scale_colour_manual(name = "Features of Conservation Interest",
			values = c("mediumorchid1", "goldenrod2", "darkgreen",
				"darkolivegreen3", "darkred", "orangered3", "dodgerblue4"),
			breaks = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR"),
			labels = c("VERT", "FLY", "PLANT", "FORM", "ELEV", "ACD", "CORR")) +
		ylab("Proportion of Feature Protected \n") +
		ylim(0, 1) +
		geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	cov_p
	
	
	
	
	