###############################################################################
#  Achievement plots
# July 30, 2018
#  Script to generate plots showing coverage of each input feature (i.e., achievement)
#   across each scenario's selected planning units.
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(raster)
library(rgdal)
library(rgeos)



# =============================================================================
#  Individual input features, grouped by categories.
# =============================================================================
	
	# ----------------------
	#  Current scenario prioritization outputs (scenarios 1 - 4)
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen1_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen1_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4_blm_sf.Rdata")

	
# =============================================================================
#  Individual input features, grouped by categories.
# =============================================================================

	# ----------------------
	#   Species ranges and SDMs
	end_fly <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/end_fly.grd")
	end_fly_m <- raster::mask(end_fly, tpa_sp, inverse = TRUE)
	end_fly_m <- raster::mask(end_fly_m, main_sabah_sp)
	end_fly_m <- stack(end_fly_m)
	end_fly_m@filename <- "end_fly"
	
	non_end_fly <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/non_end_fly.grd")
	non_end_fly_m <- raster::mask(non_end_fly, tpa_sp, inverse = TRUE)
	non_end_fly_m <- raster::mask(non_end_fly_m, main_sabah_sp)
	non_end_fly_m <- stack(non_end_fly_m)
	non_end_fly_m@filename <- "non_end_fly"
	
	cr_vert <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/cr_vert.grd")
	cr_vert_m <- raster::mask(cr_vert, tpa_sp, inverse = TRUE)
	cr_vert_m <- raster::mask(cr_vert_m, main_sabah_sp)
	cr_vert_m <- stack(cr_vert_m)
	cr_vert_m@filename <- "cr_vert"
		
	en_vert <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/en_vert.grd")
	en_vert_m <- raster::mask(en_vert, tpa_sp, inverse = TRUE)
	en_vert_m <- raster::mask(en_vert_m, main_sabah_sp)
	en_vert_m <- stack(en_vert_m)
	en_vert_m@filename <- "en_vert"
		
	vu_vert <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/vu_vert.grd")
	vu_vert_m <- raster::mask(vu_vert, tpa_sp, inverse = TRUE)
	vu_vert_m <- raster::mask(vu_vert_m, main_sabah_sp)
	vu_vert_m <- stack(vu_vert_m)
	vu_vert_m@filename <- "vu_vert"
		
	end_plant <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/end_plant.grd")
	end_plant@filename <- "end_plant"
	end_plant_m <- raster::mask(end_plant, tpa_sp, inverse = TRUE)
	end_plant_m <- raster::mask(end_plant_m, main_sabah_sp)
	end_plant_m <- stack(end_plant_m)
	end_plant_m@filename <- "end_plant"
	
	non_end_plant <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/coverage_plots/non_end_plant.grd")
	non_end_plant_m <- raster::mask(non_end_plant, tpa_sp, inverse = TRUE)
	non_end_plant_m <- raster::mask(non_end_plant_m, main_sabah_sp)
	non_end_plant_m <- stack(non_end_plant_m)
	non_end_plant_m@filename <- "non_end_plant"
	
	#  Rare plant localities
	#rare_end_plant_sf <- st_read("C:/Users/saraw/Documents/SEARRP_Analyses/feature_prep/Locked_rare_spp/Endemic/AllEndemicPlants_SquareBuffers.shp")
	#rare_non_end_plant_sf <- st_read("C:/Users/saraw/Documents/SEARRP_Analyses/feature_prep/Locked_rare_spp/Non-endemic/Dipterocarp_500mSquareBufferNonEndemic.shp")
	
	#   Elevational connectivity - Condatis
	load("C:/Users/saraw/Documents/SEARRP_Analyses/feature_prep/stackedresultsreweight.Rdata")
	elev_conn <- condstack2$reweight
	elev_conn_m <- raster::mask(elev_conn, main_sabah_sp)
	elev_conn_m <- stack(elev_conn_m)
	elev_conn_m@filename <- "elev_conn"
	
	#   Corridors
	corr <- stack("C:/Users/saraw/Desktop/feature_inputs/corr_feat_in.grd")
	corr_m <- raster::mask(corr, tpa_sp, inverse = TRUE)
	corr_m <- raster::mask(corr_m, main_sabah_sp)
	corr_m <- stack(corr_m)
	corr_m@filename <- "corr"

	#   ACD
	acd <- stack("C:/Users/saraw/Documents/SEARRP_Analyses/feature_prep/CAO_ACD_30m_unmasked.tif")
	acd@filename <- "acd"
	acd_m <- raster::mask(acd, tpa_sp, inverse = TRUE)
	acd_m <- raster::mask(acd_m, main_sabah_sp)
	acd_m <- stack(acd_m)
	acd_m@filename <- "acd"

	
	
# =============================================================================
#  Acvhievement plots for non-BLM scenarios
# =============================================================================

	# ----------------------
	# List of inputs and scenarios (non-BLM)	
	input_ls <- list(end_fly_m, non_end_fly_m, cr_vert_m, en_vert_m, vu_vert_m, end_plant_m, 
		non_end_plant_m, elev_conn_m, corr_m, acd_m)
	scen_ls <- list(scen1_sf, scen2_sf, scen3_sf, scen4_sf)

	# ----------------------
	# Loop over each scenario
	for(k in 1:length(scen_ls)){
		
		scen_out <- scen_ls[[k]]
		scen_nam <- paste("cov_scen", k, sep = "_")
		plot_nam1 <- paste("point_p_cov_scen", k, sep = "_")
		plot_nam2 <- paste("box_p_cov_scen", k, sep = "_")
		
		scen_sp <- as(st_union(st_buffer(scen_out, 10)), "Spatial")
		
		# ----------------------
		#  Loop over each input -calculate total area of each input feature and 
		#   determine coverage by selected planning units
		for(j in 1:length(input_ls)){
			stack <- input_ls[[j]]
			nam <- paste(stack@filename, "prop_cov", sep = "_")
			
			prop_cov <- rep(NA, nlayers(stack))
				
			for(i in 1:nlayers(stack)){
				range_cov_r <- raster::mask(stack[[i]], scen_sp)
				range_cov_spp_sum <- cellStats(stack[[i]], sum, na.rm = TRUE)
				range_cov_sol_sum <- cellStats(range_cov_r, sum, na.rm = TRUE)
				cov <- range_cov_sol_sum/range_cov_spp_sum
				prop_cov[i] <- cov
				}
			
			assign(nam, prop_cov)
			}
		
		# ----------------------
		#  Combine output of proportion covered across all input layers
		#   into single DF
		all_cov_scen <- c(end_fly_prop_cov, non_end_fly_prop_cov,
			end_plant_prop_cov, non_end_plant_prop_cov, 
			cr_vert_prop_cov, en_vert_prop_cov, vu_vert_prop_cov,
			elev_conn_prop_cov, corr_prop_cov, acd_prop_cov)
		cats_scen <- c(rep("EB", length(end_fly_prop_cov)), 
			rep("NB", length(non_end_fly_prop_cov)), 
			rep("EP", length(end_plant_prop_cov)),
			rep("NP", length(non_end_plant_prop_cov)),
			rep("CV", length(cr_vert_prop_cov)),
			rep("EV", length(en_vert_prop_cov)),
			rep("VV", length(vu_vert_prop_cov)),
			"ELEV", "CORR", "ACD")
		all_cov_scen_df <- as.data.frame(cats_scen) %>%
			cbind(all_cov_scen)
		names(all_cov_scen_df) <- c("Category", "Proportion")	
		all_cov_scen_df$Category <- factor(all_cov_scen_df$Category , levels = c("ACD", 
			"ELEV", "CORR", "EP", "NP", "EB", "NB", "CV", "EV", "VV"))
		
		assign(scen_nam, all_cov_scen_df)
		
		# ----------------------
		# Box plot - for current scenario
		box_cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_boxplot(aes(colour = Category)) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", 
					"Endemic Plants", "Non-Endemic Plants", 
					"Endemic Butterflies", "Non-endemic Butterflies",
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
			
		#assign(plot_nam1, point_cov_p)
		assign(plot_nam2, box_cov_p)
		}
	
	
	
# =============================================================================
#  Acvhievement plots for BLM scenarios
# =============================================================================

	# ----------------------
	# List of inputs and scenarios (BLM)
	input_ls <- list(end_fly_m, non_end_fly_m, cr_vert_m, en_vert_m, vu_vert_m, end_plant_m, 
		non_end_plant_m, elev_conn_m, corr_m, acd_m)
	scen_ls <- list(scen1_blm_sf, scen2_blm_sf, scen3_blm_sf, scen4_blm_sf)
		
	# ----------------------
	# Loop over each scenario
	for(k in 1:length(scen_ls)){
		
		scen_out <- scen_ls[[k]]
		scen_nam <- paste("cov_scen_blm", k, sep = "_")
		plot_nam1 <- paste("point_p_cov_scen_blm", k, sep = "_")
		plot_nam2 <- paste("box_p_cov_scen_blm", k, sep = "_")
		
		scen_sp <- as(st_union(st_buffer(scen_out, 10)), "Spatial")
		
		# ----------------------
		#  Loop over each input -calculate total area of each input feature and 
		#   determine coverage by selected planning units
		for(j in 1:length(input_ls)){
			stack <- input_ls[[j]]
			nam <- paste(stack@filename, "prop_cov", sep = "_")
			
			prop_cov <- rep(NA, nlayers(stack))
				
			for(i in 1:nlayers(stack)){
				range_cov_r <- raster::mask(stack[[i]], scen_sp)
				range_cov_spp_sum <- cellStats(stack[[i]], sum, na.rm = TRUE)
				range_cov_sol_sum <- cellStats(range_cov_r, sum, na.rm = TRUE)
				cov <- range_cov_sol_sum/range_cov_spp_sum
				prop_cov[i] <- cov
				}
			
			assign(nam, prop_cov)
			}
		
		# ----------------------
		#  Combine output of proportion covered across all input layers
		#   into single DF
		all_cov_scen <- c(end_fly_prop_cov, non_end_fly_prop_cov,
			end_plant_prop_cov, non_end_plant_prop_cov, 
			cr_vert_prop_cov, en_vert_prop_cov, vu_vert_prop_cov,
			elev_conn_prop_cov, corr_prop_cov, acd_prop_cov)
		cats_scen <- c(rep("EB", length(end_fly_prop_cov)), 
			rep("NB", length(non_end_fly_prop_cov)), 
			rep("EP", length(end_plant_prop_cov)),
			rep("NP", length(non_end_plant_prop_cov)),
			rep("CV", length(cr_vert_prop_cov)),
			rep("EV", length(en_vert_prop_cov)),
			rep("VV", length(vu_vert_prop_cov)),
			"ELEV", "CORR", "ACD")
		all_cov_scen_df <- as.data.frame(cats_scen) %>%
			cbind(all_cov_scen)
		names(all_cov_scen_df) <- c("Category", "Proportion")	
		all_cov_scen_df$Category <- factor(all_cov_scen_df$Category , levels = c("ACD", 
			"ELEV", "CORR", "EP", "NP", "EB", "NB", "CV", "EV", "VV"))
		
		assign(scen_nam, all_cov_scen_df)
		
		
		# ----------------------
		# Box plot - for current scenario
		box_cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_boxplot(aes(colour = Category)) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", "Endemic Butterflies", "Non-endemic Butterflies",
					"Endemic Plants", "Non-Endemic Plants", 
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
			
		#assign(plot_nam1, point_cov_p)
		assign(plot_nam2, box_cov_p)
		}
	
	
	
# =============================================================================
#  Acvhievement plots for rare plant localities
# =============================================================================
	
	# ----------------------
	#  Bar plots for rare endemic plant localities
	# inputs: rare_end_plant_sf, rare_non_end_plant_sf)
	
	# ----------------------
	#  Endemic rare plant localities input
	input_sf <- rare_end_plant_sf
	rare_spp_cov <- matrix(nrow =  length(scen_ls), ncol = nrow(input_sf))
	
	# ----------------------
	#  Loop over scenarios
	for(k in 1:length(scen_ls)){
		
		scen_out <- scen_ls[[k]]
		scen_sp <- as(st_union(st_buffer(scen_out, 10)), "Spatial")
		scen_sf <- st_as_sf(scen_sp)
		
		# ----------------------
		#  Loop over rows of input - each row is a location with potentially
		#   multiples species represented
		for(i in 1:nrow(input_sf)){
			
			num_spp <- as.numeric(input_sf$Species[i])
			pt_cov_int <- st_intersection(input_sf[i,], scen_sf)
			pt_cov_num <- ifelse(nrow(pt_cov_int) < 1, 0, 1)
			
			pt_spp_cov <- num_spp * pt_cov_num
			rare_spp_cov[k,i] <- pt_spp_cov
			}
		}
	
	rare_end_spp_cov <- rare_spp_cov
	rare_end_spp_cov <- as.data.frame(rare_end_spp_cov)
	cletters <- rep(c("Site"), times= ncol(rare_end_spp_cov))
	cindexes <-seq(1, ncol(rare_end_spp_cov), 1)
	cnames  <- paste(cletters, cindexes, sep="_")
	colnames(rare_end_spp_cov) <- cnames

	rare_end_spp_cov_tmp1 <- rare_end_spp_cov %>%
		dplyr::mutate(Tot_Num_Spp = rowSums(.[1:ncol(rare_end_spp_cov)])) %>%
		dplyr::select(Tot_Num_Spp)
	
	rare_end_spp_cov_tmp2 <- rare_end_spp_cov
	rare_end_spp_cov_tmp2[rare_end_spp_cov_tmp2 > 0] <- 1
	rare_end_spp_cov_tmp2 <- rare_end_spp_cov_tmp2 %>%
		dplyr::mutate(Tot_Num_Sites = rowSums(.[1:ncol(rare_end_spp_cov)])) %>%
		dplyr::select(Tot_Num_Sites) %>%
	
	col_tmp <- as.data.frame(seq(1, length(scen_ls), 1))
	colnames(col_tmp) <- "Scenario"
	
	rare_end_p <- cbind(col_tmp, rare_end_spp_cov_tmp1, rare_end_spp_cov_tmp2)
	rare_end_p$Prop_Sites = rare_end_p$Tot_Num_Sites/ncol(rare_end_spp_cov)
	
	# ----------------------
	#  Bar plot 
	rare_end_cov_p <- ggplot(rare_end_p, aes(as.factor(Scenario), Prop_Sites, 
		fill = as.factor(Scenario))) +
		geom_bar(stat = "identity") +
		scale_fill_brewer(palette = "Dark2", 
			name = "Scenario",
			breaks = c(1,2,3,4,5,6), 
			labels = c("1", "2", "3", "4", "5", "6")) +
		xlab("Scenario") +
		ylab("Proportion of Rare Endemic Sites Protected \n") +
		ylim(0, 1) +
		#geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	rare_end_cov_p
	
	
	# ----------------------
	#  Non-endemic rare plant localities input
	input_sf <- rare_non_end_plant_sf
	rare_spp_cov <- matrix(nrow =  length(scen_ls), ncol = nrow(input_sf))
	
	# ----------------------
	#  Loop over scenarios
	for(k in 1:length(scen_ls)){
		
		scen_out <- scen_ls[[k]]
		scen_sp <- as(st_union(st_buffer(scen_out, 10)), "Spatial")
		scen_sf <- st_as_sf(scen_sp)
		
		# ----------------------
		#  Loop over rows of input - each row is a location with potentially
		#   multiples species represented
		for(i in 1:nrow(input_sf)){
			
			num_spp <- as.numeric(input_sf$Species[i])
			pt_cov_int <- st_intersection(input_sf[i,], scen_sf)
			pt_cov_num <- ifelse(nrow(pt_cov_int) < 1, 0, 1)
			
			pt_spp_cov <- num_spp * pt_cov_num
			rare_spp_cov[k,i] <- pt_spp_cov
			}
		}
	
	rare_non_end_spp_cov <- rare_spp_cov
	rare_non_end_spp_cov <- as.data.frame(rare_non_end_spp_cov)
	cletters <- rep(c("Site"), times= ncol(rare_non_end_spp_cov))
	cindexes <-seq(1, ncol(rare_non_end_spp_cov), 1)
	cnames  <- paste(cletters, cindexes, sep="_")
	colnames(rare_non_end_spp_cov) <- cnames

	rare_non_end_spp_cov_tmp1 <- rare_non_end_spp_cov %>%
		dplyr::mutate(Tot_Num_Spp = rowSums(.[1:ncol(rare_non_end_spp_cov)])) %>%
		dplyr::select(Tot_Num_Spp)
	
	rare_non_end_spp_cov_tmp2 <- rare_non_end_spp_cov
	rare_non_end_spp_cov_tmp2[rare_non_end_spp_cov_tmp2 > 0] <- 1
	rare_non_end_spp_cov_tmp2 <- rare_non_end_spp_cov_tmp2 %>%
		dplyr::mutate(Tot_Num_Sites = rowSums(.[1:ncol(rare_non_end_spp_cov)])) %>%
		dplyr::select(Tot_Num_Sites) 
	
	col_tmp <- as.data.frame(seq(1, length(scen_ls), 1))
	colnames(col_tmp) <- "Scenario"
	
	rare_non_end_p <- cbind(col_tmp, rare_non_end_spp_cov_tmp1, rare_non_end_spp_cov_tmp2)
	rare_non_end_p$Prop_Sites = rare_non_end_p$Tot_Num_Sites/ncol(rare_non_end_spp_cov)
	
	# ----------------------
	#  Bar plot 
	rare_non_end_cov_p <- ggplot(rare_non_end_p, aes(as.factor(Scenario), Prop_Sites, 
		fill = as.factor(Scenario))) +
		geom_bar(stat = "identity") +
		scale_fill_brewer(palette = "Dark2", 
			name = "Scenario",
			breaks = c(1,2,3,4,5,6), 
			labels = c("1", "2", "3", "4", "5", "6")) +
		xlab("Scenario") +
		ylab("Proportion of Rare Non-endemic Sites Protected \n") +
		ylim(0, 1) +
		#geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
		theme_bw()
	rare_non_end_cov_p
	
	
	
# =============================================================================
#  Other plotting ideas...
# =============================================================================	
	
		point_cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_point(aes(colour = Category), shape = 15, size = 5, alpha = 0.6) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", 
					"Endemic Plants", "Non-Endemic Plants", 
					"Endemic Butterflies", "Non-endemic Butterflies",
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
		
	point_cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_point(aes(colour = Category), shape = 15, size = 5, alpha = 0.6) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", "Endemic Butterflies", "Non-endemic Butterflies",
					"Endemic Plants", "Non-Endemic Plants", 
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
	
	
		cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_point(aes(colour = Category), shape = 15, size = 5, alpha = 0.6) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", "Endemic Butterflies", "Non-endemic Butterflies",
					"Endemic Plants", "Non-Endemic Plants", 
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
		cov_p
		
		cov_p <- ggplot(all_cov_scen_df, aes(x = Category, y = Proportion)) +
			geom_boxplot(aes(colour = Category)) +
			scale_colour_manual(name = "Features of Conservation Interest",
				values = c("mediumaquamarine", "goldenrod2", "darkslateblue",
					"darkolivegreen1", "darkolivegreen3", "mediumorchid1", "mediumorchid4", 
					"orangered1", "orangered3", "orangered4"),
				breaks = c("ACD", "ELEV", "CORR", "EP", "NP", "EB", "NB", 
					"CV", "EV", "VV"), 
				labels = c("Aboveground Carbon Density", "Elevational Connectiviy",
					"Corridors", "Endemic Butterflies", "Non-endemic Butterflies",
					"Endemic Plants", "Non-Endemic Plants", 
					"Critically Endangered Vertebrates", "Endangered Vertebrates",
					"Vulnerable Vertebrates")) +
			xlab(" ") +
			ylab("Proportion of Feature Protected \n") +
			ylim(0, 1) +
			geom_hline(yintercept = 0.05, linetype="dashed", color = "grey80") +
			theme_bw()
		cov_p
	
