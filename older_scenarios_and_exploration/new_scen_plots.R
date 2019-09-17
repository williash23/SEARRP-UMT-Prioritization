###############################################################################
#  Prioritization scenario plots
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
library(ggplot2)


setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Boundaries
	load(file ="study_area_boundaries_and_pas/sabah_border.Rdata")
	load(file ="study_area_boundaries_and_pas/main_sabah_sf.Rdata")
	load(file ="study_area_boundaries_and_pas/kali_border.Rdata")
	load(file ="study_area_boundaries_and_pas/sarawak_border.Rdata")
	load(file = "study_area_boundaries_and_pas/sabah_tpa.Rdata")
	load(file = "planning_unit_grids/sa_grid_sf.Rdata")

	# ----------------------
	#  Scenario outputs
	scen1_blm_r <- raster("scenario_outputs/scen1_blm_r.tif")
	scen1_r <- raster("scenario_outputs/scen1_r.tif")
	scen1_w_pa_blm_r <- raster("scenario_outputs/scen1_w_pa_blm_r.tif")
	scen1_w_pa_r <- raster("scenario_outputs/scen1_w_pa_r.tif")

	
	
	
# =============================================================================
#  Prep data for plotting
# =============================================================================		

	# ----------------------
	#  Prep PA and planning unit grids
	sabah_tpa$CLASS <- as.character("1")
	tpa <- sabah_tpa %>%
		dplyr::select(IN = CLASS)
	sa_grid_sf$IN <- as.character("1")
	pu_in <- as(sa_grid_sf, "Spatial")
	
	# ----------------------
	#  Convert scenario outputs to sf objects
	scen1_blm_sp <- raster::rasterToPolygons(scen1_blm_r)
	scen1_blm_sf <- st_as_sf(scen1_blm_sp)
	scen1_sp <- raster::rasterToPolygons(scen1_r)
	scen1_sf <- st_as_sf(scen1_sp)
	scen1_w_pa_blm_sp <- raster::rasterToPolygons(scen1_w_pa_blm_r)
	scen1_w_pa_blm_sf <- st_as_sf(scen1_w_pa_blm_sp)
	scen1_w_pa_sp <- raster::rasterToPolygons(scen1_w_pa_r)
	scen1_w_pa_sf <- st_as_sf(scen1_w_pa_sp)

	
	
# =============================================================================
#  Plots
# =============================================================================		

	# ----------------------
	#  No PA inclusion, no BLM
	sel <- scen1_sf %>%
		mutate(IN = ifelse(scen1_r > 0.5, "2", "0")) %>%
		dplyr::filter(IN == "2") %>%
		dplyr::select(IN)
	p <- rbind(tpa, sel)

	map1 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = p,  aes(colour = IN, fill = IN), alpha = 0.8) +
		scale_fill_manual(name = "", values = c("goldenrod3", "sienna3"), 
			labels = c("Existing TPA", "Selected Planning Unit")) +
		scale_colour_manual(guide = FALSE, values = c("goldenrod3", "sienna3")) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ylab("Latitude") +
		xlab("Longitude") +
		theme_bw() +
		theme(legend.direction = "horizontal", legend.position = "bottom") +
		ggtitle("410,000 ha prioritized - existing PAs not included - no BLM")
	map1
	
	# ----------------------
	#  No PA inclusion, BLM
	sel <- scen1_blm_sf %>%
		mutate(IN = ifelse(scen1_blm_r > 0.5, "2", "0")) %>%
		dplyr::filter(IN == "2") %>%
		dplyr::select(IN)
	p <- rbind(tpa, sel)

	map2 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = p,  aes(colour = IN, fill = IN), alpha = 0.8) +
		scale_fill_manual(name = "", values = c("goldenrod3", "sienna3"), 
			labels = c("Existing TPA", "Selected Planning Unit")) +
		scale_colour_manual(guide = FALSE, values = c("goldenrod3", "sienna3")) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ylab("Latitude") +
		xlab("Longitude") +
		theme_bw() +
		theme(legend.direction = "horizontal", legend.position = "bottom") +
		ggtitle("410,000 ha prioritized - existing PAs not included - with BLM")
	map2
	
	# ----------------------
	#  PA inclusion, no BLM
	sel <- scen1_w_pa_sf %>%
		mutate(IN = ifelse(scen1_w_pa_r > 0.5, "2", "0")) %>%
		dplyr::filter(IN == "2") %>%
		dplyr::select(IN)
	p <- rbind(tpa, sel)

	map3 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = p,  aes(colour = IN, fill = IN), alpha = 0.8) +
		scale_fill_manual(name = "", values = c("goldenrod3", "sienna3"), 
			labels = c("Existing TPA", "Selected Planning Unit")) +
		scale_colour_manual(guide = FALSE, values = c("goldenrod3", "sienna3")) +
		geom_sf(data = sabah_tpa, colour = "goldenrod3", fill = "goldenrod3") +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ylab("Latitude") +
		xlab("Longitude") +
		theme_bw() +
		theme(legend.direction = "horizontal", legend.position = "bottom") +
		ggtitle("2,281,300 ha prioritized - existing PAs included - no BLM")
	map3
	
	# ----------------------
	#  PA inclusion, BLM
	sel <- scen1_w_pa_blm_sf %>%
		mutate(IN = ifelse(scen1_w_pa_blm_r > 0.5, "2", "0")) %>%
		dplyr::filter(IN == "2") %>%
		dplyr::select(IN)
	p <- rbind(tpa, sel)

	map4 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = p,  aes(colour = IN, fill = IN), alpha = 0.8) +
		scale_fill_manual(name = "", values = c("goldenrod3", "sienna3"), 
			labels = c("Existing TPA", "Selected Planning Unit")) +
		scale_colour_manual(guide = FALSE, values = c("goldenrod3", "sienna3")) +
		geom_sf(data = sabah_tpa, colour = "goldenrod3", fill = "goldenrod3") +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ylab("Latitude") +
		xlab("Longitude") +
		theme_bw() +
		theme(legend.direction = "horizontal", legend.position = "bottom") +
		ggtitle("2,281,300 ha prioritized - existing PAs included - with BLM")
	map4
	
	
	
	