###############################################################################
#  Reclassify Canopy Tree Height.
#  January 29, 2018
#  Script to reclassify cells identified as forest and non-forest from CAO canopy tree height (CTH).
#  Sara Williams
###############################################################################



# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(rasterVis)
library(RColorBrewer)



# =============================================================================
#  !!!!!! Value to be adjusted !!!!!! 
#   Set value for Canopy Tree Height that is considered non-forest.
# =============================================================================

	# ----------------------
	#  Elevation oil palm cutoff (m)
	elev_cutoff <- 54
	# Median palm oil cultivation from Brodie et al. 2015
	
	# ----------------------
	#  Canopy tree height forest cutoff value(s) (m).
	#forest_cutoff <- 10	
	forest_cutoff_h_elev <- 10	
	forest_cutoff_l_elev <- 20
	
	# ----------------------
	#  Aggregation factor. Current resolution is 30m x 30m.
	#agg_fact <- 30


# =============================================================================
#  Load input data.
# =============================================================================
	
	# ----------------------
	#  Canopy Tree Height (CTH) map from Asner et al.
	cth <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/raw_spat_data/CAO_TCH_30m_unmasked.tif")

	# ----------------------
	#  Forest reserves
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/for_res_rt_sf.Rdata")

	# ----------------------
	# Global Multi-resolution Terrain Elevation raster
	elev <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/elev_250m_m.grd")

	
	
# =============================================================================
#  Aggregate CTH to larger cell size
# =============================================================================	

	# ----------------------
	#  Aggregate based on factor selected and using 'modal' function.
	#   To me, aggregation seems to make a big change in what is considered forest.....
	#cth_agg <- aggregate(cth, fact = agg_fact, fun = mean)
	
	

# =============================================================================
#  Resample elevation to match CTH forest cover.
# =============================================================================

	# ----------------------
	#  Resample
	elev_rs <- resample(elev, cth, "bilinear")
	#elev_agg_rs <- resample(elev, cth_agg, "bilinear")
	
	rm(elev)
	
# =============================================================================
#  Reclassify elevation to identify low vs high elevation areas (low elevation can identify oil
#   palm plantation areas).
# =============================================================================

	# ----------------------
	#  Reclassification values for elevation - reset LOW areas to 0
	#   0 = below elev_cutoff
	#   1 = above elev_cutoff
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_elev_h <- c(-Inf, elev_cutoff , 0,
		elev_cutoff  + 0.0000001, Inf, 1)
	rc_mat_elev_h <- matrix(rc_val_elev_h, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	elev_rc_h <- raster::reclassify(elev_rs, rc_mat_elev_h, include.lowest = TRUE)
	#elev_agg_rc <- raster::reclassify(elev_agg_rs, rc_mat_elev, include.lowest = TRUE)
	
	# ----------------------
	#  Reclassification values for elevation - reset HIGH areas to 0
	#   1 = below elev_cutoff
	#   0 = above elev_cutoff
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_elev_l <- c(-Inf, elev_cutoff , 1,
		elev_cutoff  + 0.0000001, Inf, 0)
	rc_mat_elev_l <- matrix(rc_val_elev_l, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	elev_rc_l <- raster::reclassify(elev_rs, rc_mat_elev_l, include.lowest = TRUE)
	#elev_agg_rc <- raster::reclassify(elev_agg_rs, rc_mat_elev, include.lowest = TRUE)

	
	
	
# =============================================================================
#  Combine elevation and CTH to determine variable cutoffs for identifying forest by canopy height.
# =============================================================================
	
	# ----------------------
	#  Overlay CTH and reclassified elevation (both low and high).
	cth_elev_h <- overlay(cth, elev_rc_h, fun=function(x,y){return(x*y)})
	cth_elev_l <- overlay(cth, elev_rc_l, fun=function(x,y){return(x*y)})

	
	
	
# =============================================================================
#  Reclassify canopy tree height raster layer so that values below cut-off (decided above) are 
#   assigned as non-forest (0) and values above are assigned as forest (1).
# =============================================================================

	# ---------------------
	#  Reclassification values for CTH at LOW elevations
	#   0 = non-forest
	#   1 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_for_l <- c(0, forest_cutoff_l_elev, 0,
		forest_cutoff_l_elev + 0.0000001, Inf, 1)
	rc_mat_for_l <- matrix(rc_val_for_l, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	cth_elev_l_rc <- raster::reclassify(cth_elev_l, rc_mat_for_l, include.lowest = TRUE)
	#cth_agg_rc <- raster::reclassify(cth_agg, rc_mat, include.lowest = TRUE)
	
	# ----------------------
	#  Reclassification values for CTH at HIGH elevations
	#   0 = non-forest
	#   1 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_for_h <- c(0, forest_cutoff_h_elev, 0,
		forest_cutoff_h_elev + 0.0000001, Inf, 1)
	rc_mat_for_h <- matrix(rc_val_for_h, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	cth_elev_h_rc <- raster::reclassify(cth_elev_h, rc_mat_for_h, include.lowest = TRUE)
	#cth_agg_rc <- raster::reclassify(cth_agg, rc_mat, include.lowest = TRUE)
	
	cth_by_elev <- overlay(cth_elev_h_rc, cth_elev_l_rc, fun=function(x,y){return(x+y)})
	
	
	
	
	
# =============================================================================
# Remove forest class reserve that holds mangrove forest.
# =============================================================================	

	# ----------------------	
	#  Select only Class V (Mangroves)
	for_res_rt_sf$NOV_2016 <- as.factor(for_res_rt_sf$NOV_2016)
	mang <- for_res_rt_sf %>% 
		dplyr::filter(NOV_2016 == 5)
	mang_sp <- as(mang, "Spatial")
	
	# ----------------------	
	#  Mask out mangroves from reclassified forest raster.
	#cth_rc_no_mang <- raster::mask(cth_rc, mang_sp, updatevalue = 0, inverse = TRUE)
	cth_agg_rc_no_mang <- raster::mask(cth_agg_rc, mang_sp, updatevalue = 0, inverse = TRUE)
	
	# ----------------------
	#  Save output as polygon and raster
	# writeRaster(cth_agg_rc_no_mang, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/cth_agg_rc_no_mang.grd")

	
	
# =============================================================================
#  Polygonize to use with planning unit sp objects.
# =============================================================================	

	# ----------------------	
	# Polygonize
	cth_p <- polygonizer(cth_agg_rc_no_mang)
	cth_sf <- st_as_sf(cth_p)
	cth_sf_for <- cth_sf %>%
		dplyr::filter(DN == 1)

	# ----------------------
	#  Save output as polygon 
	# save(cth_sf_for, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/cth_sf_for.Rdata")
	


# =============================================================================
#  Plots
# =============================================================================

	# ----------------------
	#  Load border data for plotting
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sabah_pa_sp.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sabah_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sarawak_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_kali_d.Rdata")

	# ----------------------
	#  Set color palette themes
	r_theme <- rasterTheme(rev(brewer.pal(2, "Accent")))
	r_theme$panel.background$col = 'grey80' 
	r_theme$regions$col = c("#ccb388", "#008b45")

	# ----------------------
	#  Plot all forest as assigned by Gaveau forest cover layer
	gav_for_p <- levelplot(for_123, par.settings = r_theme, 
		main = "Forest assignment by Gaveau forest cover",
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(sabah_pa_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6)) 
	gav_for_p
	
	# ----------------------
	#  Plot all forest as assigned by CTH forest cover layer
	cth_for_p <- levelplot(cth, par.settings = r_theme, 
		main = "Forest assignment by CTH forest cover",
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(sabah_pa_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6)) +
		layer(sp.polygons(mang_sp, lwd = 1, col = '#d86711')) 
	cth_for_p

	# ----------------------
	#  Plot all forest as assigned by CTH forest cover layer
	cth_p <- levelplot(cth, par.settings = r_theme, 
		main = "Forest assignment by CTH forest cover",
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(sabah_pa_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6))
	cth_p

	
	
# =============================================================================	
###############################################################################






















	
	# ----------------------
	#  Gaveau forest cover layer (for comparison)
	gav <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/raw_spat_data/forest_cover/REGBorneo_ForestCover_2016_CIFOR.tif")




# ============================================================================
#  Rasters must align in all aspects to compare - but origin and extents do not align, so....
# =============================================================================

	# ----------------------
	#   Resample, crop and mask Gaveau forest cover raster to match extent of CTH raster.
	for_rs <- resample(gav, cth, "bilinear")
	ex <- extent(cth)
	for_c <- raster::crop(for_rs, ex)
	for_m <- raster::mask(for_c, cth)


# =============================================================================
#  Compare Gaveau and CTH forest cover.
#   Reclassfiy Gaveau forest cover layer so that intact, logged, and regrowth forest are all assigned
#   as forest (1) and all other cover types are assigned as non-forest (0).
# =============================================================================

	# ----------------------
	#  Reclassification values for Gaveau forest cover layer:
	#   0 = non-forest
	#   1 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_for <- c(0, 0.9999, 0,
		1, 3.00001, 1,
		3.00002, Inf, 0)
	rc_mat_for <- matrix(rc_val_for, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	for_123 <- raster::reclassify(for_m, rc_mat_for, inclulde.lowest = TRUE)


