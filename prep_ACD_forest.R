###############################################################################
#  Reconcile forest cover layers.
#  November 28, 2017; last updated February 2, 2018
#   Script to reclassify cells identified as forest and non-forest from CAO above-ground carbon
#   density (ACD).
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
#   Set value for carbon storage limit that corresponds to forest
#   So far, I've tested 60mg and 40 mg, which seem to be lower limits according to text and figures 
#   from Asner et al. publication text and figures.
# =============================================================================

	# ----------------------
	#  Cutoff value
	acd_cutoff <- 40	
	
	# ----------------------
	#  Text of cutoff value to be used in plotting
	acd_cutoff_text <- as.name(acd_cutoff)

	# ----------------------
	#  Aggregation factor
	agg_fact <- 30
	
	

# =============================================================================
#  Load input data.
# =============================================================================
	
	# ----------------------
	#  Load forest reserves, to remove mangroves
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/for_res_rt_sf.Rdata")
	
	# ----------------------
	#  Load unmasked carbon map from Asner et al.
	acd <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/raw_spat_data/CAO_ACD_30m_unmasked.tif")
	acd_agg <- raster::aggregate(acd, agg_fact)
	
	# ----------------------
	#  Forest reserves
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/for_res_rt_sf.Rdata")

	# ----------------------
	#  Boundaries
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sabah_d.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sarawak_d.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_kali_d.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/main_sabah_sf.Rdata")
	
	

# =============================================================================
#  Reclassify carbon raster layer so that values below carbon value decided above are assigned as 
#   non-forest (0) and values above are assigned as forest (1).
# =============================================================================

	# ----------------------
	#  Reclassification values for carbon layer:
	#   0 = non-forest
	#   1 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_acd <- c(-Inf, acd_cutoff, 0,
		acd_cutoff , Inf, 1)
	rc_mat_acd <- matrix(rc_val_acd, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	acd_rc <- raster::reclassify(acd, rc_mat_acd, include.lowest = TRUE)
	acd_rc_agg <- raster::reclassify(acd_agg, rc_mat_acd, include.lowest = TRUE) # aggregated



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
	acd_rc_no_mang <- raster::mask(acd_rc, mang_sp, updatevalue = 0, inverse = TRUE)
	acd_rc_no_mang_agg <- raster::mask(acd_rc_agg, mang_sp, updatevalue = 0, inverse = TRUE)

	# ----------------------
	#  Save output as polygon and raster
	# writeRaster(acd_rc, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc.grd")
	# writeRaster(acd_rc_agg, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc_agg.grd")
	# writeRaster(acd_rc_no_mang, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc_no_mang.grd")
	# writeRaster(acd_rc_no_mang_agg, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc_no_mang_agg.grd")



# =============================================================================
#  Polygonize to use with planning unit sp objects.
# =============================================================================	

	# ----------------------	
	# Polygonize
	acd_p <- polygonizer(acd_rc_no_mang)
	acd_sf <- st_as_sf(acd_p)
	acd_sf_for <- acd_sf %>%
		dplyr::filter(DN == 1)
		
	acd_agg_p <- polygonizer(acd_rc_no_mang_agg)
	acd_agg_sf <- st_as_sf(acd_agg_p)
	acd_agg_sf_for <- acd_agg_sf %>%
		dplyr::filter(DN == 1)

	# ----------------------
	#  Save output as polygon 
	#save(acd_sf_for, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_sf_for.Rdata")
	#save(acd_agg_sf_for, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_agg_sf_for.Rdata")

	
	

	
	
# =============================================================================
#  Extend forest layer to area just outside Sabah (for use in movement modeling).
# =============================================================================

	# ----------------------
	#  Load forest cover from above (if needed), forest cover from Gaveau (for areas outside Sabah),
	#	 and coast line.
	#acd_rc_no_mang_agg <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc_no_mang_agg.grd")
	add_for <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/land_forest_cover/for_cov_ssk_agg.grd")
	ocean <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/ocean_r.grd")

	# ----------------------
	#  Match the Gaveau layer to the extent and resolution of the ACD layer (for smaller forest
	#  connectivity cover layer).
	#add_for_rs <- raster::resample(add_for, acd_rc_no_mang_agg)
	
	# ----------------------
	#  OR, to extend the area of forest further down  into Sarawak use the larger buffer and resample and
	#   reclassify the ACD layer to match the Gaveau layer.
	acd_rs <- raster::resample(acd, add_for)
	
		# ----------------------
		#  Set reclassification values in matrix format.
		rc_val_acd <- c(-Inf, acd_cutoff, 0,
			acd_cutoff , Inf, 1)
		rc_mat_acd <- matrix(rc_val_acd, 
			ncol = 3, 
			byrow = TRUE)
			
		# ----------------------
		#  Reclassify function.
		acd_rs_rc <- raster::reclassify(acd_rs, rc_mat_acd, include.lowest = TRUE) # aggregated

	# ----------------------
	#  Generate desired buffer distance around Sabah. 
	#sabah_buff_sf <- st_buffer(main_sabah_sf, 15000)
	sabah_buff_sf <- st_buffer(main_sabah_sf, 100000)
	sabah_buff_sp <- as(sabah_buff_sf, 'Spatial')
	
	# ----------------------
	#  Reclassify forest from Gaveau to add to ACD forest.
	rc_val_add <- c(0, 3.1, 1,
		3.15, Inf, 0)
	rc_mat_add <- matrix(rc_val_add, 
		ncol = 3, 
		byrow = TRUE)
	#add_for_rc <- raster::reclassify(add_for_rs, rc_mat_add)
	add_for_rc <- raster::reclassify(add_for, rc_mat_add)

	# ----------------------
	#  Mask Gaveau forest cover (for area outside Sabah)
	add_for_m1 <- raster::mask(add_for_rc, sabah_buff_sp)
	
	# ----------------------
	#  Mask out area of Sabah.
	add_for_m2 <- raster::mask(add_for_m1, border_sabah_d, inverse = TRUE)
	
	# ----------------------
	#  Overlay forest rasters.
	add_for_m2[is.na(add_for_m2)] <- 0
	#acd_rc_no_mang_agg[is.na(acd_rc_no_mang_agg)] <- 0
	acd_rs_rc[is.na(acd_rs_rc)] <- 0
	#extended_for_cov <- overlay(acd_rc_no_mang_agg, add_for_m2, fun=function(x,y){return(x+y)})
	extended_for_cov <- overlay(acd_rs_rc, add_for_m2, fun=function(x,y){return(x+y)})
	
	# ----------------------
	#  Reclassify extended forest cover.
	rc_val_ext <- c(0, 0,
		1, 1,
		2, 1)
	rc_mat_ext <- matrix(rc_val_ext, 
		ncol = 2, 
		byrow = TRUE)
	extended_for_cov_rc <- raster::reclassify(extended_for_cov, rc_mat_ext)
	
	# ----------------------
	#  Resample ocean to match ACD forest cover. 
	#ocean_rs <- raster::resample(ocean, acd_rc_no_mang_agg)
	ocean_rs <- raster::resample(ocean, acd_rs_rc)
	
	# ----------------------
	#  Reclassify ocean.
	rc_val_ocean <- c(-Inf, 75, 0,
		75, Inf, 3)
	rc_mat_ocean <- matrix(rc_val_ocean, 
		ncol = 3, 
		byrow = TRUE)
	ocean_rc <- raster::reclassify(ocean_rs, rc_mat_ocean, include.lowest = TRUE)
	
	# ----------------------
	#  Mask area outside Sabah to NA except for buffer of forest.
	ext_for_cov_ocean <-overlay(extended_for_cov_rc, ocean_rc, fun=function(x,y){return(x+y)})
	
	# ----------------------
	#  Reclassify new forest and ocean.
	# 0 = non-forest
	# 1 = forest
	# 3 = ocean
	rc_val_for_oc <- c(0, 0,
		1, 1,
		3, 3, 
		4, 3)
	rc_mat_for_oc <- matrix(rc_val_for_oc, 
		ncol = 2, 
		byrow = TRUE)
	ext_for_cov_ocean_rc <- raster::reclassify(ext_for_cov_ocean, rc_mat_for_oc, include.lowest = TRUE)
	
	# ----------------------
	#  Save to raster.
	# writeRaster(ext_for_cov_ocean_rc, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ext_for_cov_ocean_rc.grd")

	
	
# =============================================================================
#  Rasterize existing protected areas and incorporate into forest cover for connectivity.
# =============================================================================
	
	
	# ----------------------
	#  "Locked in" already existing protected areas.
	#load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf.Rdata")
	ssk_pa_near_sp  <- as(ssk_pa_near_sf, "Spatial")
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	#r_mat <- matrix(NA, nrow(ext_for_cov_ocean_rc), ncol(ext_for_cov_ocean_rc))
	r_mat <- matrix(NA, nrow(acd_rs_rc), ncol(acd_rs_rc))
	r_template <- raster(r_mat)
	#extent(r_template) <- extent(ext_for_cov_ocean_rc)
	extent(r_template) <- extent(acd_rs_rc)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Rasterize existing protected areas with raster template.
	ssk_pa_near_r <- rasterize(ssk_pa_near_sp, r_template)
	rc_val_pa <- c(cellStats(ssk_pa_near_r, min), cellStats(ssk_pa_near_r, max), 1)
	rc_mat_pa <- matrix(rc_val_pa, 
		ncol = 3, 
		byrow = TRUE)
	pa_rc <- raster::reclassify(ssk_pa_near_r, rc_mat_pa)
	
	# ----------------------
	#  Mask areas of PAs on forest cover.
	# 4 = PA
	ext_for_cov_ocean_pa_rc <- raster::mask(ext_for_cov_ocean_rc, pa_rc,  maskvalue = 1, updatevalue = 4)
	
	# ----------------------
	#  Save to raster.
	# writeRaster(ext_for_cov_ocean_pa_rc, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ext_for_cov_ocean_pa_rc.grd")
	
	
	
# =============================================================================
#  Rasterize roads and incorporate into forest cover for connectivity.
# =============================================================================
	
	# ----------------------
	#  Load main roads (paved).
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/main_rds_rt_sf.Rdata")
	rds_sf <- main_rds_rt_sf %>%
		dplyr::select(value, Status) 
	rds_sf$Status <- as.factor(rds_sf$Status)
	rds_sp <- as(rds_sf, "Spatial")
	
	# ----------------------
	#  Reformat as raster
	rds_r <- raster::rasterize(rds_sp, r_template, 
		field = rds_sp$value, fun = modal)
	rds_r[is.na(rds_r)] <- 0
	
	# ----------------------
	#  Overlay buffered forest cover and roads rasters.
	ext_for_ocean_pa_rds <- overlay(ext_for_cov_ocean_pa_rc, rds_r, fun=function(x,y){return(x+y)})
	
	# ----------------------
	#  Reclassify...forest cover, roads, and ocean.
	rc_val_for_conn <- c(0, 0,
		1, 1,
		3, 3,
		4, 4,
		5, 2,
		6, 2,
		8, 3, 
		9, 2)
	rc_mat_for_conn <- matrix(rc_val_for_conn, 
		ncol = 2, 
		byrow = TRUE)
		
	for_cov_conn <- raster::reclassify(ext_for_ocean_pa_rds, rc_mat_for_conn, include.lowest = TRUE)
	
	
	# ----------------------
	#  Rename larger spatial extent forest cover (from larger Sabah buffer) to not overwrite.
	for_cov_conn_lg <- for_cov_conn
	
	
	# ----------------------
	#  Save to raster.
	# writeRaster(for_cov_conn, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/for_cov_conn.grd")
	writeRaster(for_cov_conn_lg, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/for_cov_conn_lg.grd")
	


	
# =============================================================================
#  Plots
# =============================================================================

	# ----------------------
	#  Load border data for plotting
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sabah_pa_sp.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sabah_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sarawak_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_kali_d.Rdata")
	acd_rc_no_mang_agg <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/acd_rc_no_mang_agg.grd")
	
	# ----------------------
	#  Set color palette themes
	r_theme <- rasterTheme(rev(brewer.pal(4, "Spectral")))
	r_theme$panel.background$col = 'grey80' 
	r_theme$regions$col = c("grey90", "grey90", "grey90", "grey90",
		"palegreen3", "palegreen3", "palegreen3", "palegreen3",
		"palegreen4", "palegreen4", "palegreen4", "palegreen4", 
		"darkslategray", "darkslategray", "darkslategray", "darkslategray")
	
	# ----------------------
	#  Plot forest using levelPlot.
	acd_p <- levelplot(acd_rc_no_mang_agg, par.settings = r_theme, 
		#main = paste("Differences in forest assignment between Gaveau and Carbon ",
		#"(", carb_cutoff_text, "mg cutoff)", sep =""),
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) #+
		#layer(sp.polygons(sabah_pa_sp, lwd = 1.2, col = 'grey50', fill = 'transparent', alpha = 0.6)) 
	acd_p

	# ----------------------
	#  Plot forest using ggplot2.
	fortify.Raster <- function(x, maxPixel = 1000000) {
	   
		if (ncell(x) > maxPixel) {
			x <- sampleRegular(x, maxPixel, asRaster=TRUE)
		}
		xy <- xyFromCell(x, seq_len(ncell(x)))
		out <- x %>%
			getValues() %>%
			data.frame(values = .) %>%
			cbind(xy)
		return(out)
	}

	for_gg <- fortify.Raster(acd_rc_no_mang_agg)
	for_gg[is.na(for_gg)] <- 3
	for_gg$values <- as.factor(for_gg$values)

	# ----------------------
	#  Save plot to see planning units selected over time.
	acd_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_raster(data = for_gg, aes(x, y, fill = values)) +
		scale_fill_manual(values=c("transparent, "#68B378", "transparent")) +
		#geom_sf(data = pu_grid_10k_no_pa_sf, colour = "black", fill = "transparent") +
		geom_sf(data = ssk_pa_near_sf, fill = "darkorange2", colour = "grey50", alpha = 1) +
		#geom_sf(data = all_sol_sf,  aes(fill = factor(scenario))) +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	acd_p
	
	
	
# =============================================================================	
###############################################################################






























	
	
# =============================================================================
#  To overlay ACD and CTH rasters
# =============================================================================	

	acd_elev_l <- overlay(acd_rc, elev_rc_l, fun=function(x,y){return(x*y)})
		acd_cth_by_elev <- overlay(cth_by_elev, acd_elev_l, fun=function(x,y){return(x+y)})















	# ----------------------
	#  Forest cover layer from Gaveau et al. 
	forest <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/for_cov_m.grd")




# =============================================================================
#  Rasters must align in all aspects to compare - but origin and extents do not align, so....
# =============================================================================

	# ----------------------
	#   Extend forest cover raster to match extent of carbon raster
	for_e <- raster::extend(forest, carb_r)

	# ----------------------
	#  Then, resample carbon raster so that origins align (can't resample forest cover raster
	#   because that turns the values from categorical to continuous)
	carb_r <- raster::resample(carb_r, for_e)

	# ----------------------
	#  Finally, make sure NA pixels in the carbon map raster are set to 0 so that calculations can be 
	#   completed. These correspond to the "masked out" areas of the carbon map.
	carb_r[is.na(carb_r[])] <- 0 



# =============================================================================
#  Reclassify carbon raster layer so that values below carbon value decided above are assigned as 
#   non-forest (0) and values above are assigned as forest (5).
# =============================================================================

	# ----------------------
	#  Reclassification values for carbon layer:
	#   0 = non-forest
	#   5 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_carb <- c(0, carb_cutoff - 0.00001, 0,
		carb_cutoff  - (9e-06), 501, 5)
	rc_mat_carb <- matrix(rc_val_carb, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	carb_rc <- raster::reclassify(carb_r, rc_mat_carb, inclulde.lowest = TRUE)



# =============================================================================
#  Reclassfiy Gaveau forest cover layer so that intact, logged, and regrowth forest are all assigned
#   as forest (1) and all other cover types are assigned as non-forest (0).
# =============================================================================

	# ----------------------
	#  Reclassification values for Gaveau forest cover layer:
	#   0 = non-forest
	#   1 = forest
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_for <- c(0, 0,
		1, 1,
		2, 1,
		3, 1,
		4, 0,
		5, 0)
	rc_mat_for <- matrix(rc_val_for, 
		ncol = 2, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	for_123 <- raster::reclassify(for_e, rc_mat_for)



# =============================================================================
#  Compare Gaveau forest cover and carbon map rasters and reclassify to determine where 
#   forest categories align and where they do not.
# =============================================================================

	# ----------------------
	#  Subtract Gaveau forest cover layer from carbon layer.
	#  Possible outcome values: 
	#   -5 = Gaveau non-forest - carbon forest
	#   -4 = Gaveau forest - carbon forest
	#   0 = Gaveau non-forest - carbon non-forest
	#   1 = Gaveau forest - carbon non-forest
	for_carb_tmp <- for_123 - carb_rc

	# ----------------------
	#  Reclassify outcomes of Gaveau forest cover layer - carbon layer
	#  Reclassification values:
	#   0 = Both assigned as non-forest (0 above)
	#   1 = Both assigned as forest (-4 above)
	#   2 = Gaveau assigned as non-forest but carbon assigned as forest (-5 above)
	#   3 = Gaveau assigned as forest but carbon assigned as non-forest (1 above)
	
	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_for_carb <- c(-Inf, -4.9, 2,
		-4.1, -3.9, 1,
		-0.9, 0.5, 0,
		0.9, 2, 3)
	rc_mat_for_carb <- matrix(rc_val_for_carb, 
		ncol = 3, 
		byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	for_carb <- raster::reclassify(for_carb_tmp , rc_mat_for_carb, right = FALSE)



# =============================================================================
#  Create and save new raster that sets any forest deisgnation (i.e., 1, 2, or 3 above) as forest (1). 
# =============================================================================

	# ----------------------
	#  Reclassify outcomes of Gaveau forest cover layer - carbon layer
	#  Reclassification values:
	#   0 = Non-forest (0 above)
	#   1 = Forest (1, 2, or 3 above)

	# ----------------------
	#  Set reclassification values in matrix format.
	rc_val_all_for <- c(0, 0,
		1, 1,
		2, 1,
		3, 1)
	rc_mat_all_for <- matrix(rc_val_all_for, 
		ncol = 2, 
		byrow = TRUE)
	
	# ----------------------
	#  Reclassify function.
	all_for <- raster::reclassify(for_carb , rc_mat_all_for)

	# ----------------------
	#  Save output as raster
	#writeRaster(all_for, "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/planning_units/new_forest_cover.grd")



# =============================================================================
#  Plots
# =============================================================================

	# ----------------------
	#  Load border data for plotting
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/planning_units/locked_in_sp.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sabah_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sarawak_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_kali_d.Rdata")

	# ----------------------
	#  Set color palette themes
	r_theme <- rasterTheme(rev(brewer.pal(4, "Spectral")))
	r_theme$panel.background$col = 'grey80' 
	r_theme2 <- rasterTheme(brewer.pal(2, "Spectral"))
	r_theme2$panel.background$col = 'grey80' 
	r_theme2$regions$col = c( "#2b83ba", "#abdda4")

	# ----------------------
	#  Plot differences between Gaveau and carbon map
	for_diff_p <- levelplot(for_carb, par.settings = r_theme, 
		main = paste("Differences in forest assignment between Gaveau and Carbon ",
		"(", carb_cutoff_text, "mg cutoff)", sep =""),
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(locked_in_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6)) 
	for_diff_p

	# ----------------------
	#  Plot all forest as assigned by either Gaveau or carbon map
	all_for_p <- levelplot(all_for, par.settings = r_theme2, 
		main = paste("Forest assignment in either or both Gaveau and Carbon ",
		"(", carb_cutoff_text, "mg cutoff)", sep =""),
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(locked_in_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6)) 
	all_for_p



# =============================================================================
#  Compare Gaveau forest cover and carbon map rasters in opposite direction (for confirmation 
#   that results are the same) and reclassify to determine where forest categories align and where 
#   they do not.
# =============================================================================

	# ----------------------
	# #  Subtract carbon layer from Gaveau forest cover layer 
	# #  Possible outcome values:
	# #   -1 = carbon non-forest - Gaveau forest
	# #   0 = carbon non-forest - Gaveau non-forest
	# #   4 = carbon forest - Gaveau forest
	# #   5 = carbon forest - Gaveau non-forest
	# carb_for_tmp <- carb_rc - for_123

	# ----------------------
	# #  Reclassify outcomes of carbon layer - Gaveau forest cover layer
	# #  Reclassification values:
	# #   0 = Both assigned as non-forest (0 above)
	# #   1 = Both assigned as forest (4 above)
	# #   2 = Gaveau assigned as non-forest but carbon assigned as forest (5 above)
	# #   3 = Gaveau assigned as forest but carbon assigned as non-forest (-1 above)
	
	# ----------------------
	#  Set reclassification values in matrix format.
	# rc_val_carb_for <- c(-Inf, -0.5, 3,
		# -0.4, 0.4, 0,
		# 3.5, 4.5, 1,
		# 4.6, 5.6, 2)
	# rc_mat_carb_for <- matrix(rc_val_carb_for, 
		# ncol = 3, 
		# byrow = TRUE)
		
	# ----------------------
	#  Reclassify function.
	# carb_for <- raster::reclassify(carb_for_tmp , rc_mat_carb_for, right = FALSE)

	
	
# =============================================================================
#  Plots
# =============================================================================




	# ----------------------
	#  Plot differences between Gaveau and carbon map
	# for_diff_p2 <- levelplot(carb_for, par.settings = r_theme, 
		# main = paste("Differences in forest assignment between Carbon and Gaveau ",
		# "(", carb_cutoff_text, "mg cutoff)", sep =""),
		# xlab= "Longitude (UTM)",
		# ylab="Latitude (UTM)",
		# margin = FALSE) +
		# layer(sp.polygons(border_sabah_d, lwd = 1, col = 'grey40')) + 
		# layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		# layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		# layer(sp.polygons(locked_in_sp, lwd = 1.2, col = 'grey70', fill = 'grey40', alpha = 0.6)) 
	# for_diff_p2
	
	

# =============================================================================	
###############################################################################