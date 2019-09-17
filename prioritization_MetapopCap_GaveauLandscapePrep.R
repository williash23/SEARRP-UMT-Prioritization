
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
 
no_conn_sol <- st_read("C:/Users/saraw/Documents/Prioritization/scenario_outputs/12_1_18/no_connectivity_inputs_prioritization_w_pa.shp") %>%
	mutate(value = 100)
sol <- st_read("C:/Users/saraw/Documents/Prioritization/scenario_outputs/12_1_18/all_inputs_prioritization_w_pa.shp") %>%
	mutate(value = 100)
	
for_tmp <- raster("C:/Users/saraw/Documents/Spatial_data/raw_spat_data/forest_cover/REGBorneo_ForestCover_2016_CIFOR.tif")

 
# =============================================================================
#  Prep forest cover layer.
# =============================================================================
	
	# ----------------------
	# Crop to study area 
	for_sa <- raster::crop(for_tmp, main_sabah_sf)
	for_sa <- raster::mask(for_sa, main_sabah_sf)
	
	# ----------------------
	# Aggregate by factor of 33 (so resolution is 990m x 990m)
	for_sa_a <- raster::aggregate(for_sa, fact = 33, fun = modal)
	#  after cropping to study area, possible values: 1, 2, 3, 4, 5
	# 1 = intact forest
	# 2 = logged forest
	# 3 = regrowth
	# 4 = non-forest
	# 5  = water
	
	
	
# =============================================================================
#  Set up raster template.
# =============================================================================
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(for_sa_a), ncol(for_sa_a))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(for_sa_a)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
	
	
	
# =============================================================================
#  Make raster less patchy
# =============================================================================	
	
	# ----------------------
	#  Drop small polygons  
	area_thresh <- set_units(5, km^2)
	for_sa_a_sf <- st_as_sf(rasterToPolygons(for_sa_a, dissolve = TRUE))  %>%
		st_cast("POLYGON") %>%
		drop_crumbs(area_thresh) 
		
	for_sa_a_r <- fasterize(for_sa_a_sf, r_template, field = "REGBorneo_ForestCover_2016_CIFOR")
	
	# ----------------------
	#  Fill NA raster cells  (created by dropping crumbs) with most common neighbors
	fill.na <- function(x, i=5) {
	  if( is.na(x)[i] ) {
		return(round(modal(x, na.rm=TRUE),0) )
		} else {
			return(round(x[i],0) )
			}
		}  
	
	for_sa_a_fill <- focal(for_sa_a_r, w = matrix(1, 3, 3), fun = fill.na,  pad = TRUE, na.rm = FALSE)
	
	
	


# =============================================================================
#  Add plantation shapefile to raster.
# =============================================================================
	
	# ----------------------
	#  Load plantations (polygon shapefile)
	plan_tmp <- shapefile("C:/Users/saraw/Documents/Spatial_data/raw_spat_data/REGBorneo_OriginOfLandConvertedToITPAndIOPPComplexTrajectory_1973to2016_CIFOR.shp")
	
	# ----------------------
	#  Crop to study area.
	plan_sf <- st_as_sf(plan_tmp) %>%
		dplyr::select(klas, F2016, Country) %>%
		dplyr::mutate(plan_code = 10) %>%
		st_buffer(0.001)
	plan_sf_sa <- st_intersection(plan_sf,  main_sabah_sf)
	plan_sp <- as(plan_sf_sa, "Spatial")
	plan <- st_as_sf(plan_sp)
	
	
	# ----------------------
	#  Rasterize plantation shapefile so it can be combined with forest cover.
	plan_r <- fasterize(plan, r_template, field = "plan_code")
	plan_r[is.na(plan_r)] <- 0
	# possible values: 
	# 0 = not plantation
	# 10 = plantation

	
	# ----------------------
	#  Overlay wiht forest cover.
	for_plan <- raster::overlay(for_sa_a_fill, plan_r, fun = function(x,y) {return(x+y)})
	# possible values: 1, 2, 3, 4, 11, 12, 13, 14
	# 1 = intact forest
	# 2 = logged forest
	# 3 = regrowth
	# 4 = non-forest
	# 11 = intact forest + plantation
	# 12 = logged forest + plantation
	# 13 = regrowth + plantation
	# 14 = non-forest + plantation

	
	# ----------------------
	#  Reclassify.
	rc_val <- c(1, 1,
		2, 2,
		3, 2,
		4, 0,
		11, 3,
		12, 3,
		13, 3,
		14, 3)
	rc_mat <- matrix(rc_val, 
		ncol = 2, 
		byrow = TRUE)
	sa_for_plan <- raster::reclassify(for_plan, rc_mat, include.lowest = TRUE)
	# 0 = non-forest/no data/water
	# 1 = intact forest
	# 2 = logged forest/regrowth
	# 3 = plantation

	
	# ----------------------
	# Rasterize connecticivity and no connectivity solutions
	nc_sol_r <- fasterize(no_conn_sol, r_template, field = "value")
	nc_sol_r[is.na(nc_sol_r)] <- 0
	# Included in solution = 100, not included = 0
	sol_r <- fasterize(sol, r_template, field = "value")
	sol_r[is.na(sol_r)] <- 0
	sol_r[sol_r > 1] <- 200
	# Included in solution = 200, not included = 0
 
	
	# ----------------------
	#  Overaly solution raster with land cover raster (no connectivity solution)
	sa_for_cov_nc <- raster::overlay(sa_for_plan, nc_sol_r, fun = function(x,y, na.rm = TRUE) {return(x+y)})
	# 0 = no data/non-forest                            
	# 1 = intact forest
	# 2 = logged forest/regrowth
	# 3 = plantation
	# 100 = no data/non-forest + NC solution/PA
	# 101 = intact forest + NC solution/PA
	# 102 = logged forest/regrowth + NC solution/PA
	# 103 = plantation + NC solution/PA
	
	
	sa_for_cov_nc_c <- raster::overlay(sa_for_cov_nc, sol_r, fun = function(x,y, na.rm = TRUE) {return(x+y)})
	# 0 = no data/non-forest                            
	# 1 = intact forest
	# 2 = logged forest/regrowth
	# 3 = plantation
	# 100 = no data/non-forest + NC solution/PA
	# 101 = intact forest + NC solution/PA
	# 102 = logged forest/regrowth + NC solution/PA
	# 103 = plantation + NC solution/PA
	# 200 = no data/non-forest + C solution/PA
	# 201 = intact forest + C solution/PA
	# 202 = logged forest/regrowth + C solution/PA
	# 203 = plantation + C solution/PA
	# 300 = NC solution/PA + C solution/PA
	# 301 = intact forest + NC solution/PA + C solution/PA
	# 302 = logged forest/regrowth + NC solution/PA + C solution/PA
	# 303 = plantation + NC solution/PA + C solution/PA
	
	


 
	
# =============================================================================
#  Convert to polygons and simplfy to generate habitat patches.
# =============================================================================
	
	# ----------------------
	#  Convert to polygons.
	sa_landcov_sf_tmp <- spex::polygonize(sa_for_cov_nc_c) 
	cov_types <- sort(unique(sa_landcov_sf_tmp$layer))

	area_thresh <- set_units(2, km^2)
	for(i in cov_types){
		sf_tmp <- sa_landcov_sf_tmp %>%
			dplyr::filter( layer == i)
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
	
	
	area_thresh <- set_units(5, km^2)
	sa_landcov_sf_all_tmp <- rbind(landcov_0, landcov_1, landcov_2, landcov_3, 
		landcov_100, landcov_101, landcov_102, landcov_103,
		landcov_200, landcov_201, landcov_202, landcov_203, 
		landcov_300, landcov_301, landcov_302, landcov_303) %>%
		drop_crumbs(area_thresh)
	
	
	# ----------------------
	#  Rasterize based on new simplified polygons.
	sa_landcov_r_all_tmp <- fasterize(sa_landcov_sf_all_tmp, r_template, field = "value")

	# ----------------------
	#  Fill NA raster cells  (created by dropping crumbs) with most common neighbors
	fill.na <- function(x, i=5) {
	  if( is.na(x)[i] ) {
		return(round(modal(x, na.rm=TRUE),0) )
	  } else {
		return(round(x[i],0) )
	  }
	}  
	sa_landcov_r_fill <- focal(sa_landcov_r_all_tmp, w = matrix(1, 3, 3), fun = fill.na,  pad = TRUE, na.rm = FALSE)
	
	# 0 = no data/non-forest                            
	# 1 = intact forest
	# 2 = logged forest/regrowth
	# 3 = plantation
	# 100 = no data/non-forest + NC solution/PA
	# 101 = intact forest + NC solution/PA
	# 102 = logged forest/regrowth + NC solution/PA
	# 103 = plantation + NC solution/PA
	# 200 = no data/non-forest + C solution/PA
	# 201 = intact forest + C solution/PA
	# 202 = logged forest/regrowth + C solution/PA
	# 203 = plantation + C solution/PA
	# 300 = NC solution/PA + C solution/PA
	# 301 = intact forest + NC solution/PA + C solution/PA
	# 302 = logged forest/regrowth + NC solution/PA + C solution/PA
	# 303 = plantation + NC solution/PA + C solution/PA
	
	sa_landcov_sf <- st_as_sf(rasterToPolygons(sa_landcov_r_fill, dissolve = TRUE)) %>%
		st_cast("POLYGON")
	
	# Reclassification:
	# Non-forest & Plantation = 0
	# PA and prioritized area = 1
	# Intact forest = 1
	# Logged forest/regrowth = 2

	
	# Conn Prior Landscape Patches:
	conn_prior_patches <- sa_landcov_sf
	conn_prior_patches$layer[conn_prior_patches$layer == 0] <- 0
	conn_prior_patches$layer[conn_prior_patches$layer == 1] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 2] <- 2
	conn_prior_patches$layer[conn_prior_patches$layer == 3] <- 0
	conn_prior_patches$layer[conn_prior_patches$layer == 100] <- 0
	conn_prior_patches$layer[conn_prior_patches$layer == 101] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 102] <- 2
	conn_prior_patches$layer[conn_prior_patches$layer == 103] <- 0
	conn_prior_patches$layer[conn_prior_patches$layer == 200] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 201] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 202] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 203] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 300] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 301] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 302] <- 1
	conn_prior_patches$layer[conn_prior_patches$layer == 303] <- 1

	conn_prior_r <- fasterize(conn_prior_patches, r_template, field = "layer")
	
	
	# No Conn Prior Landscape Patches:
	no_conn_prior_patches <- sa_landcov_sf
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 0] <- 0
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 1] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 2] <- 2
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 3] <- 0
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 100] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 101] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 102] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 103] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 200] <- 0
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 201] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 202] <- 2
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 203] <- 0
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 300] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 301] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 302] <- 1
	no_conn_prior_patches$layer[no_conn_prior_patches$layer == 303] <- 1

	no_conn_prior_r <- fasterize(no_conn_prior_patches, r_template, field = "layer")
	
	
	save(conn_prior_patches, file = "conn_prior_patches.Rdata")
	save(no_conn_prior_patches, file = "no_conn_prior_patches.Rdata")
	
	
	writeRaster(conn_prior_r, file = "conn_prior_r.tif", overwrite = TRUE)
	writeRaster(no_conn_prior_r, file = "no_conn_prior_r.tif", overwrite = TRUE)
	
	
	
	
