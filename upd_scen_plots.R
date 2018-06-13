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
# install.packages("C:/gurobi751/win64/R/gurobi_7.5-1.zip", repos = NULL)
library(gurobi)
#  devtools::install_github("prioritizr/priortizr") # can use official CRAN version
library(prioritizr)
library(ggplot2)


# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Boundaries
	load(file ="C:/Users/saraw/Desktop/boundaries/sabah_border.Rdata")
	load(file ="C:/Users/saraw/Desktop/boundaries/kali_border.Rdata")
	load(file ="C:/Users/saraw/Desktop/boundaries/sarawak_border.Rdata")
	load(file = "C:/Users/saraw/Desktop/boundaries/sabah_tpa.Rdata")
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/sfi_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/idris_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/deram_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/sega_grid.Rdata")
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/crock_kina_grid.Rdata")
	
	# ----------------------
	#  Prioritization outputs
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen1_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen1_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4_blm_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen5/scen5_sf.Rdata")
	load(file = "C:/Users/saraw/Desktop/scenario_outputs/scen5/scen5_blm_sf.Rdata")
	
	
	
# =============================================================================
#  Plots.
# =============================================================================		
	
	
	#  All scenarios run for both:
	#     1. No BLM constraints
	#     2. BLM constraint

	
	### Scenario 1 ###
	#  Conditions:
	#     1. No SFI/Idris rounds
	#     2. Total area for overall prioritization 410,000 ha
	#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
	#          for 410,000 ha across all Sabah

	# ----------------------
	#  Map 1
	map1 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen1_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 1 (Background) \nNo BLM") +
		theme_bw()
	map1
	
	# ----------------------
	#  Map 1 with BLM constraint
	map1_blm <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen1_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 1 (Background) \nWith BLM") +
		theme_bw()
	map1_blm
	
	
	
	### Scenario 2 ###
	#  Conditions:
	#     1. No SFI/Idris rounds
	#     2. Deramakot locked in 
	#     3. Segaliud Lokan Forest Reserve not locked in but available for selection in solution.
	#     4. Area between Crocker Range and Kinabalu Park is locked out   
	#     5. Total area for overall prioritization 410,000 ha
	#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
	#          for 410,000 ha across all Sabah
	
	# ----------------------
	#  Map 2
	map2 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen2_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 2 \nNo BLM") +
		theme_bw()
	map2
	
	# ----------------------
	#  Map 2 with BLM constraint
	map2_blm <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen2_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 2 \nWith BLM") +
		theme_bw()
	map2_blm

	

	### Scenario 3 ###
	#  Conditions:
	#     1. SFI/Idris:
	#        - SFI- 50,000 ha
	#        - Idris - 100,000 ha
	#        - 6 individual inputs prioritized separately for each area then prioritized together for each
	#          area for a total ha in each
	#     2. Together, SFI, Idris and Deramakot became locked in for overall prioritization
	#     3. Segaliud Lokan Forest Reserve not locked in but available for selection in solution.
	#     4. Area between Crocker Range and Kinabalu Park is locked out   
	#     5. Total area for overall prioritization 410,000 ha
	#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
	#          for 410,000 ha across all Sabah

	# ----------------------
	#  Map 3
	map3 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen3_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 3 \nNo BLM") +
		theme_bw()
	map3

	# ----------------------
	#  Map 3 with BLM constraint
	map3_blm <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen3_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 3 \nWith BLM") +
		theme_bw()
	map3_blm
	
	
	
	### Scenario 4 ###
	#  Conditions:
	#     1. SFI/Idris:
	#        - SFI- 50,000 ha
	#        - Idris - 100,000 ha
	#        - 6 individual inputs prioritized separately for each area then prioritized together for each
	#          area for a total ha in each
	#     2. Together, SFI, Idris and Deramakot became locked in for overall prioritization
	#     3. Segaliud Lokan Forest Reserve locked in and added to locked in area for overall prioritization
	#     4. Area between Crocker Range and Kinabalu Park is locked out   
	#     5. Total area for overall prioritization 410,000 ha
	#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
	#          for 410,000 ha across all Sabah

	# ----------------------
	#  Map 4
	map4 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = sega_grid,  fill = "tomato2",  colour = "tomato2", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen4_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 4 \nNo BLM") +
		theme_bw()
	map4
	
	# ----------------------
	#  Map 4 with BLM constraint
	map4_blm <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		geom_sf(data = sega_grid,  fill = "tomato2",  colour = "tomato2", alpha = 0.7) +
		geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen4_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 4 \nWith BLM") +
		theme_bw()
	map4_blm
	
	
	
	### Scenario 5###
	#  Conditions:
	#     1. Total area for overall prioritization 454,000 ha
	#        - 6 individual inputs prioritized separately for whole of Sabah then prioritized together 
	#          for 454,000 ha across all Sabah

	# ----------------------
	#  Map 5
	map5 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		#geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		#geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		#geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		#geom_sf(data = sega_grid,  fill = "tomato2",  colour = "tomato2", alpha = 0.7) +
		#geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen5_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 5 \nNo BLM") +
		theme_bw()
	map5
	
	# ----------------------
	#  Map 5 with BLM constraint
	map5_blm <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = sarawak_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey30", fill = "transparent") +
		geom_sf(data = sabah_tpa,  colour = "#185b27", fill = "#185b27", alpha = 0.7) +
		#geom_sf(data = sfi_grid, fill = "mediumpurple4", colour = "mediumpurple4", alpha = 0.7) +
		#geom_sf(data = idris_grid, fill = "mediumpurple3", colour = "mediumpurple3", alpha = 0.7) +
		#geom_sf(data = crock_kina_grid, fill = "darkslategray4", colour = "darkslategray4", alpha = 0.7) +
		#geom_sf(data = sega_grid,  fill = "tomato2",  colour = "tomato2", alpha = 0.7) +
		#geom_sf(data = deram_grid,  fill = "tomato4",  colour = "tomato4", alpha = 0.7) +
		geom_sf(data = scen5_blm_sf, fill = "goldenrod3", colour = "goldenrod3", alpha = 0.8) +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ggtitle("Scenario 5 \nWith BLM") +
		theme_bw()
	map5_blm
	
	
	
# =============================================================================
#  Write prioritization solutions to rasters.
# =============================================================================		

	# ----------------------
	#  Create raster following template of study area 
	fly_feat_in_single <- stack("C:/Users/saraw/Desktop/feature_inputs/fly_all.grd")
	temp <- fly_feat_in_single[[1]]
	r_mat <- matrix(0, nrow(temp), ncol(temp))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(temp)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Scenario 1 (Background)
	scen1_sp <- as(scen1_sf, "Spatial")
	scen1_blm_sp <- as(scen1_blm_sf, "Spatial")
	scen1_r <- rasterize(scen1_sp, r_template, field = scen1_sp$ras_val)
	scen1_blm_r <- rasterize(scen1_blm_sp, r_template, field = scen1_blm_sp$ras_val)
	writeRaster(scen1_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen_1.tif")
	writeRaster(scen1_blm_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen1/scen_1_blm.tif")
	
	# ----------------------
	#  Scenario 2
	scen2_sp <- as(scen2_sf, "Spatial")
	scen2_blm_sp <- as(scen2_blm_sf, "Spatial")
	scen2_r <- rasterize(scen2_sp, r_template, field = scen2_sp$ras_val)
	scen2_blm_r <- rasterize(scen2_blm_sp, r_template, field = scen2_blm_sp$ras_val)
	writeRaster(scen2_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2.tif")
	writeRaster(scen2_blm_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen2/scen2_blm.tif")
	
	# ----------------------
	#  Scenario 3
	scen3_sp <- as(scen3_sf, "Spatial")
	scen3_blm_sp <- as(scen3_blm_sf, "Spatial")
	scen3_r <- rasterize(scen3_sp, r_template, field = scen3_sp$ras_val)
	scen3_blm_r <- rasterize(scen3_blm_sp, r_template, field = scen3_blm_sp$ras_val)
	writeRaster(scen3_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3.tif")
	writeRaster(scen3_blm_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen3/scen3_blm.tif")
	
	# ----------------------
	#  Scenario 4
	scen4_sp <- as(scen4_sf, "Spatial")
	scen4_blm_sp <- as(scen4_blm_sf, "Spatial")
	scen4_r <- rasterize(scen4_sp, r_template, field = scen4_sp$ras_val)
	scen4_blm_r <- rasterize(scen4_blm_sp, r_template, field = scen4_blm_sp$ras_val)
	writeRaster(scen4_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4.tif")
	writeRaster(scen4_blm_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen4/scen4_blm.tif")
	
	# ----------------------
	#  Scenario 5
	scen5_sp <- as(scen5_sf, "Spatial")
	scen5_blm_sp <- as(scen5_blm_sf, "Spatial")
	scen5_r <- rasterize(scen5_sp, r_template, field = scen5_sp$ras_val)
	scen5_blm_r <- rasterize(scen5_blm_sp, r_template, field = scen5_blm_sp$ras_val)
	writeRaster(scen5_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen5/scen5.tif")
	writeRaster(scen5_blm_r, file = "C:/Users/saraw/Desktop/scenario_outputs/scen5/scen5_blm.tif")
	
