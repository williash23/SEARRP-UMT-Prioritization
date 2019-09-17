
# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
library(smoothr)

setwd("C:/Users/saraw/Documents/Prioritization/")


sol_s <- stack("scenario_outputs/12_1_18/all_sol_stack_181201.grd")
sol_no_conn <- stack("scenario_outputs/12_1_18/no_conn_sol.tif")  # for no connectivity outputs


load("study_area_boundaries_and_pas/main_sabah_sf.Rdata")
pu_r <- raster("mainland_sabah_planning_units.tif")
main_sabah_sp <- as(main_sabah_sf, "Spatial")

load("study_area_boundaries_and_pas/sabah_tpa.Rdata")
main_sabah_tpa <- st_intersection(sabah_tpa, main_sabah_sf)
main_sabah_tpa_sp <- as(main_sabah_tpa, "Spatial")



# With PA
tmp1 <- sol_s[[1]]  
#tmp1 <- sol_no_conn # for no connectivity outputs
tmp2 <- rasterToPolygons(sol_w_all_pa, fun = function(x){x>0}, dissolve = TRUE)


tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	smooth(method = "chaikin")
	st_buffer(400) %>%
	#st_buffer(340) %>% # for no connectivity outputs
	fill_holes(area_thresh2)
st_area(tmp4) - sum(st_area(main_sabah_tpa))

#all_inputs_prioritization_w_pa <- as(tmp4, "Spatial")
#no_connectivity_inputs_prioritization_w_pa <- as(tmp4, "Spatial")
#shapefile(all_inputs_prioritization_w_pa, file = "all_inputs_prioritization_w_pa.shp")
#shapefile(no_connectivity_inputs_prioritization_w_pa, file = "no_connectivity_inputs_prioritization_w_pa.shp")

#all_inputs_prioritization <- rasterize(tmp5, pu_r)


#  Mask out existing TPA

tmp5 <- rasterize(as(tmp4, "Spatial"), pu_r)
tmp6 <- raster::mask(tmp5,  main_sabah_tpa_sp, inverse = TRUE)
tmp7 <- rasterToPolygons(tmp6, fun = function(x){x>0}, dissolve = TRUE)

prioritization_wo_existing_PA <- tmp7
shapefile(prioritization_wo_existing_PA, file = "prioritization_wo_existing_PA.shp")



#  Removals
setwd("C:/Users/saraw/Documents/Prioritization/")

sol_s <- stack("scenario_outputs/12_1_18/all_sol_stack_181201.grd")

tmp1 <- raster::mask(sol_s[[2]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_vert_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_vert_input_prioritization, file = "no_vert_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[3]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_fly_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_fly_input_prioritization, file = "no_fly_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[4]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_plant_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_plant_input_prioritization, file = "no_plant_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[5]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_form_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_form_input_prioritization, file = "no_form_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[6]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_cond_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_cond_input_prioritization, file = "no_cond_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[7]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_acd_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_acd_input_prioritization, file = "no_acd_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[8]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_corr_input_prioritization <- as(tmp4, "Spatial")
#shapefile(no_corr_input_prioritization, file = "no_corr_input_prioritization.shp")





sol_rems <- sol_s
all_overlap <- raster::calc(sol_rems, sum)

tmp1 <- raster::mask(all_overlap, main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

prioritizations_overlap_all_val <- as(tmp4, "Spatial")
#shapefile(prioritizations_overlap_all_val, file = "prioritizations_overlap_all_val.shp")




sol_rems <- sol_s
all_overlap <- raster::calc(sol_rems, sum)
all_overlap[all_overlap < 8] =  0

tmp1 <- raster::mask(all_overlap, main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

prioritizations_overlap <- as(tmp4, "Spatial")
#shapefile(prioritizations_overlap, file = "prioritizations_overlap.shp")
















# Zonation output

setwd("C:/Users/saraw/Documents/Zon_Out/sabah1/")
tmp <- raster("outputs/sabah1.TBF_E.wrscr.compressed.tif")

tmp1 <- tmp
tmp1[tmp1 > 0.0000968] <- 1
tmp1[tmp1 < 1] <- 0
cellStats(tmp1, sum)
zon1 <- tmp1

setwd("C:/Users/saraw/Documents/Prioritization/")

load("study_area_boundaries_and_pas/main_sabah_sf.Rdata")
pu_r <- raster("mainland_sabah_planning_units.tif")
main_sabah_sp <- as(main_sabah_sf, "Spatial")

load("study_area_boundaries_and_pas/sabah_tpa.Rdata")
main_sabah_tpa <- st_intersection(sabah_tpa, main_sabah_sf)
main_sabah_tpa_sp <- as(main_sabah_tpa, "Spatial")

main_tpa_u <- main_sabah_tpa %>%st_buffer(1) %>% st_union()
main_tpa_b <- main_tpa_u %>% st_buffer(100)
main_tpa_b_sp <- as(main_tpa_b, "Spatial")

tmp1 <-  raster::mask(zon1, main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)


tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(3000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) #%>%
	#st_buffer(400) %>%
	#fill_holes(area_thresh2)
st_area(tmp4)

	sol <- tmp4
	
	sel <- sol %>%
		mutate(IN = "2") %>%
		dplyr::select(IN) %>%
		st_transform(st_crs(tpa))
	p <- rbind(tpa, sel)

	map1 <- ggplot() +
		geom_sf(data = sabah_border, colour = "grey80", fill = "transparent") +
		geom_sf(data = sarawak_border, colour = "grey80", fill = "transparent") +
		geom_sf(data = kali_border, colour = "grey80", fill = "transparent") +
		geom_sf(data = p,  aes(colour = IN, fill = IN)) +
		scale_fill_manual(name = "", values = c("grey20", "grey50"), 
			labels = c("Existing TPA", "Prioritized Area")) +
		scale_colour_manual(guide = FALSE, values = c("grey20", "grey50")) +
		geom_sf(data = tpa, color = "grey20", fill = "grey20") +
		coord_sf(crs = st_crs(32650)) +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		ylab("Latitude") +
		xlab("Longitude") +
		theme_bw() #+
		#theme(legend.direction = "horizontal", legend.position = "bottom")
	map1
	





	
	
	
	
	
	
	
	
	
	
	

all_inputs_prioritization <- as(tmp4, "Spatial")
shapefile(all_inputs_prioritization, file = "all_inputs_prioritization.shp")
#all_inputs_prioritization <- rasterize(tmp5, pu_r)


# With PA
tmp1 <- zon1
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)


tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) #%>%
	#st_buffer(400) %>%
	#fill_holes(area_thresh2)
st_area(tmp4)

all_inputs_prioritization_w_pa <- as(tmp4, "Spatial")
shapefile(all_inputs_prioritization_w_pa, file = "all_inputs_prioritization_w_pa.shp")

#all_inputs_prioritization <- rasterize(tmp5, pu_r)
setwd("C:/Users/saraw/Documents/Prioritization/")


setwd("C:/Users/saraw/Documents/Prioritization/")





sol_s <- stack("scenario_outputs/12_1_18/all_sol_stack_181201.grd")

tmp1 <- raster::mask(sol_s[[2]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_vert_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_vert_input_prioritization, file = "no_vert_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[3]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_fly_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_fly_input_prioritization, file = "no_fly_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[4]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_plant_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_plant_input_prioritization, file = "no_plant_input_prioritization.shp")


tmp1 <- raster::mask(sol_s[[5]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_form_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_form_input_prioritization, file = "no_form_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[6]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_cond_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_cond_input_prioritization, file = "no_cond_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[7]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_acd_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_acd_input_prioritization, file = "no_acd_input_prioritization.shp")



tmp1 <- raster::mask(sol_s[[8]], main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

no_corr_input_prioritization <- as(tmp4, "Spatial")
shapefile(no_corr_input_prioritization, file = "no_corr_input_prioritization.shp")





sol_rems <- sol_s
all_overlap <- raster::calc(sol_rems, sum)

tmp1 <- raster::mask(all_overlap, main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

prioritizations_overlap_all_val <- as(tmp4, "Spatial")
shapefile(prioritizations_overlap_all_val, file = "prioritizations_overlap_all_val.shp")




sol_rems <- sol_s
all_overlap <- raster::calc(sol_rems, sum)
all_overlap[all_overlap < 8] =  0

tmp1 <- raster::mask(all_overlap, main_tpa_b_sp, inverse = TRUE)
tmp2 <- rasterToPolygons(tmp1, fun = function(x){x>0}, dissolve = TRUE)

tmp3 <- st_as_sf(tmp2)
area_thresh1 <- units::set_units(1000, ha)
area_thresh2 <-  units::set_units(1000, ha)
tmp4 <- tmp3 %>%
	drop_crumbs(area_thresh1) %>%
	st_buffer(400) %>%
	fill_holes(area_thresh2)
st_area(tmp4)

prioritizations_overlap <- as(tmp4, "Spatial")
shapefile(prioritizations_overlap, file = "prioritizations_overlap.shp")

