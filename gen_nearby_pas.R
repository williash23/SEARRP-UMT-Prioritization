###############################################################################
#  Generate layer of protected areas that are near the Sabah border.
#  December 27, 2017; last updated April 9, 2018
#  Script to define which protected areas in Sarawak and Kalimantan are within a certain
#   buffer distance of the Sabah border for use within connectivity assessments.
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
library(units)
	

# =============================================================================
#  Load input data: protected areas of Malaysia and Indonesia from WDAP that are cropped
#   to the borders of Sarawak and Kalimantan.
# =============================================================================

	# ----------------------
	#  All protected areas in Sarawak and Kalimantan.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/sarawak_pa_sf.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/kali_pa_sf.Rdata")
	
	# ----------------------
	#  Region/country borders.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sarawak_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_kali_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/main_sabah_sf.Rdata")
	
	# ----------------------
	#  Existing protected areas ("locked in") in Sabah - forest reserve class I, wildlife sanctuaries, 
	#   virgin jungle reserves, and parks.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sabah_pa_sf.Rdata")

	# ----------------------
	#  All forest reserve areas (class I - VII) in Sabah.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/for_res_rt_sf.Rdata")
	
	

# =============================================================================
#  Combine all protected areas into a single object.
# =============================================================================
	
	# ----------------------
	#  Make all PA objects have same columns.
	sarawak_pa_sf_tmp <- sarawak_pa_sf %>%
		dplyr::select(geometry, name = NAME, status = DESIG) %>%
		mutate(locked_in_val = 1)  %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, c("area_h", "unit"), sep = " ", extra = "drop", fill = "right") %>%
		dplyr::select(-unit)
	sarawak_pa_sf_tmp$area_h <- as.numeric(sarawak_pa_sf_tmp$area_h)
	
	kali_pa_sf_tmp <- kali_pa_sf %>%
		dplyr::select(geometry, name = NAME, status = DESIG) %>%
		mutate(locked_in_val = 1) %>%
		mutate(area_tmp = st_area(.) * 0.0001) %>%
		separate(area_tmp, c("area_h", "unit"), sep = " ", extra = "drop", fill = "right") %>%
		dplyr::select(-unit)
	kali_pa_sf_tmp$area_h <- as.numeric(kali_pa_sf_tmp$area_h)

	# ----------------------
	#  Use rbind to compile by just adding rows together.
	ssk_pa_sf_tmp <- rbind(sabah_pa_sf, sarawak_pa_sf_tmp, kali_pa_sf_tmp)
	
	
	
# =============================================================================
#  Give each PA a unique ID.
# =============================================================================
	
	# ----------------------
	#  Union all features to dissolve internal boundaries.
	ssk_pa_sf <- st_union(ssk_pa_sf_tmp, by_feature = TRUE) %>%
		mutate(id = row_number())
	ssk_pa_sp <- as(ssk_pa_sf, "Spatial")
	
	# ----------------------
	#  Save as .Rdata objects.
	save(ssk_pa_sf, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_sf.Rdata")



# =============================================================================
#  Buffer Sabah to find which nearby PAs will be used in connectivity assessment.
# =============================================================================

	# ----------------------
	#  Set buffer distance.
	buff <- 200000
	
	# ----------------------
	#  Buffer Saba border.
	main_sabah_buff_sf <- st_buffer(main_sabah_sf, buff)
	
	# ----------------------
	#  Crop combined protected areas from Sabah, Sarawak, and Kalimantan to the 
	#   buffered border.
	ssk_pa_near_sf <- st_intersection(ssk_pa_sf, main_sabah_buff_sf)
	ssk_pa_near_sp <- as(ssk_pa_near_sf, "Spatial")
	
	# ----------------------
	#  Save as .Rdata objects.
	# save(ssk_pa_near_sf, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf.Rdata")

	# ----------------------
	#  Function to combine overlapping protected areas 
	clusterSF <- function(sfpolys, thresh){

		dmat = st_distance(sfpolys)
		hc = hclust(as.dist(dmat>thresh), method="single")
		groups = cutree(hc, h=0.5)
		
		d = st_sf(
			geom = do.call(c,
				lapply(1:max(groups), function(g){
					st_union(sfpolys[groups==g,])
				})
				)
			)
		d$group = 1:nrow(d)
		d
	}

	# ----------------------
	#  Set threshold distance and run function for each forest cover type.
	#thresh_d <- set_units(1000, m)
	thresh_d <- set_units(5, m)
 
	# ----------------------
	#  Run function.
	#ssk_pa_near_sf_clust_1k <- clusterSF(ssk_pa_near_sf, thresh_d)
	ssk_pa_near_sf_clust <- clusterSF(ssk_pa_near_sf, thresh_d)
	
	# ----------------------
	#  Save as .Rdata objects.
	# save(ssk_pa_near_sf_clust_1k, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust_1k.Rdata")
	# save(ssk_pa_near_sf_clust, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust.Rdata")

	
	
# =============================================================================
#  Make boundary around Maliau/Danum complex to avoid internal connections.
# =============================================================================

	mal_dan_bound <- ssk_pa_near_sf_clust_1k %>%
		dplyr::filter(group == 18) %>%
		st_buffer(950) %>%
		st_union()
	# save(mal_dan_bound, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/mal_dan_bound.Rdata")


	
# =============================================================================
#  Plots.
# =============================================================================	
	
	
	# ----------------------
	#  Plot.
	all_pa_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = sarawak_pa_sf, fill = "grey30", colour = "grey30") +
		#geom_sf(data = kali_pa_sf, fill = "grey30", colour = "grey30") +
		#geom_sf(data = for_res_rt_sf, aes(fill = factor(NOV_2016)),  colour = "grey50") +
 		#scale_fill_brewer(type = "qual", palette = "Greens",
			#labels = c("	Class I", "	Class II", "	Class III", "	Class IV", "	Class V", "	Class VI", "	ClassVII"), 
			#name = "Forest Reserves") +
		geom_sf(data = sabah_pa_sf, aes(fill = factor(status)), colour = "grey50") +
		#geom_sf(data = conn_sf, fill = "darkred", colour = "darkred") +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	all_pa_p
	
	