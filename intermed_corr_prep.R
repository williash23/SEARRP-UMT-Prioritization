	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	corr_conn_feat_in <- raster("feature_inputs/corr_conn_feat_in_w_pa_range.grd")
	
	moves_sm <- raster("C:/Users/saraw/Desktop/moves_sm.grd")
	moves_sm <- resample(moves_sm, corr_conn_feat_in)
	x_min <- cellStats(moves_sm, min)
	x_max <- cellStats(moves_sm, max)
	moves_01 <- raster::calc(moves_sm, range01)
	
	pt_wt_r <- raster("C:/Users/saraw/Documents/Prioritization/feature_prep/pt_wt_r.grd")
	pt_wt_r <- resample(	pt_wt_r, corr_conn_feat_in)
	x_min <- cellStats(pt_wt_r, min)
	x_max <- cellStats(pt_wt_r, max)
	wt_01 <- raster::calc(pt_wt_r, range01)
	

	
	corr_conn <- stack(wt_01, moves_01, moves_01, moves_01)
	tmp1 <- raster::calc(corr_conn, mean, na.rm = TRUE)
	tmp2 <- stack(corr_conn_feat_in, moves_01)
    corr_conn_new <- raster::calc(tmp2, mean, na.rm = TRUE)
	#corr_conn_new_sm <- raster::focal(corr_conn_new, fun = mean, na.rm = TRUE, w=matrix(1/15, nc=5, nr=5), pad = TRUE)
