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