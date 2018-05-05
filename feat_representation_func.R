feat_rep <- function(x, solution) {
    # assert valid arguments
    assertthat::assert_that(
      inherits(solution, c("SpatialPointsDataFrame", "SpatialLinesDataFrame",
                           "SpatialPolygonsDataFrame")),
      number_of_zones(x) == ncol(solution@data),
      number_of_total_units(x) == nrow(solution@data),
      class(x$data$cost)[1] == class(solution)[1],
      is.numeric(unlist(solution@data)),
      min(unlist(solution@data), na.rm = TRUE) >= 0,
      max(unlist(solution@data), na.rm = TRUE) <= 1)
    # subset planning units with finite cost values
    pos <- x$planning_unit_indices()
    solution <- as.matrix(solution@data)
    if (any(solution[setdiff(seq_len(nrow(solution)), pos), ,
                     drop = FALSE] > 0))
      stop("planning units with NA cost data have non-zero allocations in the ",
           "argument to solution")
    solution <- solution[pos, , drop = FALSE]
    # calculate amount of each feature in each planning unit
    total <- x$feature_abundances_in_total_units()
    held <- vapply(seq_len(x$number_of_zones()),
                   function(i) rowSums(
                     x$data$rij_matrix[[i]] *
                     matrix(solution[, i], ncol = nrow(solution),
                            nrow = nrow(x$data$rij_matrix[[1]]),
                            byrow = TRUE)),
                     numeric(nrow(x$data$rij_matrix[[1]])))
    out <- tibble::tibble(feature = rep(x$feature_names(), x$number_of_zones()),
                          absolute_held = c(held),
                          relative_held = c(held / total))
    if (x$number_of_zones() > 1) {
      out$zone <- rep(x$zone_names(), each = x$number_of_features())
      out <- out[, c(1, 4, 2, 3), drop = FALSE]
    }
    out
}



