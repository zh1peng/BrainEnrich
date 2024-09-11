#' Generate Null Brain Data
#'
#' This function generates null brain datasets based on permutations provided in perm_id.
#' It rearranges the brain data according to the permutations and outputs the shuffled datasets.
#'
#' @param brain_data A matrix representing brain data, where each row corresponds to a region.
#' @param perm_id A matrix of permutations, where each column represents a permutation and each row corresponds to an index in `brain_data`.
#'
#' @return A matrix of null brain data with the same dimensions as `brain_data` but with permuted rows according to `perm_id`.
#' @export

generate_null_brain_data <- function(brain_data, perm_id) {
  # Check for duplicates in perm_id
  if (any(duplicated(perm_id))) {
    stop("Duplicate entries found in perm_id.")
  }

  if (ncol(brain_data) != 1) {
    stop("This function is designed for group-level analysis and supports one column of brain data. Consider using lapply.")
  }
  if (is.data.frame(brain_data)) {
    brain_data <- as.matrix(brain_data)
  }

  # Define dimensions
  region_n <- nrow(perm_id)
  perm_n <- ncol(perm_id)

  # Check if dimensions match
  if (nrow(brain_data) != region_n) {
    stop("The number of regions in brain_data and perm_id do not match.")
  }

  # Generate null brain data
  null_brain_data <- sapply(1:perm_n, function(idx) {
    brain_data[perm_id[, idx], , drop = FALSE]
  })

  # Set row and column names
  rownames(null_brain_data) <- rownames(brain_data)
  colnames(null_brain_data) <- paste0("null_", 1:perm_n)

  return(null_brain_data)
}





#' Rotate Parcellation
#'
#' Generate a permutation map from a set of cortical regions of interest to itself,
#' while (approximately) preserving contiguity and hemispheric symmetry.
#' The function is based on a rotation of the FreeSurfer projection of coordinates
#' of a set of regions of interest on the sphere.
#' #' This function is modified from the original version available at:
#' https://github.com/frantisekvasa/rotate_parcellation. Please cite relevant papers from this source if you used this function.
#'
#' Modifications include:
#' - Added support for scenarios where only one hemisphere's coordinates are provided.
#' - Improved handling of coordinate dimensions and conditional concatenation of reference and rotation indices.
#' - Included importFrom directives for required functions from `matrixStats` and `clue`.
#' - Ensured the function generates `nrot + 100` permutations, removes duplicates, and returns exactly `nrot` unique permutations.

#' @param coord.l Coordinates of left hemisphere regions on the sphere (array of size n(LH regions) x 3). Can be NULL if only right hemisphere is used.
#' @param coord.r Coordinates of right hemisphere regions on the sphere (array of size n(RH regions) x 3). Can be NULL if only left hemisphere is used.
#' @param nrot Number of rotations (default = 5000).
#' @param method Method to match rotated and unrotated regions; options are 'vasa' (faster, can be suboptimal) or 'hungarian' (default, slower, optimal).
#' @param seed Seed for reproducibility.
#' @return Array of permutations, from set of regions to itself (array of size n(total regions) x nrot).
#' @importFrom matrixStats rowMins
#' @importFrom clue solve_LSAP
#' @importFrom stats rnorm
#' @export

rotate_parcellation <- function(coord.l = NULL,
                                coord.r = NULL, nrot = 5000,
                                method = c("hungarian", "vasa"),
                                seed = NULL) {
  message("If you use this function to create perm_ids, consider citing the original source at https://github.com/frantisekvasa/rotate_parcellation and relevant papers.")
  method <- match.arg(method)

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }


  # Check that at least one set of coordinates is provided
  if (is.null(coord.l) && is.null(coord.r)) {
    stop("At least one of coord.l or coord.r must be provided.")
  }

  if (!is.null(coord.l) & is.data.frame(coord.l)) {
    coord.l <- as.matrix(coord.l)
  }

  if (!is.null(coord.r) & is.data.frame(coord.r)) {
    coord.r <- as.matrix(coord.r)
  }



  # Check that coordinate dimensions are correct
  if (!is.null(coord.l) && dim(coord.l)[2] != 3) {
    if (dim(coord.l)[1] == 3) {
      coord.l <- t(coord.l)
    } else {
      stop("Left hemisphere coordinates must have dimensions [n(LH regions) x 3].")
    }
  }
  if (!is.null(coord.r) && dim(coord.r)[2] != 3) {
    if (dim(coord.r)[1] == 3) {
      coord.r <- t(coord.r)
    } else {
      stop("Right hemisphere coordinates must have dimensions [n(RH regions) x 3].")
    }
  }

  nroi.l <- if (!is.null(coord.l)) dim(coord.l)[1] else 0
  nroi.r <- if (!is.null(coord.r)) dim(coord.r)[1] else 0
  nroi <- nroi.l + nroi.r

  perm_id <- array(0, dim = c(nroi, round(nrot * 1.2)))
  r <- 0
  c <- 0

  I1 <- diag(3)
  I1[1, 1] <- -1

  while (r < round(nrot * 1.2)) {
    A <- matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec <- qr(A)
    TL <- qr.Q(qrdec)
    temp <- qr.R(qrdec)
    TL <- TL %*% diag(sign(diag(temp)))
    if (det(TL) < 0) {
      TL[, 1] <- -TL[, 1]
    }
    TR <- I1 %*% TL %*% I1

    coord.l.rot <- if (!is.null(coord.l)) coord.l %*% TL else NULL
    coord.r.rot <- if (!is.null(coord.r)) coord.r %*% TR else NULL


    if (!is.null(coord.l)) {
      dist.l <- array(0, dim = c(nroi.l, nroi.l))
      for (i in 1:nroi.l) { # left
        for (j in 1:nroi.l) {
          dist.l[i, j] <- sqrt(sum((coord.l[i, ] - coord.l.rot[j, ])^2))
        }
      }
    } else {
      dist.l <- NULL
    }

    if (!is.null(coord.r)) {
      dist.r <- array(0, dim = c(nroi.r, nroi.r))
      for (i in 1:nroi.r) { # right
        for (j in 1:nroi.r) {
          dist.r[i, j] <- sqrt(sum((coord.r[i, ] - coord.r.rot[j, ])^2))
        }
      }
    } else {
      dist.r <- NULL
    }

    if (method == "vasa") {
      if (!is.null(coord.l)) {
        temp.dist.l <- dist.l
        rot.l <- ref.l <- numeric(nroi.l)
        for (i in seq_len(nroi.l)) {
          ref.ix <- which(rowMins(temp.dist.l, na.rm = TRUE) == max(rowMins(temp.dist.l, na.rm = TRUE), na.rm = TRUE))
          rot.ix <- which(temp.dist.l[ref.ix, ] == min(temp.dist.l[ref.ix, ], na.rm = TRUE))
          ref.l[i] <- ref.ix
          rot.l[i] <- rot.ix
          temp.dist.l[, rot.ix] <- NA
          temp.dist.l[ref.ix, ] <- 0
        }
      }

      if (!is.null(coord.r)) {
        temp.dist.r <- dist.r
        rot.r <- ref.r <- numeric(nroi.r)
        for (i in seq_len(nroi.r)) {
          ref.ix <- which(rowMins(temp.dist.r, na.rm = TRUE) == max(rowMins(temp.dist.r, na.rm = TRUE), na.rm = TRUE))
          rot.ix <- which(temp.dist.r[ref.ix, ] == min(temp.dist.r[ref.ix, ], na.rm = TRUE))
          ref.r[i] <- ref.ix
          rot.r[i] <- rot.ix
          temp.dist.r[, rot.ix] <- NA
          temp.dist.r[ref.ix, ] <- 0
        }
      }
    } else if (method == "hungarian") {
      if (!is.null(coord.l)) {
        rot.l <- as.vector(solve_LSAP(dist.l, maximum = FALSE))
        ref.l <- seq_len(nroi.l)
      }

      if (!is.null(coord.r)) {
        rot.r <- as.vector(solve_LSAP(dist.r, maximum = FALSE))
        ref.r <- seq_len(nroi.r)
      }
    }

    # Correctly concatenate vectors without including NULL values
    ref.lr <- c()
    rot.lr <- c()
    if (!is.null(coord.l)) {
      ref.lr <- c(ref.lr, ref.l)
      rot.lr <- c(rot.lr, rot.l)
    }
    if (!is.null(coord.r)) {
      ref.lr <- c(ref.lr, nroi.l + ref.r)
      rot.lr <- c(rot.lr, nroi.l + rot.r)
    }

    b <- sort(ref.lr, index.return = TRUE)
    ref.lr.sort <- ref.lr[b$ix]
    rot.lr.sort <- rot.lr[b$ix]

    if (!all(sort(rot.lr.sort, decreasing = FALSE) == seq_len(nroi))) {
      browser("permutation error")
    }

    if (!all(rot.lr.sort == seq_len(nroi))) {
      r <- r + 1
      perm_id[, r] <- rot.lr.sort
    } else {
      c <- c + 1
      message(sprintf("map to itself n. %d", c))
    }

    if (r %% 500 == 0) message(sprintf("permutation %d of %d", r, nrot * 1.2))
  }
  # Remove duplicate columns and ensure only nrot unique permutations are returned
  perm_id <- perm_id[, !duplicated(t(perm_id))]
  if (ncol(perm_id) > nrot) {
    perm_id <- perm_id[, 1:nrot]
  }
  return(perm_id)
}
