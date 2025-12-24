#' @title Canonical Polyadic Decomposition
#' @description CP decomposition of a tensor using Alternating Least Squares
#' @importFrom methods is
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar tail txtProgressBar
#' @param tnsr Tensor with K modes
#' @param num_components the number of rank-1 K-Tensors to use in approximation
#' @param max_iter maximum number of iterations if error stays above tol 
#' @param tol relative Frobenius norm error tolerance
#' @return a list containing decomposition results
cpDecomposition <- function(tnsr, num_components = NULL, max_iter = 25, tol = 1e-5){
  
  # Helper functions
  kronecker_list <- function(L){
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- kronecker(retmat, L[[i]])
    }
    retmat
  }
  
  fnorm <- function(tnsr){
    arr <- tnsr$data
    sqrt(sum(arr * arr))
  }
  
  rs_unfold <- function(tnsr, m = NULL){
    if(is.null(m)) stop("mode m must be specified")
    num_modes <- tnsr$num_modes
    rs <- m
    cs <- (1:num_modes)[-m]
    unfold(tnsr, row_idx = rs, col_idx = cs)
  }
  
  unfold <- function(tnsr, row_idx = NULL, col_idx = NULL){
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs) || is.null(cs)) stop("row and column indices must be specified")
    num_modes <- tnsr$num_modes
    if(length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
    perm <- c(rs, cs)
    modes <- tnsr$modes
    mat <- tnsr$data
    new_modes <- c(prod(modes[rs]), prod(modes[cs]))
    mat <- aperm(mat, perm)
    dim(mat) <- new_modes
    as.tensor(mat)
  }
  
  hadamard_list <- function(L){
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- retmat * L[[i]]
    }
    retmat
  }
  
  khatri_rao_list <- function(L, reverse = FALSE){
    ncols <- unique(sapply(L, ncol))
    if(length(ncols) != 1) stop("All matrices must have same number of columns")
    ncols <- ncols[1]
    nrows <- sapply(L, nrow)
    retmat <- matrix(0, nrow = prod(nrows), ncol = ncols)
    if(reverse) L <- rev(L)
    for(j in 1:ncols){
      Lj <- lapply(L, function(x) x[, j])
      retmat[, j] <- kronecker_list(Lj)
    }
    retmat
  }
  
  superdiagonal_tensor <- function(num_modes, len, elements = 1L){
    modes <- rep(len, num_modes)
    arr <- array(0, dim = modes)
    if(length(elements) == 1) elements <- rep(elements, len)
    for(i in 1:len){
      idx <- rep(i, num_modes)
      arr[matrix(idx, nrow = 1)] <- elements[i]
    }
    as.tensor(arr)
  }
  
  ttl <- function(tnsr, list_mat, ms = NULL){
    if(is.null(ms) || !is.vector(ms)) stop("m modes must be specified as a vector")
    if(length(ms) != length(list_mat)) stop("m modes length does not match list_mat length")
    for(i in 1:length(list_mat)){
      mat <- list_mat[[i]]
      m <- ms[i]
      mat_dims <- dim(mat)
      modes_in <- tnsr$modes
      modes_out <- modes_in
      modes_out[m] <- mat_dims[1]
      tnsr_m <- rs_unfold(tnsr, m = m)$data
      retarr_m <- mat %*% tnsr_m
      tnsr <- rs_fold(retarr_m, m = m, modes = modes_out)
    }	
    tnsr	
  }
  
  rs_fold <- function(mat, m = NULL, modes = NULL){
    if(is.null(m)) stop("mode m must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    num_modes <- length(modes)
    rs <- m
    cs <- (1:num_modes)[-m]
    fold(mat, row_idx = rs, col_idx = cs, modes = modes)
  }
  
  fold <- function(mat, row_idx = NULL, col_idx = NULL, modes = NULL){
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs) || is.null(cs)) stop("row space and col space indices must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    if(!is.list(mat) && !is.matrix(mat)) stop("mat must be of class 'matrix'")
    if(is.list(mat)){
      mat <- mat$data			
    }
    num_modes <- length(modes)
    mat_modes <- dim(mat)
    iperm <- match(1:num_modes, c(rs, cs))
    as.tensor(aperm(array(mat, dim = c(modes[rs], modes[cs])), iperm))
  }
  
  # Main CP decomposition algorithm
  if(is.null(num_components)) stop("num_components must be specified")
  
  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  U_list <- vector("list", num_modes)
  unfolded_mat <- vector("list", num_modes)
  tnsr_norm <- fnorm(tnsr)
  
  # Initialize with random values
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr, m = m)$data
    U_list[[m]] <- matrix(rnorm(modes[m] * num_components), nrow = modes[m], ncol = num_components)
  }
  
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(as.tensor(est$data - tnsr$data))
    fnorm_resid[curr_iter] <<- curr_resid
    if(curr_iter == 1) return(FALSE)
    if(abs(curr_resid - fnorm_resid[curr_iter-1]) / tnsr_norm < tol) return(TRUE)
    return(FALSE)
  }	
  
  pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  
  norm_vec <- function(vec){
    sqrt(sum(vec^2))
  }
  
  # ALS iterations
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb, curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m], function(x) {t(x) %*% x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]] %*% khatri_rao_list(U_list[-m], reverse = TRUE) %*% V_inv
      lambdas <- apply(tmp, 2, norm_vec)
      U_list[[m]] <- sweep(tmp, 2, lambdas, "/")	
      Z <- superdiagonal_tensor(num_modes = num_modes, len = num_components, elements = lambdas)
      est <- ttl(Z, U_list, ms = 1:num_modes)
    }
    
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb, max_iter)
    } else {
      curr_iter <- curr_iter + 1
    }
  }
  
  if(!converged){
    setTxtProgressBar(pb, max_iter)
  }
  close(pb)
  
  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100
  
  invisible(list(
    lambdas = lambdas, 
    U = U_list, 
    conv = converged, 
    est = est, 
    norm_percent = norm_percent, 
    fnorm_resid = tail(fnorm_resid, 1),
    all_resids = fnorm_resid
  ))
}

