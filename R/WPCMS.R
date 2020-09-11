#' @title Weighted Principal Curve Metric Scaling
#'
#' @description WPCMS function calculates the Weighted Principal Curve Metric Scaling solution for a contact matrix \eqn{Z} (e.g. representing a Hi-C contact matrix after some proper transformation applied) and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the weighted Frobenius norm \deqn{ \|W * (Z - D^2(X) + \beta) \|^2_F }{|| W x (Z - D^2(X) + \beta) ||} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' Here \eqn{D(X)} refers to the matrix of pairwise distances.
#' The spatial coordiantes of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param Z a square symmetric matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalization should be done via QR decomposition.
#' @param beta an intercept. Default option \code{beta = NULL}, i.e. the algorithm finds an optimal intercept value. Alternatively, fixed value \code{beta = const} can be considered.
#' @param W a weights matrix, same dimension as \code{Z}. By default equal weights \code{W = 1} are assumed.
#' @param Theta0 an initialization for spline basis coefficient matrix \eqn{Theta}. If \code{Theta0 = NULL}, a random initialization is considered.
#' @param eps a positive convergence tolerance.
#' @param maxiter an integer giving the maximal number of iterations.
#' @param verbose logical. If \code{TRUE}, the WPCMS loss after each iteration is printed.
#' @return A list containing the WPCMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{beta} -- the resulting intercept value.
#'   \item \code{loss} -- the resulting value of the WPCMS loss.
#'   \item \code{iter} -- the total number of iterations.
#'   \item \code{plot} -- the list of WPCMS plots: 'loss' corresponds to loss vs. iteration plot, 'intercept' represents beta vs. iteration plot.
#' }
#'
#' @examples
#' data(C)
#' 
#' #transform contact counts to distances
#' Z = 1/(C+1)
#' 
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #run WPCMS with equal weights; optimize intercept 
#' WPCMS(Z, H)$X
#' 
#' #run WPCMS with random weights; fixed intercept
#' W = matrix(runif(length(Z)), dim(Z))
#' WPCMS(Z, H, W, beta = 1)$X
#'
#' @export WPCMS


WPCMS = function(Z, H, W = matrix(1, nrow(Z), ncol(Z)), beta = NULL, Theta0 = NULL, eps = 1e-6, maxiter = 100, verbose = FALSE){
  #Initialize
  n = nrow(H)
  df = ncol(H)
  if(is.null(Theta0)){
    Theta0 = matrix(rnorm(df * 3), df, 3)
  } 
  Theta = Theta0
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  if(is.null(beta)){
    update = TRUE
    beta = update_beta_WPCMS(X, Z, W)
  } else {
    update = FALSE
  } 
  betas = beta
  loss = loss_WPCMS(X, Z, W, beta)
  losses = loss
  delta = Inf
  deltaX = Inf
  iter = 0
  if(verbose) cat("Iter:", iter, "Beta:", beta, "WPCMS objective:", loss, "\n")
  #Iterate
  while(delta > eps && iter < maxiter){
    iter = iter + 1
    #Line search
    find = WPCMS_rate(Theta, X, loss, Z + beta, H, W)
    #Update solution
    Theta = find$Theta
    X_new = find$X
    X_new = scale(X_new, scale = FALSE, center = TRUE)
    deltaX = vegan::procrustes(X, X_new, scale = FALSE, symmetric = TRUE)$ss
    X = X_new
    #Update beta
    if(update) beta = update_beta_WPCMS(X, Z, W)
    betas = c(betas, beta)
    #Update loss
    loss_new = loss_WPCMS(X, Z, W, beta)
    delta = abs((loss - loss_new)/loss)
    loss = loss_new
    losses = c(losses, loss)
    if(verbose) cat("Iter:", iter, 'Beta:', beta, 'Rate:', find$rate, "WPCMS Objective:", loss, "Delta WPCMS objectives:", delta, 'Delta X:', deltaX, "\n")
  }
  data = data.frame('iter' = 0:iter, 'objective' = losses)
  plt_obj = ggplot2::ggplot(data, ggplot2::aes(x = iter, y = objective))+
    ggplot2::geom_line(color = 'red')+
    ggplot2::geom_point(color = 'red', size = 2)+
    ggplot2::ylab('WPCMS objective')+
    ggplot2::xlab('Iter')+
    ggplot2::ggtitle(paste('Degree of freedom = ', df, ', Loss value = ', round(loss, 2), ', Number of iterations = ', iter))+
    ggplot2::theme(plot.title = ggplot2::element_text(size=9, face="bold", hjust = 0.5), aspect.ratio = 0.5)
  data = data.frame('iter' = 0:iter, 'beta' = betas)
  plt_beta = ggplot2::ggplot(data, ggplot2::aes(x = iter, y = beta))+
    ggplot2::geom_line()+
    ggplot2::geom_point(size = 2)+
    ggplot2::ylab('Beta')+
    ggplot2::xlab('Iter')+
    ggplot2::ggtitle('WPCMS intercept')+ 
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5), aspect.ratio = 0.5)
  return(list(Theta = Theta, X = X, beta = beta, loss = loss, iter = iter, plot = list('objective' = plt_obj, 'intercept' = plt_beta)))
}

loss_WPCMS = function(X, Z, W, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  return(mean(W * (Z - D^2 + beta)^2))
}

update_beta_WPCMS = function(X, Z, W){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  beta = -sum(W * (Z - D^2))/sum(W)
  return(beta)
}

gradient_step = function(S, G, rate){
  G_plus = diag(colSums(G))
  return(S - rate * (G - G_plus))
}

projection_step = function(S, H){
  pcms = PCMS(S, H)
  return(list(Theta = pcms$Theta, X = pcms$X, rank = pcms$rank))
}

WPCMS_rate = function(Theta, X, loss, Z, H, W){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  S = X %*% t(X)
  G = W * (Z - D^2)
  rate = 1
  pgd = projection_step(gradient_step(S, G, rate), H)
  while(loss < loss_WPCMS(pgd$X, Z, W, 0) || pgd$rank < 3){
    rate = rate * 0.5
    pgd = projection_step(gradient_step(S, G, rate), H)
    if(rate < 1e-20) return(list(Theta = Theta, X = X, rate = 0))
  }
  return(list(Theta = pgd$Theta, X = pgd$X, rate = rate))
}

