#' @title Weighted Principal Curve Metric Scaling
#'
#' @description WPCMS function calculates the Weighted Principal Curve Metric Scaling solution for a contact matrix \eqn{C} and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the weighted Frobenius norm \deqn{ \|W * (C-X X^T) \|^2_F }{|| W x (C - XX') ||} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' The spatial coordiantes of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param C a square symmetric matrix representing a Hi-C contact matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalized via QR decomposition.
#' @param W a weights matrix, same dimension as \code{C}. If \code{W = NULL}, the Principal Curve Metric Scaling solution is calculated.
#' @param X0 an initialization for \eqn{X} of the WPCMS algorithm. If \code{X0 = NULL}, the PCMS warm start is considered.
#' @param eps a positive convergence tolerance.
#' @param maxiter an integer giving the maximal number of iterations.
#' @param verbose logical. If \code{TRUE}, the WPCMS loss after each iteration is printed.
#' @return A list containing the WPCMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{loss} -- the resulting value of the WPCMS loss.
#'   \item \code{iter} -- the total number of iterations.
#' }
#'
#' @examples
#' data(C)
#' data(H)
#' 
#' ##W = NULL, the PCMS solution is calculated
#' WPCMS(C, H)$X
#' 
#' ##The weighted PCMS solution is calculated
#' W = matrix(runif(length(C)), dim(C))
#' WPCMS(C, H, W)$X
#'
#' @export WPCMS


WPCMS = function(C, H, W = NULL, X0 = NULL, eps = 1e-6, maxiter = 100, verbose = FALSE){
  #Orthogonalize H
  QR = qr(H)
  H = qr.Q(QR)

  #Initialize
  pcms = PCMS(C, H)
  Theta = pcms$Theta
  X = pcms$X
  loss = l_WPCMS(X, C, W)
  delta = Inf
  iter = 0

  if(is.null(W)){
    return(list(Theta = Theta, X = X, loss = loss, iter = iter))
  }

  if(!is.null(X0)){
    X = X0
    loss = l_WPCMS(X, C, W)
  }

  if(verbose) cat("Iter:", iter, "WPCMS objective:", loss, "Delta WPCMS objectives:", delta, "\n")

  #Iterate
  while(delta > eps && iter < maxiter){
    iter = iter + 1

    #Line search
    find = WPCMS_rate(X, loss, C, H, W)

    #Update solution
    Theta = find$Theta
    X = find$X
    delta = abs((loss - find$loss)/loss)
    loss = find$loss

    if(verbose) cat("Iter:", iter, 'Rate:', find$rate, "WPCMS Objective:", loss, "Delta WPCMS objectives:", delta, "\n")
  }
  return(list(Theta = Theta, X = X, loss = loss, iter = iter))
}

l_WPCMS = function(X, C, W = NULL){
  C_hat = X%*%t(X)
  if(is.null(W)) return(mean((C - C_hat)^2))
  else return(mean(W * (C - C_hat)^2))
}

PCMS = function(C, H){
  ED = eigen(t(H)%*%C%*%H, symmetric = TRUE)
  U3 = ED$vectors[,1:3]
  d3 = ED$values[1:3]
  d3 = pmax(d3, 0)
  rank = sum(d3 > 1e-10)
  Theta = U3%*%diag(sqrt(d3))
  X = H%*%Theta
  return(list(Theta = Theta, X = X, loss = l_WPCMS(X, C), rank = rank))
}

WPCMS_rate = function(X, loss, C, H, W){
  rate = 1
  C_hat = X%*%t(X)

  C_star = rate * W * C + (1 - rate * W) * C_hat
  pcms = PCMS(C_star, H)

  while(loss < l_WPCMS(pcms$X, C, W)|| pcms$rank < 3){
    rate = rate * 0.5
    C_star = rate * W * C + (1 - rate * W) * C_hat
    pcms = PCMS(C_star, H)
  }

  return(list(Theta = pcms$Theta, X = pcms$X, loss = l_WPCMS(pcms$X, C, W), rate = rate))
}
