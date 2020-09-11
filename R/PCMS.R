#' @title Principal Curve Metric Scaling
#'
#' @description PCMS function calculates the Principal Curve Metric Scaling solution for some matrix \eqn{Z} (e.g. representing a Hi-C contact matrix after some proper transformation applied) and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the Frobenius norm \deqn{\|Z - S(X)\|^2_F}{||Z - S(X)||} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}. 
#' Here \deqn{S(X) = XX^T}{S(X) = XX'} refers to the matrix of inner products.
#' The spatial coordiantes of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param Z a square symmetric matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalization should be done via QR decomposition.
#' @return A list containing the PCMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{loss} -- the resulting value of the PCMS loss.
#' }
#'
#' @examples
#' data(C)
#' 
#' #transform contact counts to similarities
#' D = 1/(C+1)
#' Z = -D^2/2
#' Z = scale(Z, scale = FALSE, center = TRUE)
#' Z = t(scale(t(Z), scale = FALSE, center = TRUE))
#' 
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #compute the PCMS solution
#' PCMS(Z, H)$X
#'
#' @export PCMS

PCMS = function(Z, H){
  ED = eigen(t(H) %*% Z %*% H, symmetric = TRUE)
  U3 = ED$vectors[,1:3]
  d3 = ED$values[1:3]
  rank = sum(d3 > 1e-12)
  d3 = pmax(d3, 0)
  Theta = U3%*%diag(sqrt(d3))
  X = H%*%Theta
  return(list(Theta = Theta, X = X, loss = loss_PCMS(X, Z), rank = rank))
}

loss_PCMS = function(X, Z){
  S = X%*%t(X)
  return(mean((Z - S)^2))
}