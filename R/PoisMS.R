#' @title Poisson Metric Scaling
#'
#' @description PoisMS function calculates the Poisson Metric Scaling solution for a contact matrix \eqn{C} and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the negative log-likelihood for the Poisson probabilistic model \eqn{C~Pois(\Lambda)} with \deqn{\log(\Lambda) = -D^2(X) + \beta}{log(\Lambda) = -D^2(X) + \beta} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' Here \eqn{D(X)} refers to the matrix of pairwise distances.
#' The solution can be calculated via iterating the second order approximation of the objective (outer PoisMS loop) and applying WPCMS to optimize the obtained quadratic approximation (inner WPCMS loop).
#' The spatial coordiantes of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param C a square symmetric matrix representing a Hi-C contact matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalization should be done via QR decomposition.
#' @param beta0 an initialization for intercept \eqn{beta}. By default \code{beta0 = max(log(C))}.
#' @param Theta0 an initialization for spline basis coefficient matrix \eqn{Theta}. By defaul \code{Theta0 = matrix(rnorm(ncol(H) * 3), ncol(H), 3)}, i.e. a random initialization is considered.
#' @param update_beta If \code{update_beta = TRUE}, then the algorithm finds an optimal intercept value. If \code{update_beta = TRUE}, then the intercept is considered to be fixed and set to \code{beta0}.
#' @param eps_wpcms,eps_poisms positive convergence tolerances for WPCMS inner loop and PoisMS outer loop.
#' @param maxiter,maxepoch integers giving the maximal numbers of iterations for WPCMS inner loop and PoisMS outer loop.
#' @param verbose_wpcms,verbose_poisms logical. If \code{TRUE}, the WPCMS loss after each iteration of inner loop and PoisMS loss after each epoch of outer loop are printed.
#' @return A list containing the PoisMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{beta} -- the resulting intercept value.
#'   \item \code{loss} -- the resulting value of the PoisMS loss.
#'   \item \code{epoch} -- the total number of epochs.
#'   \item \code{iter_total} -- the total number of iterations.
#'   \item \code{plot} -- the list of PoisMS plots: 'loss' corresponds to loss vs. epoch plot, 'intercept' represents beta vs. epoch plot.
#' }
#'
#' @examples
#' data(C)
#'
#' #create spline basis matrix
#' H = splines::bs(1:ncol(C), df = 5)
#' 
#' #orthogonalize H using QR decomposition
#' H = qr.Q(qr(H))
#' 
#' #run PoisMS approach; fixed intercept
#' PoisMS(C, H, beta = 10)$X
#'
#' #run PoisMS approach; optimize intercept
#' PoisMS(C, H)$X
#' 
#' @export PoisMS


PoisMS = function(C, H, beta0 = max(log(C)), Theta0 = matrix(rnorm(ncol(H) * 3), ncol(H), 3), update_beta = TRUE, eps_wpcms = 1e-6, maxiter = 100, verbose_wpcms = FALSE, eps_poisms = 1e-6, maxepoch = 100, verbose_poisms = FALSE){
  #Initialize
  n = nrow(H)
  df = ncol(H)
  beta = beta0
  betas = beta
  Theta = Theta0
  X = H %*% Theta
  X = scale(X, scale = FALSE, center = TRUE)
  loss = loss_PoisMS(X, C, beta)
  losses = loss
  delta = Inf
  deltaX = Inf
  epoch = 0
  iter_total = 0
  if(verbose_poisms) cat("\n--------------------\nEpoch:", epoch, "Beta:", beta, "PoisMS Objective:", loss, "\n--------------------\n")
  #Iterate
  while(deltaX > eps_poisms && epoch < maxepoch){
    epoch = epoch + 1
    #SOA
    soa = SOA(X, C, beta)
    W = soa$W
    diag(W) = 0
    W = W/max(W)
    Z = soa$Z
    #WPCMS
    wpcms = WPCMS(Z, H, W, beta, Theta, FALSE, eps_wpcms, maxiter, verbose_wpcms)
    iter_total = iter_total + wpcms$iter
    #Line search
    find = PoisMS_rate(Theta, wpcms$X, X, beta, loss, C, H)
    #Update solution
    Theta = find$Theta
    X_new = find$X
    X_new = scale(X_new, scale = FALSE, center = TRUE)
    deltaX = vegan::procrustes(X, X_new, scale = TRUE, symmetric = TRUE)$ss
    X = X_new
    #Update alpha and beta
    if(update_beta) beta = update_beta_PoisMS(X, C)
    betas = c(betas, beta)
    #Update loss
    loss_new = loss_PoisMS(X, C, beta)
    delta = abs((loss - loss_new)/loss)
    loss = loss_new
    losses = c(losses, loss)
    if(verbose_poisms) cat("\n--------------------\nEpoch:", epoch, 'Beta:', beta, 'Rate:', find$rate, "WPCMS iters:", wpcms$iter, "PoisMS Objective:", loss, "Delta PoisMS objectives:", delta, "Delta X:", deltaX, "\n--------------------\n")
  }
  data = data.frame('epoch' = 0:epoch, 'objective' = losses)
  plt_obj = ggplot2::ggplot(data, ggplot2::aes(x = epoch, y = objective))+
    ggplot2::geom_line(color = 'red')+
    ggplot2::geom_point(color = 'red', size = 2)+
    ggplot2::ylab('PoisMS objective')+
    ggplot2::xlab('Epoch')+
    ggplot2::ggtitle(paste('Degree of freedom = ', ncol(H), ', Loss value = ', round(loss, 2), '\nNumber of epochs = ', epoch, ', Total number of iterations = ', iter_total))+
    ggplot2::theme(plot.title = ggplot2::element_text(size=9, face="bold", hjust = 0.5), aspect.ratio = 0.5)
  data = data.frame('epoch' = 0:epoch, 'beta' = betas)
  plt_beta = ggplot2::ggplot(data, ggplot2::aes(x = epoch, y = beta))+
    ggplot2::geom_line()+
    ggplot2::geom_point(size = 2)+
    ggplot2::ylab('Beta')+
    ggplot2::xlab('Epoch')+
    ggplot2::ggtitle('PoisMS intercept')+ 
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5), aspect.ratio = 0.5)
  return(list(Theta = Theta, X = X, beta = beta, loss = loss, epoch = epoch, iter_total = iter_total, plot = list('objective' = plt_obj, 'intercept' = plt_beta)))
}


loss_PoisMS = function(X, C, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  logL = -D^2 + beta
  return(mean(exp(logL) - C * logL))
}


update_beta_PoisMS = function(X, C){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  beta = log(sum(C)/sum(exp(-D^2)))
  return(beta)
}

SOA = function(X, C, beta){
  D = as.matrix(dist(X, diag = TRUE, upper = TRUE))
  logL = -D^2 + beta
  L = exp(logL)
  Z = D^2 - (C/L - 1)
  W = L
  return(list(W = W, Z = Z))
}

PoisMS_rate = function(Theta, X_new, X, beta, loss, C, H){
  rate = 1
  S = X %*% t(X)
  S_new = X_new %*% t(X_new)
  S_star = (1 - rate) * S + rate * S_new
  pcms = PCMS(S_star, H)
  while(loss < loss_PoisMS(pcms$X, C, beta) || pcms$rank < 3){
    rate = rate * 0.5
    S_star = (1 - rate) * S + rate * S_new
    pcms = PCMS(S_star, H)
    if(rate < 1e-20) return(list(Theta = Theta, X = X, rate = 0))
  }
  return(list(Theta = pcms$Theta, X = pcms$X, rate = rate))
}