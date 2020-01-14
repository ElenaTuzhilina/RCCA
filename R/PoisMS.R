#' @title Poisson Metric Scaling
#'
#' @description PoisMS function calculates the Poisson Metric Scaling solution for a contact matrix \eqn{C} and a spline basis matrix \eqn{H}.
#' The optimal solution is found via minimizing the negative log-likelihood for the Poisson probabilistic model \eqn{C~Pois(\Lambda)} with \deqn{\log(\Lambda) = \alpha XX^T + \beta}{log(\Lambda) = \alpha XX' + \beta} w.r.t. \eqn{\Theta} subject to the smooth curve constraint \eqn{X = H\Theta}.
#' The solution can be calculated via iterating the second order approximation of the objective (outer PoisMS loop) and applying WPCMS to optimize the obtained quadratic approximation (inner WPCMS loop).
#' The spatial coordiantes of the resulting reconstruction are presented in \eqn{X}.
#'
#' @param C a square symmetric matrix representing a Hi-C contact matrix.
#' @param H a spline basis matrix. By default assumed to have orthogonal columns. If not, orthogonalized via QR decomposition.
#' @param alpha,beta scaling and centering constants for Poisson model.
#' @param X0 an initialization for \eqn{X} of the PoisMS algorithm. If \code{X0 = NULL}, the PCMS warm start is considered.
#' @param eps_wpcms,eps_poisms positive convergence tolerances for WPCMS inner loop and PoisMS outer loop.
#' @param maxiter,maxepoch integers giving the maximal numbers of iterations for WPCMS inner loop and PoisMS outer loop.
#' @param verbose_wpcms,verbose_poisms logical. If \code{TRUE}, the WPCMS loss after each iteration of inner loop and PoisMS loss after each epoch of outer loop are printed.
#' @return A list containing the PoisMS problem solution:
#' \itemize{
#'   \item \code{Theta} -- the matrix of spline parameters.
#'   \item \code{X} -- the resulting conformation reconstruction.
#'   \item \code{loss} -- the resulting value of the PoisMS loss.
#'   \item \code{epoch} -- the total number of epochs.
#'   \item \code{plot} -- the convergence plot.
#' }
#'
#' @examples
#' data(C)
#' data(H)
#' ##PoisMS solution
#' alpha = 1
#' beta = log(mean(C))
#' PoisMS(C, H, alpha, beta)$X
#'
#' @export PoisMS

PoisMS = function(C, H, alpha = 1, beta, X0 = NULL, eps_wpcms = 1e-6, maxiter = 100, verbose_wpcms = FALSE, eps_poisms = 1e-6, maxepoch = 100, verbose_poisms = FALSE){
  #Orthogonalize H
  QR = qr(H)
  H = qr.Q(QR)

  #Initialize
  pcms = PCMS((log(C + 1) - beta)/alpha, H)
  Theta = pcms$Theta
  X = pcms$X
  if(!is.null(X0)) X = X0
  loss = l_PoisMS(X, C, alpha, beta)
  delta = Inf
  epoch = 0
  iter_total = 0
  losses = loss

  if(verbose_poisms) cat("--------------------\nEpoch:", epoch, "PoisMS Objective:", loss, "Delta PoisMS objectives:", delta, "\n--------------------\n\n")

  while(delta > eps_poisms && epoch < maxepoch){
    epoch = epoch + 1

    #SOA
    soa = SOA(X, C, alpha, beta)
    W = soa$W
    Z = soa$Z

    #WPCMS
    wpcms = WPCMS(Z, H, W, X, eps_wpcms, maxiter, verbose_wpcms)
    iter_total = iter_total + wpcms$iter

    #Line search
    find = PoisMS_rate(wpcms$X, X, loss, C, H, alpha, beta)

    #Update solution
    Theta = find$Theta
    X = find$X
    delta = abs((loss - find$loss)/loss)
    loss = find$loss
    losses = c(losses, loss)

    if(verbose_poisms) cat("--------------------\nEpoch:", epoch, 'Rate:', find$rate, "WPCMS iters:", wpcms$iter, "PoisMS Objective:", loss, "Delta PoisMS objectives:", delta, "\n--------------------\n\n")
  }

  result = data.frame('epoch' = 0:epoch, 'loss' = losses)
  p = ggplot2::ggplot(data = result, ggplot2::aes(x = epoch, y = loss))+
    ggplot2::geom_line(colour = 'red')+
    ggplot2::geom_point(colour = 'red', size = 2)+
    ggplot2::ylab('Poisson objective')+
    ggplot2::xlab('Epoch')+
    ggplot2::ggtitle(paste('Degree of freedom = ', ncol(H), ', Loss value = ', round(loss, 2), '\nNumber of epochs = ', epoch, ', Total number of iterations = ', iter_total))+
    ggplot2::theme(plot.title = ggplot2::element_text(size=9, face="bold", hjust = 0.5), aspect.ratio = 0.5)

  return(list(Theta = Theta, X = X, loss = loss, epoch = epoch, plot = p))
}

l_PoisMS = function(X, C, alpha, beta){
  C_hat = X%*%t(X)
  logL = alpha * C_hat + beta
  return(mean(exp(logL) - C * logL))
}

SOA = function(X, C, alpha, beta){
  C_hat = X%*%t(X)
  logL = alpha * C_hat + beta
  Z = C_hat + 1/alpha * (C * exp(-logL) - 1)
  W = exp(logL - max(logL))
  return(list(W = W, Z = Z))
}

PoisMS_rate = function(X_new, X, loss, C, H, alpha, beta){
  rate = 1
  C_hat = X%*%t(X)
  C_hat_new = X_new%*%t(X_new)

  C_star = (1 - rate) * C_hat + rate * C_hat_new
  pcms = PCMS(C_star, H)

  while(loss < l_PoisMS(pcms$X, C, alpha, beta)|| pcms$rank < 3){
    rate = rate * 0.5
    C_star = (1 - rate) * C_hat + rate * C_hat_new
    pcms = PCMS(C_star, H)
  }

  return(list(Theta = pcms$Theta, X = pcms$X, loss = l_PoisMS(pcms$X, C, alpha, beta), rate = rate))
}
