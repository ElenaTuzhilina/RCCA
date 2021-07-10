#' @title Canonical Correlation analysis with L2 penalty
#' 
#' @description RCCA function performs Canonical Correlation Analysis with L2 regularization and allows to conduct Canonical Correlation Analysis in high dimensions 
#' For a pair of random vectors 
#' \deqn{x = (x_1, \ldots, x_p)~~and~~y = (y_1, \ldots, y_q)}{x = (x_1, ..., x_p)  and  y = (y_1, ..., y_q)}
#' it seeks for such vectors 
#' \deqn{\alpha = (\alpha_1, \ldots, \alpha_p)~~and~~\beta = (\beta_1, \ldots, \beta_q)}{\alpha = (\alpha_1, ..., \alpha_p)  and  \beta = (\beta_1, ..., \beta_q)} 
#' that satisfy L2 constraints 
#' \deqn{\|\alpha\|\leq t_1~~and~~\|\beta\| \leq t_2}{||\alpha|| <= t_1  and  ||\beta|| <= t_2}
#' and that maximize the correlation 
#' \eqn{cor(u, v)} 
#' between the linear combnations 
#' \deqn{u = \langle x,\alpha\rangle and v = \langle y,\beta\rangle}{u = <x , \alpha> and v = <y , \beta>.} 
#' Here 
#' \eqn{<a , b>} 
#' refers to the inner product between two vectors. 
#' The optimal values for 
#' \eqn{\alpha} and \eqn{\beta} 
#' are called canonical coefficients and the resulting linear combinations 
#' \eqn{u} and \eqn{v} 
#' are called canonical variates.
#' It is actually possible to continue the process and find a sequence of canonical coefficients 
#' \deqn{\alpha^1,\ldot,\alpha^k~~and~~\beta^1,\ldots\beta^k}{\alpha[1], ..., \alpha[k]  and  \beta[1], ..., \beta[k]} 
#' that satisfy the L2 constraints and such that linear combinations 
#' \deqn{u^i = \langle x,\alpha^i\rangle~~and~~v^i = \langle y, \beta^i\rangle}{u[i] = <x , \alpha[i]>  and  v[i] = <y , \beta[i]>}  
#' form two sets 
#' \deqn{\{u^1, \ldots u^k\}~~and~~\{v^1, \ldots v^k\}}{{u[1], ..., u[k]}  and  {v[k], ..., v[k]}} 
#' of independent random variables. The maximmum possible number of such canonical variates is 
#' \eqn{k = min(p, q)}.
#' Note that the above optimization problem is equivalet to maximizing the modified correlation coefficient 
#' \deqn{\frac{cov(\langle x, \alpha\rangle, \langle y, \beta\rangle)}{\sqrt{var(\lamgle x, \alpha\rangle) + \lambda_1\|\alpha\|^2}\sqrt{var(\langle y, \beta\rangle) + \lambda_2\|\beta\|^2}}}{cov(<x , \alpha>, <y , \beta>)  /  ( cov(<x , \alpha>) + \lambda_1  ||\alpha||^2 )^1/2 ( var(<y , \beta>) + \lambda_2 ||\beta||^2 )^1/2,}
#' where \deqn{\lambda_1~~and~~\lambda_2}{\lambda_1  and  \lambda_2} control the resulting sparsity of the canonical coefficients. 

#' 
#' @param X a rectangular \eqn{n x p} matrix containing \eqn{n} observations of random vector \eqn{x}. 
#' @param Y a rectangular \eqn{n x q} matrix containing \eqn{n} observations of random vector \eqn{y}.
#' @param lambda1 a non-negative penalty factor used for regularizing \eqn{X} side coefficients. By default \code{lambda1 = 0}, i.e. no regularization is imposed. Increasing \code{lambda1} incourages sparsity of the resulting canonical coefficients.
#' @param lambda2 a non-negative penalty factor used for regularizing \eqn{Y} side coefficients. By default \code{lambda2 = 0}, i.e. no regularization is imposed. Increasing \code{lambda2} incourages sparsity of the resulting canonical coefficients. 
#' @return A list containing the PCMS problem solution:
#' \itemize{
#'   \item \code{n.comp} -- the number of computed canonical components, i.e. \eqn{k = min(p, q)}.
#'   \item \code{cors} -- the resulting \eqn{k} canonical correlations.
#'   \item \code{mod.cors} -- the resulting \eqn{k} values of modified canonical correlation.
#'   \item \code{x.coefs} -- \eqn{p x k} matrix representing \eqn{k} canonical coefficient vectors \eqn{\alpha[1]}, ..., \eqn{\alpha[k]}.
#'   \item \code{x.vars} -- \eqn{n x k} matrix representing \eqn{k} canonical variates \eqn{u[1]}, ..., \eqn{u[k]}.
#'   \item \code{y.coefs} -- \eqn{q x k} matrix representing \eqn{k} canonical coefficient vectors \eqn{\beta[1]}, ..., \eqn{\beta[k]}.
#'   \item \code{y.vars} -- \eqn{n x k} matrix representing \eqn{k} canonical variates \eqn{v[1]}, ..., \eqn{v[k]}.
#' }
#'
#' @examples
#' data(X)
#' data(Y)
#' #run RCCA 
#' rcca = RCCA(X, Y, lambda1 = 10, lambda2 = 0)
#' #check the modified canonical correlations 
#' plot(1:rcca$n.comp, rcca$mod.cors, pch = 16, xlab = 'component', 'ylab' = 'correlation', ylim = c(0, 1))
#' #check the canonical correlations
#' points(1:rcca$n.comp, rcca$cors, pch = 16, col = 'purple')
#' #compare them to cor(x*alpha, y*beta)
#' points(1:rcca$n.comp, diag(cor(X %*% rcca$x.coefs, Y %*% rcca$y.coefs)), col = 'cyan', pch = 16, cex = 0.7)
#' #check the canonical coefficients for the first canonical variates
#' barplot(rcca$x.coefs[,'can.comp1'], col = 'orange', 'xlab' = 'X feature', ylab = 'value')
#' barplot(rcca$y.coefs[,'can.comp1'], col = 'darkgreen', 'xlab' = 'Y feature', ylab = 'value')
#' 
#' @export RCCA

RCCA = function(X, Y, lambda1 = 0, lambda2 = 0){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n.comp = min(ncol(X), ncol(Y), nrow(X))
  
  #transform
  X.tr = RCCA.tr(X, lambda1)
  Y.tr = RCCA.tr(Y, lambda2)
  if(is.null(X.tr) | is.null(Y.tr)) return(invisible(NULL))
  
  #solve optimization problem
  Cxx = X.tr$cor
  Cyy = Y.tr$cor
  Cxy = cov(X.tr$mat, Y.tr$mat, use = "pairwise")
  sol = fda::geigen(Cxy, Cxx, Cyy)
  names(sol) = c("rho", "alpha", "beta")
  #modified canonical correlation
  rho.mod = sol$rho[1:n.comp]
  names(rho.mod) = paste('can.comp', 1:n.comp, sep = '')
  
  #inverse transform
  X.inv.tr = RCCA.inv.tr(X.tr$mat, sol$alpha, X.tr$tr)
  x.coefs = as.matrix(X.inv.tr$coefs[,1:n.comp])
  rownames(x.coefs) = colnames(X)
  x.vars = X.inv.tr$vars[,1:n.comp]
  rownames(x.vars) = rownames(X)
  Y.inv.tr = RCCA.inv.tr(Y.tr$mat, sol$beta, Y.tr$tr)
  y.coefs = as.matrix(Y.inv.tr$coefs[,1:n.comp])
  rownames(y.coefs) = colnames(Y)
  y.vars = Y.inv.tr$vars[,1:n.comp]
  rownames(y.vars) = rownames(Y)
  
  #canonical correlation
  rho = diag(cor(X.inv.tr$vars,  Y.inv.tr$vars))[1:n.comp]
  names(rho) = paste('can.comp', 1:n.comp, sep = '')
  
  return(list('n.comp' = n.comp, 'cors' = rho, 'mod.cors' = rho.mod, 'x.coefs' = x.coefs, 'x.vars' = x.vars, 'y.coefs' = y.coefs, 'y.vars' = y.vars))
}

RCCA.tr = function(X, lambda){
  p = ncol(X)
  n = nrow(X)
  if(lambda < 0){
    cat('please make', deparse(substitute(lambda)),'>= 0\n')
    return()
  }
  #kernel trick
  V = NULL
  if(p > n){
    if(lambda == 0){
      cat('singularity issue. please impose a penalty on', deparse(substitute(X)), 'side\n')
      return()
    }
    SVD = svd(X)
    X = SVD$u %*% diag(SVD$d)
    V = SVD$v
  }
  C = var(X, use = "pairwise")
  diag(C) = diag(C) + lambda
  return(list('mat' = X, 'tr' = V, 'cor' = C))
}

RCCA.inv.tr = function(X, alpha, V){
  n.comp = ncol(alpha)
  #find canonical variates
  u = X %*% alpha
  colnames(u) = paste('can.comp', 1:n.comp, sep = '')
  #inverse transfrom canonical coefficients 
  if(!is.null(V)) alpha = V %*% alpha
  colnames(alpha) = paste('can.comp', 1:n.comp, sep = '')
  return(list('coefs' = alpha, 'vars' = u))
}