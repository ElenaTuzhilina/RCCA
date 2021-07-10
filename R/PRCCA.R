#' @title Canonical Correlation analysis with partial L2 penalty
#' 
#' @description PRCCA function performs Canonical Correlation Analysis with partial L2 regularization and allows to conduct Canonical Correlation Analysis in high dimensions. 
#' It imposes L2 penalty only on a subset of \eqn{\alpha} and \eqn{\beta} coefficients. Specifically, if
#' \deqn{x = (x_1, \ldots, x_p)~~and~~y = (y_1, \ldots, y_q)}{x = (x_1, ..., x_p)  and  y = (y_1, ..., y_q)} 
#' are random vectors and 
#' \deqn{I = \{i_1, ..., i_m\} \subset \{1, ..., p\}~~and~~J = {j_1, ..., i_r} \subset \{1, ..., q\}}{I = {i_1, ..., i_m}  is a subset of  {1, ..., p}  and  J = {j_1, ..., i_r}  is a subset of  {1, ..., q}} 
#' then
#' PRCCA seeks for such vectors 
#' \deqn{\alpha = (\alpha_1, \ldots, \alpha_p)~~and~~\beta = (\beta_1, \ldots, \beta_q)}{\alpha = (\alpha_1, ..., \alpha_p)  and  \beta = (\beta_1, ..., \beta_q)} 
#' that satisfy partial L2 constraints 
#' \deqn{\|\alpha_I\|\leq t_1~~and~~\|\beta_J\| \leq t_2}{||\alpha_I|| <= t_1  and  ||\beta_J|| <= t_2}
#' and that maximize the correlation 
#' \eqn{cor(u, v)} 
#' between the linear combnations 
#' \deqn{u = \langle x,\alpha\rangle and v = \langle y,\beta\rangle}{u = <x , \alpha> and v = <y , \beta>.} 
#' Here 
#' \eqn{<a , b>} 
#' refers to the inner product between two vectors and 
#' \deqn{\alpha_I~~and~~\beta_J}{\alpha_I  and  \beta_J} 
#' are corresponding subvectors of \eqn{\alpha} and \eqn{\beta} with indices belonging to \eqn{I} and \eqn{J}, respectively.
#' Again, the above optimization problem is equivalet to maximizing the modified correlation coefficient 
#' \deqn{\frac{cov(\langle x, \alpha\rangle, \langle y, \beta\rangle)}{\sqrt{var(\lamgle x, \alpha\rangle) + \lambda_1\|\alpha_I\|^2}\sqrt{var(\langle y, \beta\rangle) + \lambda_2\|\beta_J\|^2}}}{cov(<x , \alpha>, <y , \beta>)  /  ( cov(<x , \alpha>) + \lambda_1  ||\alpha_I||^2 )^1/2 ( var(<y , \beta>) + \lambda_2 ||\beta_J||^2 )^1/2,}
#' where \deqn{\lambda_1~~and~~\lambda_2}{\lambda_1  and  \lambda_2} control the resulting sparsity of the canonical coefficients within \deqn{\alpha_I~~and~~\beta_J}{\alpha_I  and  \beta_J} parts of the coefficient vectors. 
#' 
#' 
#' @param X a rectangular \eqn{n x p} matrix containing \eqn{n} observations of random vector \eqn{x}. 
#' @param Y a rectangular \eqn{n x q} matrix containing \eqn{n} observations of random vector \eqn{y}.
#' @param index1 a subset of indices the penalty is imposed on while regularizing the \eqn{X} side. By default \code{index1 = 1:ncol(X)}, i.e. we include all \eqn{\alpha} coefficients in the relularization term.
#' @param index2 a subset of indices the penalty is imposed on while regularizing the \eqn{Y} side. By default \code{index2 = 1:ncol(Y)}, i.e. we include all \eqn{\beta} coefficients in the relularization term.
#' @param lambda1 a non-negative penalty factor used for regularizing \eqn{X} side coefficients \eqn{\alpha}. By default \code{lambda1 = 0}, i.e. no regularization is imposed. Increasing \code{lambda1} incourages sparsity of the resulting canonical coefficients.
#' @param lambda2 a non-negative penalty factor used for regularizing \eqn{Y} side coefficients \eqn{\beta}. By default \code{lambda2 = 0}, i.e. no regularization is imposed. Increasing \code{lambda2} incourages sparsity of the resulting canonical coefficients. 
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
#' prcca = PRCCA(X, Y, lambda1 = 100, lambda2 = 0, index1 = 1:(ncol(X) - 10))
#' #check the modified canonical correlations 
#' plot(1:prcca$n.comp, prcca$mod.cors, pch = 16, xlab = 'component', 'ylab' = 'correlation', ylim = c(0, 1))
#' #check the canonical correlations
#' points(1:prcca$n.comp, prcca$cors, pch = 16, col = 'purple')
#' #compare them to cor(x*alpha, y*beta)
#' points(1:prcca$n.comp, diag(cor(X %*% prcca$x.coefs, Y %*% prcca$y.coefs)), col = 'cyan', pch = 16, cex = 0.7)
#' #check the canonical coefficients for the first canonical variates
#' barplot(prcca$x.coefs[,'can.comp1'], col = 'orange', 'xlab' = 'X feature', ylab = 'value')
#' barplot(prcca$y.coefs[,'can.comp1'], col = 'darkgreen', 'xlab' = 'Y feature', ylab = 'value')
#' 
#' @export PRCCA


PRCCA = function(X, Y, index1 = 1:ncol(X), index2 = 1:ncol(Y), lambda1 = 0, lambda2 = 0){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n.comp = min(ncol(X), ncol(Y), nrow(X))
  
  #transform
  X.tr = PRCCA.tr(X, index1, lambda1)
  Y.tr = PRCCA.tr(Y, index2, lambda2)
  if(is.null(X.tr) | is.null(Y.tr)) return(invisible(NULL))
  
  #solve optimization problem
  Cxx = X.tr$cor
  Cyy = Y.tr$cor
  Cxy = cov(X.tr$mat, Y.tr$mat, use = "pairwise")
  sol = fda::geigen(Cxy, Cxx, Cyy)
  names(sol) = c("rho", "alpha", "beta")
  #find modified correlation
  rho.mod = sol$rho[1:n.comp]
  names(rho.mod) = paste('can.comp', 1:n.comp, sep = '')
  
  #inverse transform
  X.inv.tr = PRCCA.inv.tr(X.tr$mat, X.tr$p.pen, sol$alpha, X.tr$tr)
  x.coefs = as.matrix(X.inv.tr$coefs[,1:n.comp])
  rownames(x.coefs) = colnames(X)
  x.vars = X.inv.tr$vars[,1:n.comp]
  rownames(x.vars) = rownames(X)
  Y.inv.tr = PRCCA.inv.tr(Y.tr$mat, Y.tr$p.pen, sol$beta, Y.tr$tr)
  y.coefs = as.matrix(Y.inv.tr$coefs[,1:n.comp])
  rownames(y.coefs) = colnames(Y)
  y.vars = Y.inv.tr$vars[,1:n.comp]
  rownames(y.vars) = rownames(Y)
  #canonical correlation
  rho = diag(cor(X.inv.tr$vars,  Y.inv.tr$vars))[1:n.comp]
  names(rho) = paste('can.comp', 1:n.comp, sep = '')
  return(list('n.comp' = n.comp, 'cors' = rho, 'mod.cors' = rho.mod, 'x.coefs' = x.coefs, 'x.vars' = x.vars, 'y.coefs' = y.coefs, 'y.vars' = y.vars))
}



PRCCA.tr = function(X, index, lambda){
  p = ncol(X)
  n = nrow(X)
  if(lambda < 0){
    cat('please make', deparse(substitute(lambda)),'>= 0\n')
    return()
  }
  if(is.null(index)){
    lambda = 0
    index = 1:p
  } 
  if(lambda == 0 & p > n){
      cat('singularity issue. please impose a penalty on', deparse(substitute(X)), 'side\n')
      return()
  } 
  penalty = rep(0, p)
  permute = NULL
  B = NULL
  V = NULL
  p1 = length(index)
  if(lambda > 0){
    X1 = X[, index, drop = FALSE]
    X2 = X[, -index, drop = FALSE]
    permute = order(c(index, (1:p)[-index]))
    if(ncol(X2) >= n){
      cat('singularity issue. please penalize more',deparse(substitute(X)) ,'features\n')
      return()
    }
    if(ncol(X2) > 0){
      #regress X2 from X1
      B = solve(t(X2) %*% X2) %*% t(X2) %*% X1
      X1 = X1 - X2 %*% B
    }
    if(ncol(X1) > n){
      #apply kernel trick
      SVD = svd(X1)
      X1 = SVD$u %*% diag(SVD$d)
      V = SVD$v
      p1 = ncol(X1)
    }
    p1 = ncol(X1)
    X = cbind(X1, X2)
    penalty = rep(0, ncol(X))
    penalty[1:p1] = lambda
  }
  C = var(X, use = "pairwise")
  diag(C) = diag(C) + penalty
  return(list('mat' = X, 'p.pen' = p1, 'tr' = list('reg' = B, 'ker' = V, 'perm' = permute), 'cor' = C))
}



PRCCA.inv.tr = function(X, p1, alpha, tr){
  n.comp = ncol(alpha)
  colnames(alpha) = paste('can.comp', 1:n.comp, sep = '')
  V = tr$ker
  B = tr$reg 
  permute = tr$perm
  #find canonical variates
  u = X %*% alpha
  colnames(u) = paste('can.comp', 1:n.comp, sep = '')
  #inverse transfrom canonical coefficients
  if(!is.null(p1)){
    alpha1 = alpha[1:p1, , drop = FALSE]
    alpha2 = alpha[-(1:p1), , drop = FALSE]
    if(!is.null(V)) alpha1 = V %*% alpha1
    if(!is.null(B)) alpha2 = -B %*% alpha1 + alpha2
    alpha = rbind(alpha1, alpha2)
    if(!is.null(permute)) alpha = alpha[permute, ]
  }
  return(list('coefs' = alpha, 'vars' = u))
}
