#' @title Canonical Correlation analysis with group L2 penalty
#' 
#' @description GRCCA function performs Canonical Correlation Analysis with group L2 regularization and allows to conduct Canonical Correlation Analysis in high dimensions. 
#' It imposes group L2 penalty on the coefficient vectors \eqn{\alpha} and \eqn{\beta} coefficients. Specifically, if
#' \deqn{x = (x_1, \ldots, x_p)~~and~~y = (y_1, \ldots, y_q)}{x = (x_1, ..., x_p)  and  y = (y_1, ..., y_q)} 
#' are random vectors and 
#' \deqn{I_1, ..., I_K  \subset \{1, ..., p\}~~and~~J_I, ... , J_L \subset \{1, ..., q\}}{I_1, ..., I_K is a partition of {1, ..., p}  and  J_1, ..., J_L  is a partition of  {1, ..., q}} 
#' then
#' GRCCA seeks for such vectors 
#' \deqn{\alpha = (\alpha_1, \ldots, \alpha_p)~~and~~\beta = (\beta_1, \ldots, \beta_q)}{\alpha = (\alpha_1, ..., \alpha_p)  and  \beta = (\beta_1, ..., \beta_q)} 
#' that satisfy two within group constraints 
#' \deqn{\|\alpha\|^2_{within} = \sum_k\|\alpha_{I_k} - \bar{\alpha_{I_k}}\|^2\leq t_1}{||\alpha||_w = var(\alpha_I_1) + ... + var(\alpha_I_K) <= t_1}
#' and
#' \deqn{\|\beta\|^2_{within} = \sum_l\|\beta_{J_l} - \bar{\beta_{J_l}}\|^2\leq t_2}{||\beta||_w = var(\beta_J_1) + ... + var(\beta_J_L) <= t_2}
#' as well as two between group constraints
#' \deqn{\|\alpha\|^2_{between} = \sum_k |I_k| \bar{\alpha_{I_k}}^2\leq s_1}{||\alpha||_b = |I_1| mean(\alpha_I_1)^2 + ... + |I_K| mean(\alpha_I_K)^2 <= s_1}
#' and
#' \deqn{\|\beta\|^2_{between} = \sum_l |J_l| \bar{\beta_{J_l}}^2\leq s_2}{||\beta||_b = |J_1| mean(\beta_J_1)^2 + ... + |J_L| mean(\beta_J_L)^2 <= s_2}
#' and that maximize the correlation 
#' \eqn{cor(u, v)} 
#' between the linear combnations 
#' \deqn{u = \langle x,\alpha\rangle and v = \langle y,\beta\rangle}{u = <x , \alpha> and v = <y , \beta>.} 
#' Here 
#' \eqn{<a , b>} 
#' refers to the inner product between two vectors;
#' \deqn{\alpha_{I_k}~~and~~\beta_{J_l}}{\alpha_I_k  and  \beta_J_l} 
#' are corresponding subvectors of \eqn{\alpha} and \eqn{\beta} with indices belonging to \eqn{I_k} and \eqn{J_l}, respectively;
#' and \eqn{|A|} referes to the set candinality.
#' The above optimization problem is equivalet to maximizing the modified correlation coefficient 
#' \deqn{\frac{cov(\langle x, \alpha\rangle, \langle y, \beta\rangle)}{\sqrt{var(\lamgle x, \alpha\rangle) + \lambda_1\|\alpha\|^2_{within} + \mu_1\|\alpha\|^2_{between}}\sqrt{var(\langle y, \beta\rangle) + \lambda_2\|\beta\|^2_{within} + \mu_2\|\beta\|^2_{between}}}}{cov(<x , \alpha>, <y , \beta>)  /  ( cov(<x , \alpha>) + \lambda_1  ||\alpha||_w + \mu_1 ||\alpha||_b )^1/2 ( var(<y , \beta>) + \lambda_2 ||\beta||_w + \mu_2 ||\beta||_b )^1/2,}
#' where \deqn{\lambda_1~~and~~\lambda_2}{\lambda_1  and  \lambda_2} 
#' control the resulting within group variation of the coefficiens and
#' \deqn{\mu_1~~and~~\mu_2}{\mu_1  and  \mu_2}
#' control the sparsity on a group level of the canonical coefficients \eqn{\alpha} and \eqn{\beta}. 
#' 
#' 
#' @param X a rectangular \eqn{n x p} matrix containing \eqn{n} observations of random vector \eqn{x}. 
#' @param Y a rectangular \eqn{n x q} matrix containing \eqn{n} observations of random vector \eqn{y}.
#' @param group1 an integer valued vector representing the group assignment of \eqn{\alpha} coefficients. By default \code{group1 = rep(1, ncol(X))}, i.e. we include all \eqn{\alpha} coefficients in the same group.
#' @param group2 an integer valued vector representing the group assignment of \eqn{\beta} coefficients. By default \code{group2 = rep(1, ncol(Y))}, i.e. we include all \eqn{\beta} coefficients in the same group.
#' @param lambda1 a non-negative penalty factor used for controlling the within variation of \eqn{\alpha} coefficients. By default \code{lambda1 = 0}, i.e. no regularization is imposed. Increasing \code{lambda1} shrinks each coefficient toward it's group mean.
#' @param lambda2 a non-negative penalty factor used for controlling the within variation of \eqn{\beta} coefficients. By default \code{lambda2 = 0}, i.e. no regularization is imposed. Increasing \code{lambda2} shrinks each coefficient toward it's group mean.
#' @param mu1 a non-negative penalty factor used for controlling the between variation of \eqn{\alpha} coefficients. By default \code{mu1 = 0}, i.e. no regularization is imposed. Increasing \code{mu1} shrinks each group mean toward zero.
#' @param mu2 a non-negative penalty factor used for controlling the between variation of \eqn{\beta} coefficients. By default \code{mu2 = 0}, i.e. no regularization is imposed. Increasing \code{mu2} shrinks each group mean toward zero.
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
#' n.groups = 5
#' #run GRCCA with no sparsity on a group level
#' group1 = rep(1:n.groups, rep(ncol(X)/n.groups, n.groups))
#' grcca = GRCCA(X, Y, group1, group2 = NULL, lambda1 = 1000, lambda2 = 0, mu1 = 0, mu2 = 0)
#' #check the modified canonical correlations 
#' plot(1:grcca$n.comp, grcca$mod.cors, pch = 16, xlab = 'component', 'ylab' = 'correlation', ylim = c(0, 1))
#' #check the canonical coefficients for the first canonical variates
#' barplot(grcca$x.coefs[,'can.comp1'], col = 'orange', 'xlab' = 'X feature', ylab = 'value')
#' n.groups = 50
#' #run GRCCA with sparsity on a group level
#' group1 = rep(1:n.groups, rep(ncol(X)/n.groups, n.groups))
#' grcca = GRCCA(X, Y, group1, group2 = NULL, lambda1 = 10000, lambda2 = 0, mu1 = 100, mu2 = 0)
#' #check the modified canonical correlations 
#' plot(1:grcca$n.comp, grcca$mod.cors, pch = 16, xlab = 'component', 'ylab' = 'correlation', ylim = c(0, 1))
#' #check the canonical coefficients for the first canonical variates
#' barplot(grcca$x.coefs[,'can.comp1'], col = 'orange', 'xlab' = 'X feature', ylab = 'value')
#' 
#' @export GRCCA
#' 

GRCCA = function(X, Y, group1 = rep(1, ncol(X)), group2 = rep(1, ncol(Y)), lambda1 = 0, lambda2 = 0, mu1 = 0, mu2 = 0){
  X = as.matrix(X)
  Y = as.matrix(Y)
  n.comp = min(ncol(X), ncol(Y), nrow(X))
  
  #transform
  X.tr = GRCCA.tr(X, group1, lambda1, mu1)
  Y.tr = GRCCA.tr(Y, group2, lambda2, mu2)
  if(is.null(X.tr) | is.null(Y.tr)) return(invisible(NULL))
  
  #solve optimization problem
  sol = PRCCA(X.tr$mat, Y.tr$mat, X.tr$ind, Y.tr$ind, 1, 1)
  mod.cors = sol$mod.cors[1:n.comp]
  
  #inverse transform
  X.inv.tr = GRCCA.inv.tr(sol$x.coefs, group1, lambda1, mu1)
  x.coefs = as.matrix(X.inv.tr$coefs[,1:n.comp])
  rownames(x.coefs) = colnames(X)
  x.vars = sol$x.vars[,1:n.comp]
  rownames(x.vars) = rownames(X)
  Y.inv.tr = GRCCA.inv.tr(sol$y.coefs, group2, lambda2, mu2)
  y.coefs = as.matrix(Y.inv.tr$coefs[,1:n.comp])
  rownames(y.coefs) = colnames(Y)
  y.vars = sol$y.vars[,1:n.comp]
  rownames(y.vars) = rownames(Y)
  rho = diag(cor(x.vars,  y.vars))[1:n.comp]
  names(rho) = paste('can.comp', 1:n.comp, sep = '')
  return(list('n.comp' = n.comp, 'cors' = rho, 'mod.cors' = mod.cors, 'x.coefs' = x.coefs, 'x.vars' = x.vars, 'y.coefs' = y.coefs, 'y.vars' = y.vars))
}
  

GRCCA.tr = function(X, group, lambda, mu){
  if(is.null(group)){
    lambda = 0
    mu = 0
  } 
  if(lambda < 0){
    cat('please make', deparse(substitute(lambda)),'> 0\n')
    return(NULL)
  }
  if(mu < 0){
    cat('please make', deparse(substitute(mu)),'> 0\n')
    return(NULL)
  }
  index = NULL
  if(lambda > 0){
    #extend matrix
    group.names = unique(sort(group))
    ps = table(group)
    agg = aggregate(t(X), by = list(group), FUN = mean)
    X.mean = t(agg[, -1])
    colnames(X.mean) = agg[, 1] 
    if(mu == 0){
      mu = 1
      index = 1:ncol(X)
    } else {
      index = 1:(ncol(X) + ncol(X.mean))
    }
    X1 = 1/sqrt(lambda) * (X - X.mean[,group])
    X2 = scale(X.mean[,group.names], center = FALSE, scale = sqrt(mu/ps[group.names]))
    X = cbind(X1, X2)
  }
  return(list('mat' = X, 'ind' = index))
}

GRCCA.inv.tr = function(alpha, group, lambda, mu){
  if(is.null(group)){
    lambda = 0
    mu = 0
  } 
  if(lambda > 0){
    p = length(group)
    group.names = unique(sort(group))
    ps = table(group)
    alpha1 = alpha[1:p,]
    alpha2 = alpha[-(1:p),]
    agg = aggregate(alpha1, by = list(group), FUN = mean)
    alpha1.mean = agg[, -1]
    rownames(alpha1.mean) = agg[, 1]
    alpha1 = 1/sqrt(lambda) * (alpha1 - alpha1.mean[group,])
    if(mu == 0) mu = 1
    alpha2 = t(scale(t(alpha2[group.names,]), center = FALSE, scale = sqrt(mu*ps[group.names])))
    alpha = alpha1 + alpha2[group,]
  }
  return(list('coefs' = alpha))
}
