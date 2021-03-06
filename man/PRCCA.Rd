% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRCCA.R
\name{PRCCA}
\alias{PRCCA}
\title{Canonical Correlation analysis with partial L2 penalty}
\usage{
PRCCA(X, Y, index1 = 1:ncol(X), index2 = 1:ncol(Y), lambda1 = 0, lambda2 = 0)
}
\arguments{
\item{X}{a rectangular \eqn{n x p} matrix containing \eqn{n} observations of random vector \eqn{x}.}

\item{Y}{a rectangular \eqn{n x q} matrix containing \eqn{n} observations of random vector \eqn{y}.}

\item{index1}{a subset of indices the penalty is imposed on while regularizing the \eqn{X} side. By default \code{index1 = 1:ncol(X)}, i.e. we include all \eqn{\alpha} coefficients in the relularization term.}

\item{index2}{a subset of indices the penalty is imposed on while regularizing the \eqn{Y} side. By default \code{index2 = 1:ncol(Y)}, i.e. we include all \eqn{\beta} coefficients in the relularization term.}

\item{lambda1}{a non-negative penalty factor used for regularizing \eqn{X} side coefficients \eqn{\alpha}. By default \code{lambda1 = 0}, i.e. no regularization is imposed. Increasing \code{lambda1} incourages sparsity of the resulting canonical coefficients.}

\item{lambda2}{a non-negative penalty factor used for regularizing \eqn{Y} side coefficients \eqn{\beta}. By default \code{lambda2 = 0}, i.e. no regularization is imposed. Increasing \code{lambda2} incourages sparsity of the resulting canonical coefficients.}
}
\value{
A list containing the PCMS problem solution:
\itemize{
  \item \code{n.comp} -- the number of computed canonical components, i.e. \eqn{k = min(p, q)}.
  \item \code{cors} -- the resulting \eqn{k} canonical correlations.
  \item \code{mod.cors} -- the resulting \eqn{k} values of modified canonical correlation.
  \item \code{x.coefs} -- \eqn{p x k} matrix representing \eqn{k} canonical coefficient vectors \eqn{\alpha[1]}, ..., \eqn{\alpha[k]}.
  \item \code{x.vars} -- \eqn{n x k} matrix representing \eqn{k} canonical variates \eqn{u[1]}, ..., \eqn{u[k]}.
  \item \code{y.coefs} -- \eqn{q x k} matrix representing \eqn{k} canonical coefficient vectors \eqn{\beta[1]}, ..., \eqn{\beta[k]}.
  \item \code{y.vars} -- \eqn{n x k} matrix representing \eqn{k} canonical variates \eqn{v[1]}, ..., \eqn{v[k]}.
}
}
\description{
PRCCA function performs Canonical Correlation Analysis with partial L2 regularization and allows to conduct Canonical Correlation Analysis in high dimensions. 
It imposes L2 penalty only on a subset of \eqn{\alpha} and \eqn{\beta} coefficients. Specifically, if
\deqn{x = (x_1, \ldots, x_p)~~and~~y = (y_1, \ldots, y_q)}{x = (x_1, ..., x_p)  and  y = (y_1, ..., y_q)} 
are random vectors and 
\deqn{I = \{i_1, ..., i_m\} \subset \{1, ..., p\}~~and~~J = {j_1, ..., i_r} \subset \{1, ..., q\}}{I = {i_1, ..., i_m}  is a subset of  {1, ..., p}  and  J = {j_1, ..., i_r}  is a subset of  {1, ..., q}} 
then
PRCCA seeks for such vectors 
\deqn{\alpha = (\alpha_1, \ldots, \alpha_p)~~and~~\beta = (\beta_1, \ldots, \beta_q)}{\alpha = (\alpha_1, ..., \alpha_p)  and  \beta = (\beta_1, ..., \beta_q)} 
that satisfy partial L2 constraints 
\deqn{\|\alpha_I\|\leq t_1~~and~~\|\beta_J\| \leq t_2}{||\alpha_I|| <= t_1  and  ||\beta_J|| <= t_2}
and that maximize the correlation 
\eqn{cor(u, v)} 
between the linear combnations 
\deqn{u = \langle x,\alpha\rangle and v = \langle y,\beta\rangle}{u = <x , \alpha> and v = <y , \beta>.} 
Here 
\eqn{<a , b>} 
refers to the inner product between two vectors and 
\deqn{\alpha_I~~and~~\beta_J}{\alpha_I  and  \beta_J} 
are corresponding subvectors of \eqn{\alpha} and \eqn{\beta} with indices belonging to \eqn{I} and \eqn{J}, respectively.
Again, the above optimization problem is equivalet to maximizing the modified correlation coefficient 
\deqn{\frac{cov(\langle x, \alpha\rangle, \langle y, \beta\rangle)}{\sqrt{var(\lamgle x, \alpha\rangle) + \lambda_1\|\alpha_I\|^2}\sqrt{var(\langle y, \beta\rangle) + \lambda_2\|\beta_J\|^2}}}{cov(<x , \alpha>, <y , \beta>)  /  ( cov(<x , \alpha>) + \lambda_1  ||\alpha_I||^2 )^1/2 ( var(<y , \beta>) + \lambda_2 ||\beta_J||^2 )^1/2,}
where \deqn{\lambda_1~~and~~\lambda_2}{\lambda_1  and  \lambda_2} control the resulting sparsity of the canonical coefficients within \deqn{\alpha_I~~and~~\beta_J}{\alpha_I  and  \beta_J} parts of the coefficient vectors.
}
\examples{
data(X)
data(Y)
#run RCCA 
prcca = PRCCA(X, Y, lambda1 = 100, lambda2 = 0, index1 = 1:(ncol(X) - 10))
#check the modified canonical correlations 
plot(1:prcca$n.comp, prcca$mod.cors, pch = 16, xlab = 'component', 'ylab' = 'correlation', ylim = c(0, 1))
#check the canonical correlations
points(1:prcca$n.comp, prcca$cors, pch = 16, col = 'purple')
#compare them to cor(x*alpha, y*beta)
points(1:prcca$n.comp, diag(cor(X \%*\% prcca$x.coefs, Y \%*\% prcca$y.coefs)), col = 'cyan', pch = 16, cex = 0.7)
#check the canonical coefficients for the first canonical variates
barplot(prcca$x.coefs[,'can.comp1'], col = 'orange', 'xlab' = 'X feature', ylab = 'value')
barplot(prcca$y.coefs[,'can.comp1'], col = 'darkgreen', 'xlab' = 'Y feature', ylab = 'value')

}
