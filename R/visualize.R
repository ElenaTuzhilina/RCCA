#' @title Reconstruction vizualization
#'
#' @description This function allows to plot reconstructed chromatin conformation \emph{X} and corresponding contact matrix approximation \emph{XX'}.
#'
#' @param X a matrix representing spatial coordinates of resulting chromatin reconstruction.
#' @param index points where spline basis is evaluated; each corresponds to a particular genomic loci.
#' @param type the type of plot returned. Set \code{type = 'projection'} and \code{type = '3D'} to output the projection and 3D model of chromatin conformation reconstruction, respectively.
#' @param title optional, adds title to the plot. Default value \code{title = NULL}.
#'
#' @return Reconstruction plots.
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
#' #run WPCMS with equal weights 
#' wpcms = WPCMS(Z, H)
#' 
#' #plot projection of reconstructed chromatin conformation
#' visualize(wpcms$X, type = 'projection')
#' 
#' #plot 3D model of reconstructed chromatin conformation
#' visualize(wpcms$X, type = '3D')
#'
#' @export visualize

visualize = function(X, type = 'projection', index = 1:nrow(X), title = NULL){
  n = nrow(X)
  before_centromere = which(index < (n * 0.45))
  after_centromere = which(index >= (n * 0.45))
  col = c(rep('orange', length(before_centromere)), rep('darkturquoise', length(after_centromere)))
  colnames(X) = c('x', 'y', 'z')
  par(mfrow = c(1,1), oma = c(0, 0, 2, 0))
  panelf = function(x, y){
    col = c(rep('orange', length(before_centromere)), rep('darkturquoise', length(after_centromere)))
    points(x, y, pch = 19, cex = 1, col = col)
    lines(x, y, col = 'orange', lwd = 2)
    lines(x[after_centromere], y[after_centromere], col = 'darkturquoise', lwd = 2)
  }
  if(type == 'projection') return(pairs(X, panel = panelf, cex.labels = 5, main = title))

  if(type == '3D') return(plotly::layout(plotly::plot_ly(x = X[,1], y = X[,2], z = X[,3], type = 'scatter3d', mode = 'lines+markers',
              line = list(width = 6, color = col), marker = list(size = 3.5, color = col)), title = title))
  
}

