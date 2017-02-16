
#' Plot \code{miSum}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot miSum
#'
#' @param x Object \code{miSum}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param covs Data frame consisting of a single row indicating the covariate values to be used in plots. If none are specified, the means of any covariates appearing in the model are used (unless covariate is a factor, in which case the first factor is used).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim Parameter \code{ylim} for the step length histograms.
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values.
#' @param sepAnimals If \code{TRUE}, the data is split by individuals in the histograms.
#' Default: \code{FALSE}.
#' @param sepStates If \code{TRUE}, the data is split by states in the histograms.
#' Default: \code{FALSE}.
#' @param col Vector or colors for the states (one color per state).
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @details The state-dependent densities are weighted by the frequency of each state in the most
#' probable state sequence (decoded with the function \code{\link{viterbi}}). For example, if the
#' most probable state sequence indicates that one third of observations correspond to the first
#' state, and two thirds to the second state, the plots of the densities in the first state are
#' weighted by a factor 1/3, and in the second state by a factor 2/3.
#'
#' @examples
#' # m is a miSum object (as returned by MIfitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)
#'
#'
#' @export

plot.miSum <- function(x,animals=NULL,covs=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                        sepStates=FALSE,col=NULL,plotCI=FALSE,alpha=0.95,plotEllipse=TRUE,...)
{
  m <- x # the name "x" is for compatibility with the generic method
  m$mle <- lapply(x$Par$real,function(x) x$est)
  m$mle$beta <- x$Par$beta$beta$est
  m$mle$delta <- x$Par$real$delta$est
  m$mod <- list()
  m$mod$estimate <- x$MIcombine$coefficients
  m$mod$hessian <- ginv(x$MIcombine$variance)
  m$plotEllipse <- plotEllipse
  
  class(m) <- append("momentuHMM",class(m))
  
  plot.momentuHMM(m,animals,covs,ask,breaks,hist.ylim,sepAnimals,
                  sepStates,col,plotCI,alpha,...)
}
