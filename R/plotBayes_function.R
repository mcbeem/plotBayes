#' Function for exploring Bayesian inference for the sample mean for Normal data.
#'
#' This function plots the prior, likelihood, and posterior distribution of the sample mean. It reports the MAP (maximum a posteriori) and EAP (expected a posteriori) estimates of the mean as well as the 95\% credible interval around the MAP.
#' 
#' @param data A vector of data
#' @param prior.type The type of distribution for the prior. Can be "uniform" or "normal"
#' @param prior.parameters A vector of length two providing the parameter values for the prior. When prior.type="normal", the values are the mean and sd. When prior.type="uniform", the values are the min and the max.
#' @param min The minimum possible value of mu to consider.
#' @param max The maximum possible value of mu to consider.
#' @param credible The width of the credible interval specified as a number between zero and one. Defaults to 0.95 for a 95\% credible interval.
#' @param points The number of values of mu to be calculated, higher means more precision. Defaults to 1000.
#' @export
#' @return Returns a list with the following components:  
#' \itemize{
#'   \item \code{mean}    The sample mean of the data.  \cr 
#'   \item \code{map}    The maximum a posteriori estimate of the mean.  \cr 
#'   \item \code{eap}    The expected a posteriori estimate of the mean.  \cr
#'   \item \code{credible.interval} The values defining the endpoints of the Bayesian credible interval.
#'  }
#' @examples
#' set.seed(1)
#' data <- rnorm(n=10, mean=1, sd=1)
#' plotBayes(data, prior.type="normal", prior.parameters=c(0, .5), min=-2, max=2) 
#' plotBayes(data, prior.type="uniform", prior.parameters=c(.7, 1.5), min=-2, max=2)
#' plotBayes(data, prior.parameters=c(.7, 1.5), min=-2, max=2, credible=.9)

plotBayes <- function(data, prior.type="normal", prior.parameters, min, max, points=1000,
                      credible=.95) {
  
  # do some checks
  if (!(prior.type %in% c("normal", "uniform"))) {stop("prior.type must be one of \"normal\" or \"uniform\"")}
 
  if (length(prior.parameters) != 2) {stop("prior.parameters must be a vector of length two.")}

  
  if (prior.type=="uniform" & (min>prior.parameters[1])) {warning("min is greater than the minimum of the uniform prior")}
  if (prior.type=="uniform" & (max<prior.parameters[2])) {warning("max is less than the maximum of the uniform prior")}
  
  if (prior.type=="normal" & (min>prior.parameters[1])) {stop("min must be substantially less than the mean parameter of the prior (prior.parameters[1])")}
  if (prior.type=="normal" & (max<prior.parameters[1])) {stop("max must be substantially greater than the mean parameter of the prior (prior.parameters[1])")}
  if (prior.type=="normal" & (prior.parameters[2]<= 0)) {stop("The standard deviation parameter of the prior (prior.parameters[2]) must be greater than zero.")}
  
  if (credible <= 0 | credible >= 1) {stop("credible must be between zero and one")}
  
  # declare function for calculating the loglikelihood of the data given mu
  calc.LL <- function(data, mu, sd) {
    LL <- sum(log(dnorm(data, mu, sd)))
    return(LL)
  }
  
  # create a sequence of possible values of mu, ranging from MIN to MAX
  mus <- seq(min, max, length.out=points)
  
  # calculate distance between adjacent values of the possible mus
  width <- (max - min)/ points
  
  # calculate the log likelihood of the data given each value of mu
  LL <- sapply(mus, calc.LL, data=data, sd=sd(data))
  # convert LL to L
  L <- exp(LL)
           
  # prior distribution
  if (prior.type=="normal") {
    prior <- dnorm(mus, prior.parameters[1], prior.parameters[2])
  }
  
  if (prior.type=="uniform") {
    prior <- dunif(mus, min=prior.parameters[1], max=prior.parameters[2])
  }
  
  # normalized posterior (to area of 1, probability distribution)
  posterior <- L*prior
  posterior.area <- sum(posterior*width)
  normalized.posterior <- posterior / posterior.area
  
  # normalized likelihood for plot
  normalized.L <- L / sum(L*width)
  
  # find global max to scale yaxis
  ymax <- max(c(prior, normalized.posterior, normalized.L))

  plot(mus, prior, type = "l", 
       xlim = c(min-((max-min)/10), max+((max-min)/10)), ylim = c(0, ymax), lty = "dotted", 
       col = "gray48", xlab = bquote(mu), ylab = "Density", las = 1)
  lines(mus, normalized.posterior, type = "l", col = "skyblue", lwd = 3)
  lines(mus, normalized.L, type = "l", col = "red", lty = 2)
  lines(rev(mus), rev(prior), type="l", lty="dotted", col="gray48")
  legend("topleft", c("Prior", "Likelihood", "Posterior"), lty = c(2, 2, 1), 
         col = c("gray48", "red", "skyblue"), lwd = c(1, 1, 4), bty = "n")

  # find credible interval locations
  
  # define cumulative density function
  calc.cum.dens <- function(i) {
    # calculate the cumulative density below the index value
    return(sum(normalized.posterior[1:i])*width)
  }
  
  # indexes is a vector of all the indexes of the posterior
  indexes <- 1:length(normalized.posterior)
  
  # calculate the area below each value
    cumulative.density <- sapply(indexes, calc.cum.dens)
    
  # find the limits
  # percentile bounds given the supplied value of argument 'credible'
  percentiles <- c((1-credible)/2, 1-(1-credible)/2)
    # index of the credible interval lower limit
  ll.index <- which(abs(cumulative.density-percentiles[1])==min(abs(cumulative.density-percentiles[1])))
   # index of the credible interval upper limit
  ul.index <- which(abs(cumulative.density-percentiles[2])==min(abs(cumulative.density-percentiles[2])))
  
  return(list(
    data.mean=mean(data),
    map=mus[which(normalized.posterior == max(normalized.posterior))],
    eap=sum((normalized.posterior/sum(normalized.posterior)*mus)),
    credible.intercal=c(mus[ll.index], mus[ul.index])
  ))
}

