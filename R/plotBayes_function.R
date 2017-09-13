#' Function for exploring Bayesian inference for the sample mean for continuous data.
#'
#' This function plots the prior, likelihood, and posterior distribution of the sample mean. It reports the MAP (maximum a posteriori) estimate of the mean as well as the 95\% credible interval.
#' 
#' @param data A vector of data
#' @param prior.type The type of distribution for the prior. Can be "uniform" or "normal"
#' @param prior.parameters A vector of length two providing the parameter values for the prior. When prior.type="normal", the values are the mean and sd. When prior.type="uniform", the values are the min and the max.
#' @param min The minimum possible value of mu to consider.
#' @param max The maximum possible value of mu to consider.
#' @param points The number of values of mu to be calculated, higher means more precision. Defaults to 1000.
#' @export
#' @examples
#' set.seed(1)
#' data <- rnorm(n=10, mean=1, sd=1)
#' plotBayes(data, prior.type="normal", prior.parameters=c(0, .5), min=-2, max=2) 
#' plotBayes(data, prior.type="uniform", prior.parameters=c(.7, 1.5), min=-2, max=2)

plotBayes <- function(data, prior.type="normal", prior.parameters, min, max, points=1000) {
  
  # do some checks
  if (!(prior.type %in% c("normal", "uniform"))) {stop("prior.type must be one of \"normal\" or \"uniform\"")}
 
  if (length(prior.parameters) != 2) {stop("prior.parameters must be a vector of length two.")}

  
  if (prior.type=="uniform" & (min>prior.parameters[1])) {stop("min must be less than the minimum of the uniform prior, prior.parameters[1]")}
  if (prior.type=="uniform" & (max<prior.parameters[2])) {stop("max must be greater than the maximim of the uniform prior, prior.parameters[2]")}
  
  if (prior.type=="normal" & (min>prior.parameters[1])) {stop("min must be substantially less than the mean parameter of the prior (prior.parameters[1])")}
  if (prior.type=="normal" & (max<prior.parameters[1])) {stop("max must be substantially greater than the mean parameter of the prior (prior.parameters[1])")}
  if (prior.type=="normal" & (prior.parameters[2]<= 0)) {stop("The standard deviation parameter of the prior (prior.parameters[2]) must be greater than zero.")}
  
  
  # declare function for calculating the likelihood of the data given mu
  calc.L <- function(data, mu, sd) {
    LL <- sum(log(dnorm(data, mu, sd)))
    return(LL)
  }
  
  # create a sequence of possible values of mu, ranging from MIN to MAX
  mus <- seq(min, max, length.out=points)
  
  # calculate distance between adjacent values of the possible mus
  width <- (max - min)/ points
  
  # calculate the log likelihood of the data given each value of mu
  LL <- sapply(mus, calc.L, data=data, sd=sd(data))
  L <- exp(LL - mean(LL))
  
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

  
  return(list(
    data.mean=mean(data),
    map=mus[which(normalized.posterior == max(normalized.posterior))],
    eap=sum((normalized.posterior/sum(normalized.posterior)*mus)),
    cred.95.LL= mus[which(abs(cumsum(normalized.posterior*width)-.05) == 
                            min(abs(cumsum(normalized.posterior*width)-.05)))],
    cred.95.UL= mus[which(abs(cumsum(normalized.posterior*width)-.95) == 
                            min(abs(cumsum(normalized.posterior*width)-.95)))]))
}