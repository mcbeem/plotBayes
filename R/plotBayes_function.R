#' Function for exploring Bayesian inference for the sample mean for normal data.
#'
#' This function plots the prior, likelihood, and posterior distribution of the sample mean. It reports the MAP (maximum a posteriori) estimate of the mean as well as the 95\% credible interval.
#' 
#' @param data A vector of data.
#' @param prior.mean The mean of the prior distribution for the mean parameter mu.
#' @param prior.sd The SD of the prior distribution for the mean. Smaller means more certainty.
#' @param min The minimum possible value of mu to consider.
#' @param max The maximum possible value of mu to consider.
#' @param points The number of values of mu to be calculated, higher means more precision. Defaults to 1000.
#' @export
#' @examples
#' set.seed(1)
#' data <- rnorm(n=10, mean=1, sd=1)
#' plotBayes(data, prior.mean=0, prior.sd=.2, min=-2, max=2)

plotBayes <- function(data, prior.mean, prior.sd, min, max, points=1000) {
  
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
  prior <- dnorm(mus, prior.mean, prior.sd)
  
  # normalized posterior (to area of 1, probability distribution)
  posterior <- L*prior
  posterior.area <- sum(posterior*width)
  normalized.posterior <- posterior / posterior.area
  
  # normalized likelihood for plot
  normalized.L <- L / sum(L*width)
  
  # find global max to scale yaxis
  ymax <- max(c(prior, normalized.posterior, normalized.L))
  
  plot(mus, prior, type = "l", 
       xlim = c(min-((max-min)/10), max), ylim = c(0, ymax), lty = 2, 
       col = "gray48", xlab = bquote(theta), ylab = "Density", las = 1)
  lines(mus, normalized.posterior, type = "l", col = "skyblue", lwd = 3)
  lines(mus, normalized.L, type = "l", col = "red", lty = 2)
  legend("topleft", c("Prior", "Likelihood", "Posterior"), lty = c(2, 2, 1), 
         col = c("gray48", "red", "skyblue"), lwd = c(1, 1, 4), bty = "n")
  
  return(list(
    data.mean=mean(data),
    map=mus[which(normalized.posterior == max(normalized.posterior))],
    cred.95.LL= mus[which(abs(cumsum(normalized.posterior*width)-.05) == 
                            min(abs(cumsum(normalized.posterior*width)-.05)))],
    cred.95.UL= mus[which(abs(cumsum(normalized.posterior*width)-.95) == 
                            min(abs(cumsum(normalized.posterior*width)-.95)))]))
}