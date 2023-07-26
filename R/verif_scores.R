#' CRPS of WMO-reported visibility forecasts based on the estimated probabilities
#' and verifying observations
#'
#' @description Calculates the CRPS
#' @param prob Estimated probabilities. Dimensions: \code{[84, 21]}, where \code{84} is the number of observation categories and 21 is the number of lead times.
#' @param obs A \code{[1:21]}-sized array, containing the verifying observations.
#' @return CRPS for a specific location.

crps.vis.wmo <- function(prob, obs) {
  Z <-
    c(seq(0, 5000, by = 100),
      seq(6000, 30000, by = 1000),
      seq(35000, 70000, by = 5000))
  crps1 <-
    apply(t(prob) * abs(obs - matrix(
      Z,
      nrow = length(obs),
      ncol = 84,
      byrow = T
    )), 1, sum)
  Zmat <- matrix(Z,
                 nrow = 84,
                 ncol = 84,
                 byrow = T)
  crps2 <- apply((abs(Zmat - t(Zmat)) %*% prob) * prob, 2, sum)
  crps <- crps1 - crps2 / 2

  return(crps)
}


#' Logarithmic Score of WMO-reported visibility forecasts based on the estimated probabilities
#' and verifying observations
#'
#' @param prob Estimated probabilities. Dimensions: \code{[84, 21]}, where \code{84} is the number of observation categories and 21 is the number of lead times.
#' @param obs A \code{[1:21]}-sized array, containing the verifying observations.
#'
#' @return Logarithmic score for a specific location

logs.vis.wmo <- function(prob, obs) {
  good.ind <- !is.na(obs)
  N <- length(obs)
  logs <- rep(NA, N)
  obsind <- as.character(obs[good.ind])
  prob <- prob[, good.ind]

  if (sum(good.ind) == 1) {
    logs[good.ind] <- -log(prob[obsind])
  } else if (sum(good.ind) > 1) {
    logs[good.ind] <- -log(diag(prob[obsind,]))
  }

  return(logs)
}

#' Estimated variance of the discrete probability distribution of WMO-reported visibility forecasts
#'
#' @param prob Estimated probabilities. Dimensions: \code{[84, 21]}, where \code{84} is the number of observation categories and 21 is the number of lead times.
#'
#' @return Estimated variances for the 21 different lead times

var.vis.wmo <- function(prob) {
  Z <-
    c(seq(0, 5000, by = 100),
      seq(6000, 30000, by = 1000),
      seq(35000, 70000, by = 5000))

  MEAN <- Z %*% prob
  sqMEAN <- Z ^ 2 %*% prob
  VAR <- sqMEAN - MEAN ^ 2

  return(VAR)
}

#' Quantiles of the discrete probability distribution of WMO-reported visibility forecasts
#'
#' @param prob Estimated probabilities. Dimensions: \code{[84, 21]}, where \code{84} is the number of observation categories and 21 is the number of lead times.
#' @param q level of quantile
#'
  #' @return Estimated quantiles of the discrete distribution for the 21 different lead times

quant.vis.wmo <- function(prob, q) {
  Z <-
    c(seq(0, 5000, by = 100),
      seq(6000, 30000, by = 1000),
      seq(35000, 70000, by = 5000))

  sum.prob <- apply(prob, 2, cumsum)
  ind.q <- apply(sum.prob < q, 2, sum) + 1

  return(Z[ind.q])
}

#' PIT values based on the estimated probabilities
#' and verifying observations
#'
#' @param probs Estimated probabilities. Dimensions: \code{[84, 21]}, where \code{84} is the number of observation categories and 21 is the number of lead times.
#' @param obs A \code{[1:21]}-sized array, containing the verifying observations.
#'
#' @return PIT values for the 21 different lead times

pit.vis.wmo <- function(prob, obs) {
  N <- length(obs)

  good.ind <- !is.na(obs)

  obsind <- as.character(obs[good.ind])
  prob <- prob[, good.ind]

  Ngood <- sum(good.ind)

  if (Ngood == 1) {
    M <- length(prob)
    auxInd <- c(1:M)
    names(auxInd) <- names(prob)
    auxInd <- auxInd[obsind]
    prob <- c(0, prob)
    sum.prob <- cumsum(prob)
    boundL <- sum.prob[auxInd]
    boundU <- sum.prob[auxInd + 1]
  } else if (Ngood > 1) {
    M <- dim(prob)[1]
    auxInd <- c(1:M)
    names(auxInd) <- dimnames(prob)[[1]]
    auxInd <- auxInd[obsind]
    prob <- rbind(rep(0, Ngood), prob)
    sum.prob <- apply(prob, 2, cumsum)
    boundL <- diag(sum.prob[auxInd,])
    boundU <- diag(sum.prob[auxInd + 1,])
  }

  PIT <- rep(NA, N)
  if (Ngood > 0) {
    PIT[good.ind] <- runif(rep(1, Ngood), min = boundL, max = boundU)
  }

  return(PIT)
}
