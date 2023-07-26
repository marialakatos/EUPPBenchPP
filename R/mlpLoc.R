require(doParallel)

#' Fitting MLP-L model to the EUPPBench visibility ensemble
#'
#' @name mlpProbsVisLoc
#' @description Fitting MLP-L model to the EUPPBench datasets' visibility forecasts.
#' Architecture: 2 hidden layers (25x25 neurons), Activation function: Logistic
#' Input features:
#' \itemize{
#'  \item control forecast (normalized)
#'  \item ensemble mean (normalized)
#'  \item ensemble variance (normalized)
#'  \item proportions of forecasts predicting visibility between 0 and 1000 meters
#'  \item proportions of forecasts predicting visibility between 1001 and 2000 meters
#'  \item proportions of forecasts predicting visibility above 30000 meters
#' }
#'
#' @param visData Visibility data. Dimensions:
#' @param varNames Feature names
#' @param nStats Number of stations
#' @param nObsNames Number of observation names
#' @param nLeadTimes Number of lead times
#' @param stationIDs A character vector of station ids.
#' @param obsNames A character vector of observation names
#' @param leadTimes A character vector of lead times
#' @param start.day Index of the first day of the verification period
#' @param end.day Index of the last day of the verification period
#' @param train Length of training period (default value is 30 days)
#'
#' @return Estimated probabilities
#' @export
mlpProbsVisLoc <-
  function(visData,
           varNames,
           nStats,
           nObsNames,
           nLeadTimes,
           stationIDs,
           obsNames,
           leadTimes,
           start.day,
           end.day,
           train = 350) {


    start.day <- 366
    end.day <- 730

    forcDates <- obsDates[start.day:end.day]
    nForcDates <- length(forcDates)

    warning("In case of force stopping, the opened cluster of connections will not be closed.")
    cores <- detectCores()
    c1 <- makeCluster(cores[1] - 4)
    registerDoParallel(c1)

    probs <- array(
      data = NA,
      dim = c(nForcDates, nStats, nObsNames, nLeadTimes),
      dimnames = list(forcDates, stationIDs, obsNames, leadTimes)
    )

    for (lead in leadTimes[1:nLeadTimes]) {
      if (as.numeric(lead) == 0) {
        lt <- 1
      } else {
        lt <- ceiling(as.numeric(lead) / 24)
      }

      cat(
        "Modelling for lead time",
        lead, "hours is in progress ...\n"
      )

      predTmp <-
        foreach(
          d = (start.day):(end.day),
          .combine = "acomb",
          .multicombine = TRUE,
          .packages = "RSNNS"
        ) %dopar% {
          train.ind <- (d - train - lt + 1):(d - lt)
          data.train <- visData[train.ind, , , lead]
          per1 <- sin(2 * pi * train.ind / 365)
          per2 <- cos(2 * pi * train.ind / 365)
          data.test <- visData[d, , , lead]

          probTmp <-
            array(
              data = NA,
              dim = c(nStats, nObsNames),
              dimnames = list(stationIDs, obsNames)
            )

          for (act_gr in c(1:nStats)) {
            train.aux <- cbind(data.train[, act_gr, ], per1, per2)
            colnames(train.aux) <- varNames
            mlpTrain <- data.frame(train.aux)

            mlpTrain <- mlpTrain[complete.cases(mlpTrain), ]
            mlpTrain$obs <- decodeClassLabels(mlpTrain$obs)
            test.aux <-
              c(data.test[act_gr, ], sin(2 * pi * d / 365), cos(2 * pi * d / 365))
            names(test.aux) <- varNames
            mlpTest <- data.frame(t(test.aux))
            mlpTest$obs <- decodeClassLabels(mlpTest$obs)

            if (!anyNA(mlpTest[, -1])) {
              mlpMod <-
                mlp(
                  x = mlpTrain[, -1],
                  y = mlpTrain$obs,
                  size = c(25, 25),
                  maxit = 200
                )
              probTmp[act_gr, dimnames(mlpTrain$obs)[[2]]] <-
                predict(mlpMod, mlpTest[-1])
              missingProbs <- is.na(probTmp[act_gr, ])
              probTmp[act_gr, missingProbs] <- 0
            }
          }
          probTmp
        }

      probs[, , , lead] <- aperm(predTmp, c(3, 1, 2))
      x <- apply(probs, c(1, 2, 4), sum, na.rm = T)
      probs <- sweep(probs, c(1, 2, 4), x, "/")

    }
    stopCluster(c1)
    gc()

    return(probs)
  }
