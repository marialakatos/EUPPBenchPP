
#' Fitting POLR-L model to the EUPPBench visibility ensemble
#'
#' @name polrProbsVisLoc
#' @description Fitting POLR-L model to the EUPPBench datasets' visibility forecasts.
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
#' @param startZeta
#'
#' @return Estimated probabilities
#' @export
polrProbsVisLoc <-
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
           train = 30,
           startZeta = seq(-1.7, 6.5, by = .1)) {

    start.day <- 366
    #end.day <- 730
    end.day <- 368

    forcDates <- obsDates[start.day:end.day]
    nForcDates <- length(forcDates)

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

      start.time <- Sys.time()
      for (d in (start.day):(end.day)) {
        # print(obsDates[d])
        train.ind <- (d - train - lt + 1):(d - lt)
        data.train <- visData[train.ind, , , lead]
        data.test <- visData[d, , , lead]
        per1 <- sin(2 * pi * train.ind / 365)
        per2 <- cos(2 * pi * train.ind / 365)

        for (act_gr in c(1:nStats)) {
          train.aux <- cbind(data.train[, act_gr, ], per1, per2)
          colnames(train.aux) <- varNames
          polrTrain <- data.frame(train.aux)
          test.aux <-
            c(data.test[act_gr, ], sin(2 * pi * d / 365), cos(2 * pi * d / 365))
          names(test.aux) <- varNames
          polrTest <- data.frame(t(test.aux[-1]))
          ### Keeping coefficients of ctrl and exmean nonnegative
          coefCont <- TRUE
          posCoef <- 3
          nVals <- length(table(polrTrain$obs))
          while (coefCont) {
            err <-
              try(
                polrMod <-
                  polr(as.factor(obs) ~ ., data = polrTrain),
                silent = TRUE
              )
            if (class(err) == "try-error") {
              startPars <- c(rep(1, (6 + posCoef)), startZeta[1:(nVals - 1)])
              polrMod <-
                polr(as.factor(obs) ~ .,
                  data = polrTrain,
                  start = startPars
                )
            }
            if (posCoef > 0) {
              A <- polrMod$coefficients[1:posCoef] < 0
              if (sum(A) > 0) {
                posCoef <- posCoef - sum(A)
                badInd <- which(A)
                polrTrain <- polrTrain[, -(badInd + 1)]
              } else {
                coefCont <- FALSE
              }
            } else {
              coefCont <- FALSE
            }
          }
          probTmp <- predict(polrMod, polrTest, type = "p")
          probs[(d - start.day + 1), act_gr, , lead] <- 0
          probs[(d - start.day + 1), act_gr, names(probTmp), lead] <-
            probTmp
        }
      }
    }
    return(probs)
  }
