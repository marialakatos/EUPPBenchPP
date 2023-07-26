#' Fitting POLR-C model to the EUPPBench visibility ensemble
#'
#' @name polrProbsVisGroup
#' @description Fitting POLR-C model to the EUPPBench datasets' visibility forecasts.
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
#' @param nForcDates Length of verification period
#' @param nStats Number of stations
#' @param nObsNames Number of observation names
#' @param nLeadTimes Number of lead times
#' @param forcDates A character vector of dates in the verification period
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
polrProbsVisGroup <-
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
           startZeta = seq(-1.7, 6.5, by = .1),
           min_stats = 4,
           nGr = 8) {

    start.day <- 366
    end.day <- 730

    forcDates <- obsDates[start.day:end.day]
    nForcDates <- length(forcDates)

    probs <- array(
      data = NA,
      dim = c(nForcDates, nStats, nObsNames, nLeadTimes),
      dimnames = list(forcDates, stationIDs, obsNames, leadTimes)
    )
    shift <- 0

    for (lead in leadTimes[1:nLeadTimes]) {
      if (as.numeric(lead) == 0) {
        lt <- 1
      } else {
        lt <- ceiling(as.numeric(lead) / 24)
      }

      catNames <- c("<5km", "5-30km", ">30km")
      nCatNames <- length(catNames)

      cat(
        "Modelling for lead time",
        lead, "hours is in progress ...\n"
      )


      for (d in (start.day):(end.day)) {
        train.ind <- (d - train - lt + 1):(d - lt)
        data.train <- visData[train.ind, , , lead]
        data.test <- visData[d, , , lead]
        statDist <-
          array(
            data = 0,
            dim = c(nStats, nCatNames),
            dimnames = list(stationIDs, catNames)
          )
        statDist[, 1] <-
          apply(data.train[, , 1] <= 5000, 2, sum, na.rm = T)
        statDist[, 2] <-
          apply(((data.train[, , 1] > 5000) &
                   (data.train[, , 1] <= 30000)), 2, sum, na.rm = T)
        statDist[, 3] <-
          apply(data.train[, , 1] > 30000, 2, sum, na.rm = T)

        statCl <- kmeans(statDist, nGr)
        nGrnew <- nGr
        while (min(table(statCl$cluster)) < min_stats) {
          nGrnew <- nGrnew - 1
          statCl <- kmeans(statDist, nGrnew)
        }
        probs[(d - start.day + 1), , , lead] <- 0
        for (gr in c(1:nGrnew)) {
          act_gr <- (statCl$cluster == gr)
          obs.train <- as.vector(data.train[, act_gr, 1])
          hres.train <- as.vector(data.train[, act_gr, 2])
          ctrl.train <- as.vector(data.train[, act_gr, 3])
          exmean.train <- as.vector(data.train[, act_gr, 4])
          sd.train <- as.vector(data.train[, act_gr, 5])
          p1.train <- as.vector(data.train[, act_gr, 6])
          p2.train <- as.vector(data.train[, act_gr, 7])
          p70.train <- as.vector(data.train[, act_gr, 8])

          per1 <- rep(sin(2 * pi * train.ind / 365 + shift), sum(act_gr))
          per2 <- rep(cos(2 * pi * train.ind / 365 + shift), sum(act_gr))

          polrTrain <-
            data.frame(
              cbind(
                obs.train,
                hres.train,
                ctrl.train,
                exmean.train,
                sd.train,
                p1.train,
                p2.train,
                p70.train,
                per1,
                per2
              )
            )
          dimnames(polrTrain)[[2]] <- varNames

          test.aux <-
            cbind(data.test[act_gr, -1], rep(sin(2 * pi * d / 365 + shift), sum(act_gr)), rep(cos(2 *
                                                                                                    pi * d / 365 + shift), sum(act_gr)))
          colnames(test.aux) <- varNames[-1]
          polrTest <- data.frame(test.aux)

          ### Keeping coefficients of hres, ctrl and exmean nonnegative
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
          ###
          probTmp <- predict(polrMod, polrTest, type = "p")
          probs[(d - start.day + 1), act_gr, dimnames(probTmp)[[2]], lead] <-
            probTmp
        }
      }

    }
    return(probs)
  }
