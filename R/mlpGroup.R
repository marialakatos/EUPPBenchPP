require(doParallel)

#' Fitting MLP-C model to the EUPPBench visibility ensemble
#'
#' @name mlpProbsVisGroup
#' @description Fitting MLP-C model to the EUPPBench datasets' visibility forecasts.
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
mlpProbsVisGroup <-
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
           nGr = 8,
           min_stats = 4) {


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

    catNames <- c("<5km", "5-30km", ">30km")
    nCatNames <- length(catNames)

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

          probTmp <-
            array(
              data = NA,
              dim = c(nStats, nObsNames),
              dimnames = list(stationIDs, obsNames)
            )

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
            per1 <- rep(sin(2 * pi * train.ind / 365), sum(act_gr))
            per2 <- rep(cos(2 * pi * train.ind / 365), sum(act_gr))

            mlpTrain <-
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
            dimnames(mlpTrain)[[2]] <- varNames
            mlpTrain <- mlpTrain[complete.cases(mlpTrain), ]
            mlpTrain$obs <- decodeClassLabels(mlpTrain$obs)

            test.aux <-
              cbind(data.test[act_gr, ], rep(sin(2 * pi * d / 365), sum(act_gr)), rep(cos(2 *
                                                                                            pi * d / 365), sum(act_gr)))
            colnames(test.aux) <- varNames
            mlpTest <- data.frame(test.aux)
            mlpTest$obs <- decodeClassLabels(mlpTest$obs)

            if (anyNA(mlpTest[, -1])) {
              probTmp <- array(data = NA, dim = c(nStats, nObsNames))
            } else {
              mlpMod <-
                mlp(
                  x = mlpTrain[, -1],
                  y = mlpTrain$obs,
                  size = c(25, 25),
                  maxit = 200
                )
              probTmp[act_gr, dimnames(mlpTrain$obs)[[2]]] <-
                predict(mlpMod, mlpTest[, -1])
            }
          }
          missingProbs <- is.na(probTmp)
          probTmp[missingProbs] <- 0
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
