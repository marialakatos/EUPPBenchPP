acomb <- function(...) {
  abind(..., along = 3)
}

#' Recoding raw visibility ensemble forecasts according to the reported observation
#' categories of the WMO
#'
#' @param ens Raw ensemble forecasts
#'
#' @return Recoded ensemble forecast

recode.wmo <- function(ens) {
  ens[ens < 100] <- 0
  ens[(ens >= 100) & (ens < 200)] <- 100
  ens[(ens >= 200) & (ens < 300)] <- 200
  ens[(ens >= 300) & (ens < 400)] <- 300
  ens[(ens >= 400) & (ens < 500)] <- 400
  ens[(ens >= 500) & (ens < 600)] <- 500
  ens[(ens >= 600) & (ens < 700)] <- 600
  ens[(ens >= 700) & (ens < 800)] <- 700
  ens[(ens >= 800) & (ens < 900)] <- 800
  ens[(ens >= 900) & (ens < 1000)] <- 900
  ens[(ens >= 1000) & (ens < 1100)] <- 1000
  ens[(ens >= 1100) & (ens < 1200)] <- 1100
  ens[(ens >= 1200) & (ens < 1300)] <- 1200
  ens[(ens >= 1300) & (ens < 1400)] <- 1300
  ens[(ens >= 1400) & (ens < 1500)] <- 1400
  ens[(ens >= 1500) & (ens < 1600)] <- 1500
  ens[(ens >= 1600) & (ens < 1700)] <- 1600
  ens[(ens >= 1700) & (ens < 1800)] <- 1700
  ens[(ens >= 1800) & (ens < 1900)] <- 1800
  ens[(ens >= 1900) & (ens < 2000)] <- 1900
  ens[(ens >= 2000) & (ens < 2100)] <- 2000
  ens[(ens >= 2100) & (ens < 2200)] <- 2100
  ens[(ens >= 2200) & (ens < 2300)] <- 2200
  ens[(ens >= 2300) & (ens < 2400)] <- 2300
  ens[(ens >= 2400) & (ens < 2500)] <- 2400
  ens[(ens >= 2500) & (ens < 2600)] <- 2500
  ens[(ens >= 2600) & (ens < 2700)] <- 2600
  ens[(ens >= 2700) & (ens < 2800)] <- 2700
  ens[(ens >= 2800) & (ens < 2900)] <- 2800
  ens[(ens >= 2900) & (ens < 3000)] <- 2900
  ens[(ens >= 3000) & (ens < 3100)] <- 3000
  ens[(ens >= 3100) & (ens < 3200)] <- 3100
  ens[(ens >= 3200) & (ens < 3300)] <- 3200
  ens[(ens >= 3300) & (ens < 3400)] <- 3300
  ens[(ens >= 3400) & (ens < 3500)] <- 3400
  ens[(ens >= 3500) & (ens < 3600)] <- 3500
  ens[(ens >= 3600) & (ens < 3700)] <- 3600
  ens[(ens >= 3700) & (ens < 3800)] <- 3700
  ens[(ens >= 3800) & (ens < 3900)] <- 3800
  ens[(ens >= 3900) & (ens < 4000)] <- 3900
  ens[(ens >= 4000) & (ens < 4100)] <- 4000
  ens[(ens >= 4100) & (ens < 4200)] <- 4100
  ens[(ens >= 4200) & (ens < 4300)] <- 4200
  ens[(ens >= 4300) & (ens < 4400)] <- 4300
  ens[(ens >= 4400) & (ens < 4500)] <- 4400
  ens[(ens >= 4500) & (ens < 4600)] <- 4500
  ens[(ens >= 4600) & (ens < 4700)] <- 4600
  ens[(ens >= 4700) & (ens < 4800)] <- 4700
  ens[(ens >= 4800) & (ens < 4900)] <- 4800
  ens[(ens >= 4900) & (ens < 5000)] <- 4900
  ens[(ens >= 5000) & (ens < 6000)] <- 5000
  ens[(ens >= 6000) & (ens < 7000)] <- 6000
  ens[(ens >= 7000) & (ens < 8000)] <- 7000
  ens[(ens >= 8000) & (ens < 9000)] <- 8000
  ens[(ens >= 9000) & (ens < 10000)] <- 9000
  ens[(ens >= 10000) & (ens < 11000)] <- 10000
  ens[(ens >= 11000) & (ens < 12000)] <- 11000
  ens[(ens >= 12000) & (ens < 13000)] <- 12000
  ens[(ens >= 13000) & (ens < 14000)] <- 13000
  ens[(ens >= 14000) & (ens < 15000)] <- 14000
  ens[(ens >= 15000) & (ens < 16000)] <- 15000
  ens[(ens >= 16000) & (ens < 17000)] <- 16000
  ens[(ens >= 17000) & (ens < 18000)] <- 17000
  ens[(ens >= 18000) & (ens < 19000)] <- 18000
  ens[(ens >= 19000) & (ens < 20000)] <- 19000
  ens[(ens >= 20000) & (ens < 21000)] <- 20000
  ens[(ens >= 21000) & (ens < 22000)] <- 21000
  ens[(ens >= 22000) & (ens < 23000)] <- 22000
  ens[(ens >= 23000) & (ens < 24000)] <- 23000
  ens[(ens >= 24000) & (ens < 25000)] <- 24000
  ens[(ens >= 25000) & (ens < 26000)] <- 25000
  ens[(ens >= 26000) & (ens < 27000)] <- 26000
  ens[(ens >= 27000) & (ens < 28000)] <- 27000
  ens[(ens >= 28000) & (ens < 29000)] <- 28000
  ens[(ens >= 29000) & (ens < 30000)] <- 29000
  ens[(ens >= 30000) & (ens < 35000)] <- 30000
  ens[(ens >= 35000) & (ens < 40000)] <- 35000
  ens[(ens >= 40000) & (ens < 45000)] <- 40000
  ens[(ens >= 45000) & (ens < 50000)] <- 45000
  ens[(ens >= 50000) & (ens < 55000)] <- 50000
  ens[(ens >= 55000) & (ens < 60000)] <- 55000
  ens[(ens >= 60000) & (ens < 65000)] <- 60000
  ens[(ens >= 65000) & (ens < 70000)] <- 65000
  ens[(ens >= 70000)] <- 70000

  return(ens)
}
