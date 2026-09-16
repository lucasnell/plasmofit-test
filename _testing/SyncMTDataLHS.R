


####################################
# Code to fit PMRs from human data #
####################################

# This dataset is shared under an Attribution 4.0 International (CC BY 4.0)
# license (https://creativecommons.org/licenses/by/4.0/). The material can
# be shared and built upon, but attribution to the original authors and a
# statement of changes made is required.

options(warn=1)

library(deSolve)
library(parallel)
library(lhs)

parmIter = 1e3
optimizationRuns = 100
if(!file.exists(paste("sampleLHS",parmIter,".txt",sep = ''))){
  sampleLHS = maximinLHS(parmIter,4)
  write.table(sampleLHS, file = paste("sampleLHS",parmIter,".txt",sep = ''))}
if(file.exists(paste("sampleLHS",parmIter,".txt",sep = ''))){
  sampleLHS = read.table(paste("sampleLHS",parmIter,".txt",sep = ''))}

# calculate PMRs
## @knitr pmrfx
pmrfx=function(timeSeries, alpha){
  # function to calculate parasite multiplication rates

  # 1. build a matrix to store results
  pmrRes = as.data.frame(matrix(NA, nrow = length(timeSeries$day), ncol = 3))
  colnames(pmrRes) <- c("day", "pmr.obs", "para.obs")

  # 2a. store all time points encompassed in experiment
  pmrRes$day = timeSeries$day
  # 2b. store all parasite counts recorded
  pmrRes$para.obs = timeSeries$Asex

  # 3. calculate parasite multiplication rate (PMR) from observed time series
  dayPlusAlpha = timeSeries$day + alpha
  index.t = which(signif(24*dayPlusAlpha, digits = 6) %in% signif(24*timeSeries$day, digits = 6))
  index.tPlusAlpha = which(signif(24*timeSeries$day, digits = 6) %in% signif(24*dayPlusAlpha, digits = 6))
  obsPMR = timeSeries$Asex[index.tPlusAlpha]/timeSeries$Asex[index.t]
  pmrRes$pmr.obs[index.t] = obsPMR
  # only return rows that correspond to entries in time series
  toReturn = pmrRes

  return(toReturn)
} # end pmrfx function

## @knitr findPeakInfection
# timing of peak: use first local maximum that is
# preceded and followed by 6 consecutive lower values
# after Eichner et al. 2001
# all waves (modified from Eichner)
findPeakInfection = function(infData){
  vecVal = approx(x = infData$day[is.finite(infData$Asex)],
                  y = infData$Asex[is.finite(infData$Asex)],
                  xout = seq(min(infData$day), max(infData$day),by=1))$y
  # first locate candidate local maxima:
  candidateLMLocalOnly = c(which(diff(sign(diff(vecVal))) == -2) + 1,length(vecVal))
  candidateLM = candidateLMLocalOnly
  if(!which.max(vecVal)%in%candidateLMLocalOnly){
    candidateLM = sort(c(which.max(vecVal), candidateLMLocalOnly))}
  isLM = rep(F, length(candidateLM))
  for (lmIndex in 1:length(candidateLM)){
    lmVal = candidateLM[lmIndex]
    startIndex = lmVal - 6; if(startIndex < 1){startIndex = 1}
    endIndex = lmVal + 6; if(endIndex > length(vecVal)){endIndex = length(vecVal)}
    if(all(c(vecVal[startIndex:(lmVal-1)],vecVal[(lmVal+1):endIndex]) < vecVal[lmVal],
           na.rm = T)){isLM[lmIndex] = T} # end conditional: are all six values before and after ostensible peak lower than ostensible peak value? (or missing)
  } # end loop through candidate peak values
  eichnerLM = candidateLM[which(isLM)]
  return(eichnerLM[1])
} # end findPeakInfection fx

# load MT data
asexLimit = 10 # detection limit based on Simpson et al. 2002 Parasitology
# read in table of days listing when drugs were administered:
drugdays=read.csv("DrugDays.csv", header=T)
# read in table of circulating iRBC counts
asex=read.csv("Asexual.csv",header=F)
# add in strain & route info
strainInfo = read.csv("StrainALL.csv", header=F)
routeInfo = read.csv("RouteALL.csv", header=F)

## @knitr betaOffsetStart
# betaOffsetStart returns a vector of the initial iRBC cohort
#   sorted into age compartments according to a Beta distribution
#   where the median age of infection is adjusted from the default
#   (halfway through IED) according to the offset parameter that can
#   take on any value between 0 and 1 (note the function uses the CDF
#   to sort the initial cohort into compartments)
betaOffsetStart = function(compartments, shape, offset, initialI){
  times = seq(0, 1, length.out = round(compartments + 1)) + offset
  correctedTimes = times
  correctedTimes[times > 1] = times[times > 1] - 1
  pbetaValues = pbeta(correctedTimes, shape1 = shape, shape2 = shape)
  correctedpbetaValues = pbetaValues
  correctedpbetaValues[times > 1] = pbetaValues[times > 1] + 1
  iVector = initialI * diff(correctedpbetaValues)
  if(round(sum(iVector)) != round(initialI)){
    warning(paste("betaOffsetStart magnitude error (",
                  round(sum(iVector)) - round(initialI), ", offset = ",
                  round(offset, digits = 2), ", shape = ",
                  round(shape, digits = 2), sep = ''), immediate. = T)}
  if(round(length(iVector)) != round(compartments)){
    warning(paste("betaOffsetStart length error, offset = ",
                  round(offset, digits = 2), ", shape = ",
                  shape, sep=''), immediate. = T)}
  return(iVector)
}

## @knitr offsetToPeak
# function to convert offset value to median parasite age
# (as fraction through IED)
# multiply by cycle length in hours to get age in hours
offsetToPeak = function(offset){
  res = rep(NA, length(offset))
  for (i in 1:length(offset)){
    res[i] = 1.5 - offset[i]
    if(offset[i]<0.5){res[i] = 0.5 - offset[i]}
  }
  return(res)
}

## @knitr burstSpan
burstSpan=function(shape,percentBurst=0.99,cycleLength=1){
  # convenience function to take fitted beta shape parameter
  # and return the hours required for (by default) 99% of bursting
  return(diff(24 * cycleLength * qbeta(c((1 - percentBurst) / 2,
                                         1 - (1 - percentBurst) / 2),
                                       shape, shape)))
  }

## @knitr allowableInitSync
# get burstSpans associated with possible sP values
possiblesP = 1:3e5
burstSpanValues = unlist(lapply(as.list(possiblesP),burstSpan, cycleLength=2))
maxAllowedsP = approxfun(x=burstSpanValues, y=possiblesP, rule=2)

## @knitr gfx
# function parameterized from Kriek et al. 2003
# for sequestration as a function of hpi
gfx = function(age, p1 = 11.3869/467.6209,
               p2 = 1,
               p3 = 18.5802,
               p4 = 0.2242){
  gVal = p1+((p2-p1)/(1+10^(p4*(p3-age))))
  return(gVal)
} # end gfx function

## @knitr constantPMR
constantPMR.gammaN=function(t, y, parms){
  # Archer model with constant PMR, variable number of age compartments
  cycleLength = parms[1] # in hours
  mu = parms[2] # mortality of non-sequestered iRBCs
  museq = parms[3] # mortality of sequestered iRBCs
  R = parms[4] # parasite multiplication rate (constant)
  n = parms[5] # number of compartments

  lambdaN = n/cycleLength
  lambdaS = lambdaN

  # generate vector of ages through IED:
  ageValues = seq(1/lambdaN,cycleLength,by=1/lambdaN)
  gValues = gfx(ageValues)
  yValues = 1-gValues
  qValues = c(0, 1 - yValues[2:length(yValues)]/yValues[1:(length(yValues)-1)])

  N = y[1:n]
  S = y[(n+1):(2*n)]

  dN = rep(0, n)
  dS = rep(0, n)
  dN[1] = R*(lambdaN * N[n] + lambdaS * S[n]) - (lambdaN + mu) * N[1]
  # Correction from 5/30/24 meeting was
  # dS[1] = - (lambdaS + mu) * S[1]
  # but I believe it should be:
  dS[1] = - (lambdaS + museq) * S[1]
  dN[2:n] = (1-qValues[1:(n-1)]) * lambdaN * N[1:(n-1)] - (lambdaN + mu) * N[2:n]
  dS[2:n] = qValues[1:(n-1)] * lambdaN * N[1:(n-1)] + lambdaS * S[1:(n-1)] - (lambdaS + museq) * S[2:n]

  res=c(dN,dS)
  list(res)
} # end simplified gamma chain model

## @knitr archer.fixedN
# parms to fit: sP, offset, R, I0
archer.fixedN = function(parms, data, optimization = T, n){
  # return predicted circulating iRBC abundance
  # where data represents the observed circulating iRBC abundance
  # fit variation in beta distribution, and convert to sP
  # where sP = (1/8)*((1/var)-4) and var ranges from zero
  # to 1/12 (1/12 corresponding to a uniform distribution)
  varBetaDist = exp(-exp(parms[1]))*(1/12)
  betaShape = (1/8)*((1/varBetaDist)-4) # must be positive,
  # need not be an integer
  # set limit because larger shape parameters become
  # computationally intensive
  pfCycleLength = 48 # hours
  ageClassDuration = pfCycleLength/n
  sPLimit = maxAllowedsP(ageClassDuration) # set upper limit on sP to max allowed
  if(betaShape > sPLimit){betaShape = sPLimit}
  offset = exp(-exp(parms[2])) # so that mean bursting
  # need not fall at 0.5
  R = exp(parms[3])
  # assume that initial abundance is equal to the sum of the first two counts
  I0 = exp(parms[4])

  if(is.na(optimization)){sse = c(betaShape, offset, R, I0, n)}
  if(is.na(optimization) == F){
    mu = 0
    museq = 0
    parmValues = c(pfCycleLength, mu, museq, R, n)

    ages = seq(pfCycleLength/n,pfCycleLength,by=pfCycleLength/n)
    gs = gfx(ages)
    ys = 1-gs
    qs = c(0, 1 - ys[2:length(ys)]/ys[1:(length(ys)-1)])

    # choose offset value equivalent to initial
    startI0All = betaOffsetStart(n,betaShape,offset,I0)
    startI0 = c((1-gs)*startI0All,gs*startI0All)

    timeValues = 24*data$day # in hours
    archer.output = as.data.frame(lsoda(startI0, times = timeValues,
                                        func = constantPMR.gammaN,
                                        parms = parmValues))
    # calculate circulating iRBC abundance
    circ.iRBC = apply(archer.output[ , 2:(n + 1)],1,sum)
    # calculate sequestered iRBC abundance
    seq.iRBC = apply(archer.output[ , (n + 2):(n*2 + 1)],1,sum)

    if(optimization == F){
      res = data
      res$predCirc = circ.iRBC
      res$predTotal = circ.iRBC + seq.iRBC
      sse = res}
    if(optimization == T){
      # use sum squared error (SSE)
      sse = sum(((log10(data$Asex+1) - log10(circ.iRBC+1))^2), na.rm = T)
      if(is.na(sse) | sse > 1e10){sse = 1e10}}
  } # end conditional: optimization != to NA
  return(sse)
} # end function: archer.fit

## @knitr archer.fitN
# parms to fit: sP, offset, R, I0, n
archer.fitN = function(parms, data, optimization = T, nLimit){
  # return predicted circulating iRBC abundance
  # where data represents the observed circulating iRBC abundance
  cv.cycleLength = exp(-exp(parms[5]))
  n = round(1/(cv.cycleLength^2)) # must be an integer
  # restrict n to be nLimit or less for computational efficiency
  if(n > nLimit){n = nLimit}
  # ensure n is at least 4
  if(n < 4){n = 4}
  # fit variation in beta distribution, and convert to sP
  # where sP = (1/8)*((1/var)-4) and var ranges from zero
  # to 1/12 (1/12 corresponding to a uniform distribution)
  varBetaDist = exp(-exp(parms[1]))*(1/12)
  betaShape = (1/8)*((1/varBetaDist)-4) # must be positive,
  # need not be an integer
  # set limit because larger shape parameters become
  # computationally intensive
  pfCycleLength = 48 # hours
  ageClassDuration = pfCycleLength/n
  sPLimit = maxAllowedsP(ageClassDuration) # set upper limit on sP to max allowed
  if(betaShape > sPLimit){betaShape = sPLimit}
  offset = exp(-exp(parms[2])) # so that mean bursting
  # need not fall at 0.5
  R = exp(parms[3])
  # assume that initial abundance is equal to the sum of the first two counts
  I0 = exp(parms[4])

  if(is.na(optimization)){sse = c(betaShape, offset, R, I0, n)}
  if(is.na(optimization) == F){
    mu = 0
    museq = 0
    parmValues = c(pfCycleLength, mu, museq, R, n)

    ages = seq(pfCycleLength/n,pfCycleLength,by=pfCycleLength/n)
    gs = gfx(ages)
    ys = 1-gs
    qs = c(0, 1 - ys[2:length(ys)]/ys[1:(length(ys)-1)])

    # choose offset value equivalent to initial
    startI0All = betaOffsetStart(n,betaShape,offset,I0)
    startI0 = c((1-gs)*startI0All,gs*startI0All)

    timeValues = 24*data$day # in hours
    archer.output = as.data.frame(lsoda(startI0, times = timeValues,
                                        func = constantPMR.gammaN,
                                        parms = parmValues))
    # calculate circulating iRBC abundance
    circ.iRBC = apply(archer.output[ , 2:(n + 1)],1,sum)
    # calculate sequestered iRBC abundance
    seq.iRBC = apply(archer.output[ , (n + 2):(n*2 + 1)],1,sum)

    if(optimization == F){
      res = data
      res$predCirc = circ.iRBC
      res$predTotal = circ.iRBC + seq.iRBC
      sse = res}
    if(optimization == T){
      # use sum squared error (SSE)
      sse = sum(((log10(data$Asex+1) - log10(circ.iRBC+1))^2), na.rm = T)
      if(is.na(sse) | sse > 1e10){sse = 1e10}}
  } # end conditional: optimization != to NA
  return(sse)
} # end function: archer.fit

## @knitr archerNfx
archerNfx = function(startGuess, dataForFitting, nValue){
  fit <- optim(startGuess, fn = archer.fixedN,
               data = dataForFitting, n = nValue,
               control = list(trace = 4, maxit = 5000))
  res = c(fit$par, fit$convergence, fit$value)
  return(res)
  } # end archerNfx (optimization function)

## @knitr summaryMetrics
getInfectionMetrics = function(caseIndex, returnDataToFit = F){
  workingData=as.data.frame(cbind(1:480, asex[, caseIndex]))
  colnames(workingData)<-c("day", "Asex")
  lastPositiveDay = max(which(workingData$Asex>0))
  dayPeak = findPeakInfection(workingData[1:lastPositiveDay,])
  peakHeight = workingData$Asex[dayPeak]

  toFit = T
  # determine strain & inoculation route
  # will not fit to time series from infections
  # generated by strains that are not El Limon, McLendon, or Santee Cooper
  strainValue = strainInfo[caseIndex,]
  routeValue = routeInfo[caseIndex,]
  if(!strainValue %in% c("El Limon", "McLendon", "Santee Cooper")){toFit = F}

  # record first day of drug treatment
  firstDrugDay = NA
  # exclude any patients treated in first week
  if(any(!is.na(unlist(drugdays[caseIndex,])))){
    firstDrugDay = min(unlist(drugdays[caseIndex,]), na.rm=T)
    if(firstDrugDay<=7){toFit = F}
    } # end conditional are any drugdays listed for this patient

  # check that first seven days have nonmissing counts
  # and peak is no earlier than day 6
  # (only if infection was not treated in first 7 days)
  # allow fits if some counts are below presumed detection limit
  # but record minimum during first week to allow these to be
  # excluded later
  minFirstWeek = NA
  if(toFit){
    dataIndex = 1:7
    dataToFit = data.frame(day = workingData$day[dataIndex],
                           Asex = workingData$Asex[dataIndex])
    asexCounts = workingData$Asex[dataIndex]
    minFirstWeek = min(dataToFit$Asex, na.rm = T)
    if(any(is.na(asexCounts))){toFit = F}
    if(dayPeak < 6){toFit = F}
    } # end conditional: is the peak of infection on day six or after,
      # and are the first seven days' of data nonmissing?

  metrics = c(caseIndex, strainValue, routeValue, firstDrugDay, lastPositiveDay,
              dayPeak, peakHeight, minFirstWeek, toFit)
  if(returnDataToFit==T){metrics = dataToFit}
  return(metrics)
  } # end function: get summary statistics for each infection

allSummaryInfo = lapply(as.list(1:322), getInfectionMetrics)

summaryInfoMatrix = as.data.frame(matrix(unlist(allSummaryInfo), nrow = 322,
                           ncol = round(length(unlist(allSummaryInfo)) / 322),
                           byrow = T), stringsAsFactors = F)

colnames(summaryInfoMatrix) <- c("id", "strain", "route", "firstDrugDay", "lastPositiveDay",
                           "peakDay", "peakHeight", "minFirstWeek", "toFit")

summaryInfo = data.frame(id = as.numeric(summaryInfoMatrix$id),
                         strain = summaryInfoMatrix$strain,
                         route = summaryInfoMatrix$route,
                         firstDrugDay = as.numeric(summaryInfoMatrix$firstDrugDay),
                         lastPositiveDay = as.numeric(summaryInfoMatrix$lastPositiveDay),
                         peakDay = as.numeric(summaryInfoMatrix$peakDay),
                         peakHeight = as.numeric(summaryInfoMatrix$peakHeight),
                         minFirstWeek = as.numeric(summaryInfoMatrix$minFirstWeek),
                         toFit = summaryInfoMatrix$toFit, stringsAsFactors = F)

## @knitr latinHypercube
# function to convert calculate starting guesses from latin hypercube
convertParms = function(lhsObject, initAbund, transformParms = T){
  resSampled = as.data.frame(matrix(unlist(lapply(1:length(lhsObject[1,]), function(x)
    lhsObject[,x]*diff(parmRange[,x]) + parmRange[1,x])),
    nrow = length(lhsObject[,1]), ncol = length(lhsObject[1,])))
  colnames(resSampled) <- c("varBetaDist", "offset", "R", "I0")
  resSampled$I0 = 10^(log10(initAbund)+resSampled$I0)
  res = resSampled
  if(transformParms){
    res = resSampled
    res[,1:2] = log(-log(resSampled[,1:2]))
    res[,3:4] = log(resSampled[,3:4])
  }
  return(res)
}

## @ knitr parallelComputations
infectionsToFit = summaryInfo$id[summaryInfo$toFit=="TRUE"]

optimizationTable = expand.grid(id = infectionsToFit, nValue = c(96, 192, 288))

rowIndex = as.integer(commandArgs(trailingOnly = TRUE)[1])
if(is.na(rowIndex)){rowIndex = 1}
caseIndex = optimizationTable$id[rowIndex]
nValue = optimizationTable$nValue[rowIndex]

toFit = F
if(caseIndex %in% infectionsToFit){toFit = T}

fileTag = paste("MT", caseIndex, "N", nValue, ".txt", sep = '')

if(toFit){
  dataToFit = getInfectionMetrics(caseIndex, returnDataToFit = T)

  # create Latin hypercube to sample value of objective function across range
  # of four parameters to be optimized: 1. varBetaDist; 2. offset; 3. R; 4. I0
  # repeat sampling for three n values encompassing different variance in
    # developmental duration:
      # n = 96,  var = 24 hours
      # n = 192, var = 12 hours
      # n = 288, var = 8 hours

  parmRange = data.frame(varBetaDist = c(0,1),
                         offset = c(0,1),
                         R = c(0,32),
                         I0 = c(-1,1))
  row.names(parmRange) <- c("min", "max")

  startParms = convertParms(sampleLHS, initAbund = sum(dataToFit$Asex[1:2]))

  startParmsTranformed = as.data.frame(t(apply(startParms, 1,
                                               archer.fixedN, data = NA,
                                               optimization = NA, n = nValue)))
  evalAtStartParms = apply(startParms, 1, archer.fixedN, data = dataToFit, n = nValue)
  startParmsTranformed$sse = evalAtStartParms
  write.table(startParmsTranformed, file = paste("startEvaluations",fileTag, sep = ''))
  bestParmIndex = (order(evalAtStartParms))[1:optimizationRuns]

  cl <- makeCluster(80)

  clusterExport(cl, varlist = c("archerNfx", "archer.fixedN", "constantPMR.gammaN",
                                "gfx", "betaOffsetStart", "burstSpanValues",
                                "possiblesP", "maxAllowedsP","nValue"))
  clusterCall(cl, library, package = "deSolve", character.only = TRUE)

  detest = clusterApplyLB(cl, as.list(as.data.frame(t(startParms[bestParmIndex,]))),
                          archerNfx, dataForFitting = dataToFit, nValue = nValue)

  allResults = as.data.frame(matrix(unlist(detest), nrow = optimizationRuns,
                                    ncol = round(length(unlist(detest)) / optimizationRuns),
                                    byrow = T))

  colnames(allResults) <- c("p1", "p2", "p3", "p4", "convergence", "sse")

  write.table(allResults, file = paste("allResultsWinnow", fileTag, sep = ''))

  on.exit(stopCluster(cl))

  } # end conditional: dataToFit?
