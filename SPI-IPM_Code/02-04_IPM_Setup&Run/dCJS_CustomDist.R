##########################################################
#### SPECIFICATION OF CUSTOM DISTRIBUTION FOR CMR DATA ###
##########################################################

dCJS_vv_sum <- nimbleFunction(
  # It is assumed that the individual has already been captured.
  # Therefore, the first entry in x represents the first possible recapture event.
  # probSurvive[t] represents survival from t-1 to t.
  # probCapture[t] represents capture probability at time t.
  run = function(x = double(1),    ## standard name for the "data"
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0), #! NEWLY ADDED: argument stating number of occurences of same capture history in entire dataset
                 len = double(0, default = 0),
                 log = integer(0, default = 0) ## required log argument
  ) {
    if (len != 0) {
      if (len != length(x)) stop("Argument len must match length of data, or be 0.")
    }
    if (length(probSurvive) < length(x) - 1)
      stop("Length of probSurvive must be at least length of data minus 1.")
    if (length(x) != length(probCapture)) stop("Length of probCapture does not match length of data.")
    if (x[1] != 1) stop("dCJS requires specifying first capture. x[1] must equal 1.")

    ## Note the calculations used here are actually in hidden Markov model form.
    probAliveGivenHistory <- 1
    ## logProbData will be the final answer
    logProbData <- 0
    if (len == 0) {  ## l<1 should not occur, but just in case:
      len <- length(x)
    }
    for (t in 2:len) {
      ## probAlive is P(Alive(t) | x(1)...x(t-1))
      ## probAliveGivenHistory is (Alive(t-1) | x(1)...x(t-1))
      probAlive <- probAliveGivenHistory * probSurvive[t - 1]
      if (!is.na(x[t])) {
        if (x[t] == 1) {
          ## ProbThisObs = P(x(t) | x(1)...x(t-1))
          probThisObs <- probAlive * probCapture[t]
          probAliveGivenHistory <- 1
        } else {
          probAliveNotSeen <- probAlive * (1 - probCapture[t])
          probThisObs <- probAliveNotSeen + (1 - probAlive)
          probAliveGivenHistory <- probAliveNotSeen / probThisObs
        }
        logProbData <- logProbData + log(probThisObs) * mult #! NEWLY ADDED: "mult"
      }
    }
    if (log) {
      return(logProbData)
    }
    return(exp(logProbData))
    returnType(double())
  }
)

rCJS_vv_sum <- nimbleFunction(
  run = function(n = integer(),
                 probSurvive = double(1),
                 probCapture = double(1),
                 mult = double(0),
                 len = double(0, default = 0)) {
    if (n != 1) stop("rCJS only works for n = 1")
    if (len < 2)
      stop("len must be greater than 1.")
    if(length(probSurvive) != len - 1)
      stop("Length of probSurvive is not the same as len - 1.")
    if(length(probCapture) != len)
      stop("Length of probCapture is not the same as len.")
    ans <- numeric(length = len, init = FALSE)
    ans[1] <- 1
    alive <- 1
    if (len <= 0) return(ans)
    for (i in 2:len) {
      if (alive)
        alive <- rbinom(1, size = 1, prob = probSurvive[i - 1])
      if (alive) {
        ans[i] <- rbinom(1, size = 1, prob = probCapture[i])
      } else {
        ans[i] <- 0
      }
    }
    return(ans)
    returnType(double(1))
  }
)

registerDistributions(list(
  dCJS_vv_sum = list(
    BUGSdist = 'dCJS_vv_sum(probSurvive, probCapture, mult, len)',
    types = c('value = double(1)', 'probSurvive = double(1)', 'probCapture = double(1)', 'mult = double(0)', 'len = double()'),
    discrete = TRUE
  )
))
