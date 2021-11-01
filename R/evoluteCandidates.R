#########################################################################
#
# Package: BBMateAllocation
#
# File: evoluteCandidates.R
# Contains: evoluteCandidates, evol
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Pick candidates to reduce relatedness and maximize criterion
#'
#' Evolutionary algorithm to pick candidates for mating plan to reduce relatedness and maximize criterion
#'
#' @param criterion data frame with variable to maximize
#' @param K relationship matrix
#' @param n.candidates number of crosses in the mating plan
#' @param steps numeric vector of steps in each evolution phase
#' @param n.try attempts to replace the current number of candidates (current step) to improve candidates list before moving for next phase
#' @param n.burnIn attempts to initialize the chain
#'
#' @return A list with final candidate list, and relationship and criterion evolution
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'

evoluteCandidates = function(criterion,
                             K,
                             n.candidates,
                             step=c(5,4,3,2,1),
                             n.steps=200,
                             n.try=10000,
                             n.burnIn=10){

  y=criterion[,2]

  cat("Burn-In...\n")
  step = c(step[1],step)
  iterations = vector("list",length(step))
  iterations[[1]] = evol(y,K,n.candidates=n.candidates,step=step[1], n.steps=round(n.steps/2), n.try=round(n.try/2))
  for(i in 2:n.burnIn){
    burnInTry = evol(y,K,n.candidates=n.candidates,step=step[1],n.steps=round(n.steps/2), n.try=round(n.try/2))
    if(tail(burnInTry$gain,1) > tail(iterations[[1]]$gain,1))
      iterations[[1]] = burnInTry
  }

  cat("Evolving...")
  for(i in 2:length(step))
    iterations[[i]] = evol(y,K,iterations[[i-1]]$candidates,n.candidates=n.candidates,step=step[i],n.steps=n.steps,n.try=n.try)

  output = data.frame(gain=unlist(lapply(iterations, "[[", 3)),
                      rel=unlist(lapply(iterations, "[[", 2)))
  return(list(candidates=rownames(K)[iterations[[length(step)]]$candidates],steps=output))
}

evol = function(y,
                K,
                candidates=NULL,
                n.candidates=200,
                step,
                n.steps,
                n.try){

  rel = gain = NULL
  if(is.null(candidates))
    candidates = sort(sample(1:length(y), n.candidates))

  pool = c(1:length(y))[-candidates]

  for(i in 1:n.steps){
    rel[i] = mean(K[candidates,candidates],na.rm=TRUE)
    gain[i] = mean(y[candidates])
    evolute = TRUE
    tries = 1
    while(evolute){
      candidatesMutate = c(candidates[-sample(1:200,step)],sample(pool,step))
      gainMutate = mean(y[candidatesMutate])
      relMutate = mean(K[candidatesMutate,candidatesMutate],na.rm=TRUE)
      if(relMutate < rel[i]){
        if(gainMutate > gain[i]){
          candidates = candidatesMutate
          evolute = FALSE
        }
      }
      tries = tries+1
      if(tries == n.try){
        #cat("reach n.tries")
        return(list(candidates=candidates,rel=rel,gain=gain))
      }
    }
  }
  return(list(candidates=candidates,rel=rel,gain=gain))
}
