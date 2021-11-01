#########################################################################
#
# Package: BBMateAllocation
#
# File: mutatePlan.R
# Contains: mutatePlan
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Evolute mating plan to reduce relatedness
#'
#' Evolute mating plan to reduce relatedness
#'
#' @param initialPlan initial plan to begin to mutate
#' @param data data frame with possible crosses to filter out
#' @param K relationship matrix
#' @param rate mutation rate (how many lines to substitute)
#' @param memory.mutation after how many successes to reset the mutation memory
#' @param tries maximum number of attempts to mutate
#' @param perc.reduce.Y fraction allowed to reduce Y during mutation
#' @param n.cross number of crosses in the mating plan
#' @param max.cross maximum number of crosses per individual in the plan
#' @param min.cross minimum number of crosses per individual in the plan
#' @param max.cross.to.search maximum number of rows in the input to search for the solution
#' @param culling.pairwise.k don't allow crosses with more than this threshold of relatedness. Default is 0.25.
#'
#' @return A list with a data frame with the steps and a data frame with the mating plan at the last iteration
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'

mutatePlan = function(initialPlan,
                      data,
                      K,
                      rate = 10,
                      memory.mutation = 3,
                      tries = 1000,
                      perc.reduce.Y = 0,
                      n.cross = 100,
                      max.cross = 4,
                      min.cross = 2,
                      max.cross.to.search=100000,
                      culling.pairwise.k = 0.0675){

  step=NULL
  current.Y = initialPlan[[1]]$target.Y
  current.K = initialPlan[[1]]$target.K
  currentPlan = initialPlan
  mutations = NULL
  memory.mutation.track=0
  for(i in 1:tries){
    allCrossesTry = data
    tryMutate = c(mutations,sample(as.numeric(rownames(currentPlan$plan)),ifelse(length(rate)==1,rate,sample(rate,1))))
    allCrossesTry[tryMutate,3] = 0

    ## Mutate candidates?
    allCrossesTry[sample(1:nrow(allCrossesTry),round(nrow(allCrossesTry)/5)),3] = 0

    tryPlan = selectCrosses(data=allCrossesTry,
                            K=Amat,
                            n.cross=n.cross,
                            max.cross=max.cross,
                            min.cross=min.cross,
                            max.cross.to.search=max.cross.to.search,
                            culling.pairwise.k = culling.pairwise.k)
    if((tryPlan[[1]]$target.Y > current.Y*(1-perc.reduce.Y)) & (tryPlan[[1]]$target.K < current.K)){
      allCrossesFull= allCrossesTry
      step = rbind(step,data.frame(tryPlan[[1]],step=i))
      currentPlan = tryPlan
      current.K = tryPlan[[1]]$target.K
      current.Y = tryPlan[[1]]$target.Y
      if(nrow(step)==1){
        write.table(tail(step,1), col.names = TRUE, quote = FALSE)
      }else{
        write.table(tail(step,1), col.names = FALSE, quote = FALSE)
      }
      memory.mutation.track = memory.mutation.track+1
      if(memory.mutation.track == memory.mutation){
        mutations = NULL
        memory.mutation.track = 0
        cat("erasing mutation memory...\n")
      }else{
        mutations = tryMutate
      }
    }
    if(i%%100==0)
      cat(paste0(i,"/",tries,"\n"))
  }
  return(list(plan=currentPlan,step=step))
}
