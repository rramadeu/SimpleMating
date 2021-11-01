#########################################################################
#
# Package: BBMateAllocation
#
# File: devolutePlan.R
# Contains: devolutePlan, devolute
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Devolute mating plan to reduce relatedness
#'
#' Devolute mating plan to reduce relatedness
#'
#' @param data data frame with possible crosses to filter out
#' @param K relationship matrix
#' @param n.cross number of crosses in the mating plan
#' @param max.cross maximum number of crosses per individual in the plan
#' @param min.cross minimum number of crosses per individual in the plan
#' @param max.cross.to.search maximum number of rows in the input to search for the solution
#' @param culling.pairwise.k don't allow crosses with more than this threshold of relatedness. Default is 0.25.
#' @param n.threads number of threads to run in parallel
#' @param n.steps number of steps in the evolution
#' @param min.target.k if target K is lower than this threshold, the algorithm stops.
#'
#' @return A list with a data frame with the steps and a data frame with the mating plan at the last iteration
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'
#' @import parallel dplyr

devolutePlan =  function(data,
                         K,
                         n.cross=200,
                         max.cross=4,
                         min.cross=2,
                         max.cross.to.search=100000,
                         culling.pairwise.k=0.25,
                         n.threads=1,
                         n.steps=200,
                         min.target.k = NULL){

  step = NULL

  for(i in 1:(n.steps+1)){
    if(i>1){
      if(i==2){
        write.table(tail(step,1), col.names = TRUE, quote = FALSE)
      }else{
        write.table(tail(step,1), col.names = FALSE, quote = FALSE)
      }
      data = data %>% filter(P1 != parent2rm, P2 != parent2rm)
    }

    maxGainNow = selectCrosses(data=data,
                               K=K,
                               n.cross=n.cross,
                               max.cross=max.cross,
                               min.cross=min.cross,
                               max.cross.to.search=max.cross.to.search,
                               culling.pairwise.k=culling.pairwise.k)

    if(i>n.steps)
      break

    if(!is.null(min.target.k))
      if(out$target.K <= min.target.k)
        return(list(steps=step, plan=maxGainNow))

    parents = unique(c(maxGainNow[[2]]$P1,maxGainNow[[2]]$P2))
    cl <- makeCluster(n.threads)
    clusterExport(cl=cl, varlist = NULL)
    x=1:length(parents)
    out <- parLapply(cl, x, devolute,
                     parents=parents,
                     K=K,
                     data=data,
                     n.cross=n.cross,
                     min.cross=min.cross,
                     max.cross=max.cross,
                     max.cross.to.search=max.cross.to.search,
                     culling.pairwise.k=culling.pairwise.k)

    stopCluster(cl)
    out = do.call(rbind, out)
    out = data.frame(index=x,out)
    which2rm = which(out$target.K < maxGainNow[[1]]$target.K)
    if(length(which2rm)>0){
      out = out %>% filter(target.K <  maxGainNow[[1]]$target.K)
    }else{
      out = out %>% filter(target.K < sort(out$target.K)[round(nrow(out)/10)])
    }
    parent2rm = parents[out$index[which.max(out$target.Y)]]
    step = rbind(step,data.frame(parent=parent2rm,out[which.max(out$target.Y),2:4]))
  }
  return(list(steps=step, plan=maxGainNow))
}

devolute = function(x, parents, K, data, n.cross, min.cross, max.cross, max.cross.to.search, culling.pairwise.k){
  library(dplyr)
  data = data %>% filter(P1 != parents[x], P2 != parents[x])
  return(selectCrosses(data=data,
                       K=K,
                       n.cross=n.cross,
                       max.cross=max.cross,
                       min.cross=min.cross,
                       max.cross.to.search=max.cross.to.search,
                       culling.pairwise.k=culling.pairwise.k)[[1]])
}
