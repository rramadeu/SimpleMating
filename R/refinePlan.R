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
#' @param plan current mating plan
#' @param K relationship matrix
#' @param n.cross number of crosses in the mating plan
#' @param max.cross maximum number of crosses per individual in the plan
#' @param min.cross minimum number of crosses per individual in the plan
#' @param culling.pairwise.k don't allow crosses with more than this threshold of relatedness. Default is 0.25.
#' @param max.attempts maximum number of attempts to improve the plan
#'
#' @return A list with a data frame with the attempts and a data frame with the improved mating plan
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#'
refinePlan = function(data,
                      plan,
                      K,
                      n.cross = 50,
                      max.cross = 4,
                      min.cross = 2,
                      culling.pairwise.k = .25,
                      max.attempts = 100){

  output = NULL
  moms.init = moms = plan$P1
  dads.init = dads = plan$P2
  allCrossesSub = allCrosses[which(!is.na(match(allCrosses$P1,moms))),]
  allCrossesSub = allCrossesSub[which(!is.na(match(allCrossesSub$P2,dads))),]

  initPlan = selectCrosses(data=allCrossesSub,
                           K=K,
                           n.cross=n.cross,
                           max.cross=max.cross,
                           min.cross=min.cross,
                           max.cross.to.search=nrow(allCrossesSub),
                           culling.pairwise.k = culling.pairwise.k)

  target.K = initPlan$summary$target.K
  target.Y = initPlan$summary$target.Y

  moms2try = setdiff(unique(allCrosses$P1),unique(moms))
  dads2try = setdiff(unique(allCrosses$P2),unique(dads))

  for(j in 1:max.attempts){
    tryPlan = NULL
    indicator = FALSE
    for(i in 1:length(moms2try)){
      allCrossesSub = allCrosses[which(!is.na(match(allCrosses$P1,c(moms,moms2try[i])))),]
      allCrossesSub = allCrossesSub[which(!is.na(match(allCrossesSub$P2,dads))),]
      tryCross = try(selectCrosses(data=allCrossesSub,
                                   K=K,
                                   n.cross=n.cross,
                                   max.cross=max.cross,
                                   min.cross=min.cross,
                                   max.cross.to.search=nrow(allCrossesSub),
                                   culling.pairwise.k = culling.pairwise.k),silent = TRUE)
      if(class(tryCross)=="list"){
        tryPlan = rbind(tryPlan,data.frame(tryCross[[1]],mom.tried=moms2try[i]))
      }else{
        tryPlan = rbind(tryPlan,data.frame(culling.pairwise.k=culling.pairwise.k,target.Y=NA,target.K=NA,mom.tried=moms2try[i]))
      }
    }

    check = which((tryPlan$target.Y > target.Y) & (tryPlan$target.K < target.K))

    if(length(check)>0){
      indicator = TRUE
      moms = c(moms,tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$mom.tried)
      target.K = tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$target.K
      target.Y = tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$target.Y
      output = rbind(output,data.frame(target.K=target.K,target.Y=target.Y,mom=tail(moms,1),dad=0))
    }else{ #try dad
      tryPlan = NULL
      for(i in 1:length(dads2try)){
        allCrossesSub = allCrosses[which(!is.na(match(allCrosses$P1,moms))),]
        allCrossesSub = allCrossesSub[which(!is.na(match(allCrossesSub$P2,dads2try))),]
        tryCross = try(selectCrosses(data=allCrossesSub,
                                     K=K,
                                     n.cross=n.cross,
                                     max.cross=max.cross,
                                     min.cross=min.cross,
                                     max.cross.to.search=nrow(allCrossesSub),
                                     culling.pairwise.k = culling.pairwise.k),silent = TRUE)
        if(class(tryCross)=="list"){
          tryPlan = rbind(tryPlan,data.frame(tryCross[[1]],dad.tried=dads2try[i]))
        }else{
          tryPlan = rbind(tryPlan,data.frame(culling.pairwise.k=culling.pairwise.k,target.Y=NA,target.K=NA,dad.tried=dads2try[i]))
        }
      }
      check = which((tryPlan$target.Y > target.Y) & (tryPlan$target.K < target.K))
      if(length(check)>0){
        indicator = TRUE
        dads = c(dads,tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$dad.tried)
        target.K = tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$target.K
        target.Y = tryPlan[check,][which.max(tryPlan[check,]$target.Y),]$target.Y
        output = rbind(output,data.frame(target.K=target.K,target.Y=target.Y,mom=0,dad=tail(dads,1)))
      }
    }
    if(!indicator){
      break
    }else{
      if(j==1){
        write.table(tail(output,1), col.names = TRUE, quote = FALSE)
      }else{
        write.table(tail(output,1), col.names = FALSE, quote = FALSE)
      }
    }
  }
  if(!is.null(output)){
    allCrossesSub = allCrosses[which(!is.na(match(allCrosses$P1,moms))),]
    allCrossesSub = allCrossesSub[which(!is.na(match(allCrossesSub$P2,dads))),]
    plan = selectCrosses(data=allCrossesSub,
                         K=Amat,
                         n.cross=n.cross,
                         max.cross=max.cross,
                         min.cross=min.cross,
                         max.cross.to.search=nrow(allCrossesSub),
                         culling.pairwise.k = culling.pairwise.k)
    return(list(steps=output,plan=tryCross))
  }
}
