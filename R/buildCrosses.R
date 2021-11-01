#########################################################################
#
# Package: BBMateAllocation
#
# File: buildCrosses.R
# Contains: buildCrosses, meltK, bv
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Create a data frame with all possible crosses
#'
#' Create a data frame based on given moms, dads, maximum relatedness, and maximum of individuals per relatedness cluster, return all crosses possible.
#'
#' @param moms vector of possible moms to use
#' @param dads vector of possible dads to use
#' @param keep vector of genotypes to keep
#' @param criterion data frame with variable to maximize
#' @param K relationship matrix between all genotypes
#'
#' @return A data frame with all possible crosses.
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

buildCrosses = function(moms,dads,keep,criterion,K){

  moms = moms[which(moms %in% keep)]
  dads = dads[which(dads %in% keep)]
  bvCriterion = bv(criterion)
  KCriterion = meltK(K)
  fullPedigree = expand.grid(P1=moms,P2=dads)
  fullPedigree = fullPedigree[-which(as.character(fullPedigree$P1)==as.character(fullPedigree$P2)),]

  bvPed = paste0(bvCriterion$P1,"_x_",bvCriterion$P2)
  KPed = paste0(KCriterion$P1,"_x_",KCriterion$P2)
  cross2keep = paste0(fullPedigree$P1,"_x_",fullPedigree$P2)
  bvCriterion = bvCriterion[match(cross2keep,bvPed),]
  KCriterion = KCriterion[match(cross2keep,KPed),]
  output = data.frame(bvCriterion,K=KCriterion$K)
  output = output[order(output$Y,decreasing=TRUE),]
  if(any(which(duplicated(cbind(output$P1,output$P2)))))
    output = output[-which(duplicated(cbind(output$P1,output$P2))),] #removing duplicated
  rownames(output) = 1:nrow(output)
  cat(paste0(nrow(output)," possible crosses were assigned\n"))
  return(output)
}

bv = function(data){
  addeffect = as.vector(data[,2])
  X=kronecker(t(addeffect),matrix(1,length(addeffect),1))
  X=(X+t(X))/2
  #   X[upper.tri(X)] <- NA
  diag(X) = NA #no inbreeding
  X=cbind(which(!is.na(X),arr.ind = TRUE),na.omit(as.vector(X)))
  X=as.data.frame(X)
  X[,1]=data[X[,1],1]
  X[,2]=data[X[,2],1]
  colnames(X)=c("P1","P2","Y")
  rownames(X)=NULL
  return(X)
}


meltK = function(X){
  namesK = rownames(X)
  #    X[upper.tri(X)] = NA
  diag(X) = NA #no inbreeding
  X=cbind(which(!is.na(X),arr.ind = TRUE),na.omit(as.vector(X)))
  X=as.data.frame(X)
  X[,1]=namesK[X[,1]]
  X[,2]=namesK[X[,2]]
  colnames(X)=c("P1","P2","K")
  rownames(X)=NULL
  return(X)
}
