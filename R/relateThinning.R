#########################################################################
#
# Package: BBMateAllocation
#
# File: relateThinning.R
# Contains: relateThinning
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Return a vector of individuals that maximizes criteria per cluster in the K matrix
#'
#' Return a vector of individuals that maximizes criteria per cluster in the K matrix
#'
#' @param K relationship matrix
#' @param criterion data frame with variable to maximize
#' @param threshold K threshold to cluster
#' @param max.per.cluster maximum number of entries per cluster to keep
#'
#' @return A vector with genotype names
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

relateThinning = function(K, criterion, threshold=0.5, max.per.cluster=2){
  diag(K) = 0
  K[K>threshold] = 1
  K[K<1] = 0
  size.init = nrow(K)

  for(i in 1:nrow(K)){
    index = colSums(K)
    index = sort(index,decreasing=TRUE)
    index = index[which(index>max.per.cluster)]
    if(length(index)==0)
      break

    tmp = criterion[match(names(K[which(K[names(index)[1],]==1),names(index)[1]]),
                          criterion[,1]),]
    tmp = tmp[order(tmp$Criterion,decreasing=TRUE),]
    tmp.rm = match(tmp[-c(1:max.per.cluster),1],rownames(K))
    K = K[-tmp.rm,-tmp.rm]
  }

  cat(paste0(nrow(K)," genotypes kept out of ", size.init,"\n"))
  return(rownames(K))
}
