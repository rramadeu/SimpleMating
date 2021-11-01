#########################################################################
#
# Package: BBMateAllocation
#
# File: pastThinning.R
# Contains: pastThinning
#
# Written by Rodrigo Amadeu
#
# First version: Oct-2021
# Last update: 25-Oct-2021
#
#########################################################################

#' Given a data frame and a pedigree, filter out crosses in the data frame that are already present in the pedigree
#'
#' Given a data frame and a pedigree, filter out crosses in the data frame that are already present in the pedigree
#'
#' @param data data frame with possible crosses to filter out
#' @param pedigree pedigree to not repeat
#'
#' @return Filtered data frame
#'
#' @examples
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

pastThinning = function(data, ped){
  ped1 = paste0(ped[,1],"_crossed_",ped[,2])
  ped1 = unique(c(ped1,paste0(ped[,2],"_crossed_",ped[,1])))

  pednew = paste0(data[,1],"_crossed_",data[,2])
  torm = which(!is.na(match(pednew,ped1)))
  if(length(torm)>0)
    data = data[-torm,]
  cat(paste0(length(torm)," lines were removed\n"))
  rownames(data)=1:nrow(data)
  return(data)
}
