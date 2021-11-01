# SimpleMating

R package to make simple mating allocation for small breeding programs

This package does culling of the possible parents based on relatedness and select the best mating based on a criterion given by the user. It follows the analysis pipeline:

# Data preparation
The user should have a list of genotypes they have available to mate. Each genotype should have a criterion (breeding/genetic value, selection/economic index) that should be maximized. Additionally, it should contain a relationship matrix between all the genotypes. Example of the input data is a data frame with two columns:

```{r, eval=FALSE}
print(criterion)
       Genotype Criterion
1 Genotype_0001  49.37977
2 Genotype_0002  50.82541
3 Genotype_0003  49.97405
4 Genotype_0004  50.84229
5 Genotype_0005  51.61311
6 Genotype_0006  49.23767
...
```

## Contact
Rodrigo R Amadeu  
rramadeu at gmail dot com  
Horticultural Sciences Department  
University of Florida  
https://rramadeu.github.io

## Installing and loading
Within R:
```
install.packages("devtools")
devtools::install_github("rramadeu/SimpleMating")

## Load 
library(SimpleMating)

## Data example
data("example_SimpleMating")

head(criterion)
Kmat[1:5,1:5]
## Keeping parents per cluster based on a K threshold
## Take the individuals with highest criteria to represent the cluster
keep = relateThinning(Kmat, criterion, threshold=0.5, max.per.cluster=4)

moms = criterion$Genotype[1:50] #first 100 genotypes
dads = criterion$Genotype #every genotype

## Creating data frame to do optimization
allCrosses = buildCrosses(moms, dads, keep, criterion, Kmat)

## Removing crosses already made in the past
## Let's say that you already did the following crosses in the past
pedigree = data.frame(P1=c("Genotype_0019","Genotype_0198"),
                      P2=c("Genotype_0043","Genotype_0019"))
head(pedigree)

## Removing them from the current plan
allCrosses = pastThinning(allCrosses,pedigree)

## Max gain, no restriction of max number of crosses
maxGainPlan = selectCrosses(data=allCrosses,
                            K=Kmat,
                            n.cross=100,
                            max.cross = 50,
                            min.cross = 1,
                            max.cross.to.search=100000,
                            culling.pairwise.k = 2)
maxGainPlan[[1]] #statistics of the plan
maxGainPlan[[2]] #mating plan

## Max gain not using duos with k>0.125 (1st cousins)
plan1 = selectCrosses(data=allCrosses,
                      K=Kmat,
                      n.cross=50,
                      max.cross=4,
                      min.cross=2,
                      max.cross.to.search=100000,
                      culling.pairwise.k = 0.125)
plan1[[1]] #statistics of the plan
plan1[[2]] #mating plan

## Max gain not using duos with k>0.0675 (2nd cousins)
plan2 = selectCrosses(data=allCrosses,
                      K=Kmat,
                      n.cross=50,
                      max.cross=4,
                      min.cross=2,
                      max.cross.to.search=100000,
                      culling.pairwise.k = 0.0675)
plan2[[1]] #statistics of the plan
plan2[[2]] #mating plan

```
