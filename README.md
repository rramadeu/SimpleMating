# SimpleMating

R package to make simple mating allocation for small breeding programs. It works fine for a few thousands of genotypes (<5K).

This package does culling of the possible parents based on relatedness and select the best mating based on a criterion given by the user. It follows the analysis pipeline:

# 0) Data preparation
The user should have a list of genotypes they have available to mate. Each genotype should have a criterion (breeding/genetic value, selection/economic index) that should be maximized. Additionally, it should contain a relationship matrix between all the genotypes. Example of the input data is a data frame with two columns:

```{r, eval=FALSE}
> print(criterion[1:5,])
       Genotype Criterion
1 Genotype_0001  49.37977
2 Genotype_0002  50.82541
3 Genotype_0003  49.97405
4 Genotype_0004  50.84229
5 Genotype_0005  51.61311
...
```

And a matrix with the same dimension having the genotypes in order:
```{r, eval=FALSE}
> print(Kmat[1:5,1:5])
              Genotype_0001 Genotype_0002 Genotype_0003 Genotype_0004 Genotype_0005
Genotype_0001     1.0298019     0.5688419     0.5688419     0.5688419     0.5688419
Genotype_0002     0.5688419     1.0298019     0.5688419     0.5688419     0.5688419
Genotype_0003     0.5688419     0.5688419     1.0298019     0.5688419     0.5688419
Genotype_0004     0.5688419     0.5688419     0.5688419     1.0298019     0.5688419
Genotype_0005     0.5688419     0.5688419     0.5688419     0.5688419     1.0298019
```

# 1) Thinning by relatedness
This step selects the top `n` individuals within a family clustered based on a relatednesses `threshold`. For example, if `n=2` and `threshold=0.50` the algorithm will identify all the group of individuals that shared has more than 0.50 in the relationship matrix and select the top 2 individuals with highest criterion to be used as possible parents in the mating plan

# 2) Building crosses
Here, the use needs to define possible moms, dads, and genotypes to keep based on step 1. Then, the algorithm will return a data frame with all the possible crosses, the average of the criteria of each parent duo, and the relationship between the parents for each mating allocation.

# 3) Thinning by pedigree
In this step the user can provide a historial pedigree of the breeding program, and then the algorithm will remove from the actual list crosses already made in the past.

# 4) Select crosses
Core algorithm, the user provides the number of crosses, the maximum number of crosses per parent, the minimum number of crosses per parent (expected), and a maximum value of relationship between parents that can be used to create a cross. The algorithm sort by criterion (maximum to minimum) and build a mating plan adding lines from top to bottom until the number of crosses is reached. During the process, it doesn't add crosses out of the range of the criteria. If a parent already reached its maximum, it is not add anymore in the plan.

## Contact
Rodrigo R Amadeu  
rramadeu at gmail dot com  
Horticultural Sciences Department  
University of Florida  
https://rramadeu.github.io

## Installing
Within R:
```
install.packages("devtools")
devtools::install_github("rramadeu/SimpleMating")
```

## Usage
```
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

# Comments and future directions
For small datasets, this simple solution appeared reasonable when contrasted with AlphaMate solution. However, SimpleMating solution gives the user more flexibility and control during the process. I tried different evolutionary algorithms to improve the mating design after selectCrosses function, however they just resulted in marginal gains, I let such attempts in the package for future reference (`mutatePlan`, `evoluteCandidates`, and `devolutePlan` functions). From my perspective, a good filtering and reduction of the possible candidates based on relatedness and criteria is a simple solution for small breeding populations.

I wouldn't recommend the usage of the genomic relationship matrix for this strategy in the case of the Blueberry Breeding Program (BBP). The numerator relationship matrix has clear values that the user can use to filter (for example, filter by first cousing, by full-sib etc). On the other hand, the GRM doesn't have clear-cut values and can mislead some of the decisions. For example, in the BBP the GRM value for fullsib families can have average from .15 to .60, if it is in one or in the other end it depends on how close or far the family is from the "average" relatedness in the population. As the BBP uses within and between familiy selection during the breeding process, it end-up having several full-sib families as possible parents to avoid such inbreeding, a clear-cut is necessary.

