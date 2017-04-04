##Ayesha Rasool



## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----Q1, message=FALSE, warning=FALSE------------------------------------
## Load example data, from `Hiiragi2013`.
if (FALSE) {  # install packages if not available
source("https://bioconductor.org/biocLite.R")
biocLite("Hiiragi2013")
install.packages("mlr")
}

install.packages("knitr")
library("knitr")
library("Biobase")
library("Hiiragi2013")
library("mlr")   
library("ggplot2")
library("matrixStats")

## ------------------------------------------------------------------------
data( "x", package = "Hiiragi2013" )
table( x$sampleGroup )

## ------------------------------------------------------------------------
rowV <- data.frame( v = rowVars(exprs(x)) )

## ----nonspecific_filtering-----------------------------------------------
selectionThreshold <- 10^(0.5) 
selectedFeatures  <- ( rowV$v > selectionThreshold )
print(sum(selectedFeatures))
embryoSingleCells <- data.frame( t(exprs(x)[selectedFeatures, ]), check.names = TRUE )

## ------------------------------------------------------------------------
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

## ------------------------------------------------------------------------
embryoSingleCells$tg = sample(embryoSingleCells$tg) # randomise

## ------------------------------------------------------------------------
with( embryoSingleCells, table( tg ) )   # show counts in each class

## ------------------------------------------------------------------------
task <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

## ------------------------------------------------------------------------
fv = generateFilterValuesData(task, method = "anova.test")
plotFilterValues(fv)  # look at the distribution of the gene p-values
# the values on y-axis are between 0-15.

## ------------------------------------------------------------------------
#select just the top 30 genes based on significance of fold change.
filtered.task = filterFeatures(task, method = "anova.test", abs=30) 

## ------------------------------------------------------------------------
#build a simple LDA classifier.
lrn = makeLearner("classif.lda")

## ---- results='hide',message=FALSE---------------------------------------
rdesc <- makeResampleDesc( method = "CV", stratify = TRUE, iters = 2) 
r <- resample(learner = lrn, task = filtered.task, resampling = rdesc)
r$aggr   # mean error rate
#mmce.test.mean = 0.228
#predictive accuracy = 0.772

## ---- final exam --------------------------------------------------------

# Is your result surprising. Why? Explain this result and the cause, and a solution.
# The result is surprising because I expected the mean error rate to be negligible or much closer to zero.
# This also affects the predictive accuracy, which at 0.772  is farther away from 1 (desired accuracy).
# Ambroise et al, (2002), have explained similar results in their paper where they were expecting that:
# it is possible to construct a prediction rule from only a few genes such that it has a negligible prediction error rate.
# Same strategy was applied here by filtered task where top 30 genes were selected based on their fold change.
# The authors recommended to use 10-fold CV which we did here as well.
# But they also  suggest using the so-called .632+ bootstrap error estimate designed to handle overfitted prediction rules with small number of genes.
# They also proved that when correction is made for the selection bias, the cross-validated error is no longer zero for a subset of only a few genes.
# Therefore, I would recommend in our case, not to use such a small subset of genes with filtered task.
# Also, doing the filtering by p-value globally causes a bias, so in order to avoid that bias, we cannot this
# error rate for the whole dataset and it would be useful to apply unfiltered task and see if the 
# predictive accuracy improves.
# Moreover it would be fruitful to not randomize the data this time, or compare the rates for shuffled 
# and unshuffled data. If the randomized error rate is not close to the unshuffled one, then we know that 
# we have introduced a selection bias in the data, that needs to be removed. 


#Now write the R code to do a cross-validation using a single combined classifier 
#that combines both the gene selection filter and the LDA prediction model into a 
#single procedure. 

#Then do a 10-fold cross-validation, so this combined model will be applied separately 
#to each fold, and so the filtering is done separately on each training set. 
# (hint refer to the help for the makeFilterWrapper() function). 


# with shuffling

# Note: The labels have been shuffled so the data is essentially random, 
# and so we should not be able to predict the class labels from expression

lrn = makeFilterWrapper(lrn, fw.method = "anova.test", fw.abs=30)
rdesc <- makeResampleDesc( method = "CV", stratify = TRUE, iters = 2) 
r <- resample(learner = lrn, task = task, resampling = rdesc)
r$aggr   # mean error rate
#mmce.test.mean = 0.475

#predictive accuracy :(accuracy = 1-error rate)
#predictive accuracy =  0.525

# What error rate do you get now and why?
# The error rate we got is twice as much as the original one. In lrn, we are performing the 
# anova.test for 30 genes with highest fold change, and since the data is shuffled, our result is 
# biased towards those selective genes / probes. In order to reduce this rate, it would be 
# better to perform CV on unshuffled data and get an unbiased rate.


# without shuffling

## ------------------------------------------------------------------------
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

## ------------------------------------------------------------------------
with( embryoSingleCells, table( tg ) )   # show counts in each class

## ------------------------------------------------------------------------
task <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

## ------------------------------------------------------------------------
fv = generateFilterValuesData(task, method = "anova.test")
plotFilterValues(fv)  # look at the distribution of the gene p-values
#y-axis values range between 0-150.

## ------------------------------------------------------------------------
filtered.task = filterFeatures(task, method = "anova.test", abs=30) 

rdesc <- makeResampleDesc( method = "CV", stratify = TRUE, iters = 2) 
r <- resample(learner = lrn, task = task, resampling = rdesc)
r$aggr   # mean error rate went down after unshuffling.
#mmce.test.mean = 0.1192

#predictive accuracy :(accuracy = 1-error rate)
#predictive accuracy = 0.8808. The accuracy has improved as compared to the one with shuffled data.
# Unshuffling the data and using the non-filtered task has removed the bias from the data.

##*************************************************************************************************##

## SVM - Question 2 - without shuffling
#re-run the script above lines 85 - 99 (until filtered task is defined).
#found the following commands from:
#https://mlr-org.github.io/mlr-tutorial/devel/html/tune/index.html


# Note: sigma was defined to be (2^(-2:2)) and also (2^(-12:12)), since the latter improved the rates
# and is said to be tuned more accurately, it will be used, throughout the whole question # 2.
ps = makeParamSet(
  makeDiscreteParam("C", values = 2^(-12:12)),
  makeDiscreteParam("sigma", values = 2^(-12:12))
)

ctrl = makeTuneControlGrid()
rdesc <- makeResampleDesc( method = "CV", stratify = TRUE, iters = 3) 
res = tuneParams("classif.ksvm", task = filtered.task, resampling = rdesc, par.set = ps, control = ctrl)

outer = makeResampleDesc("CV", iters = 3)
r = resample(lrn, filtered.task, resampling = outer)

# outer loop: [Resample] Result: mmce.test.mean=0.0695
#predictive accuracy with SVM = 0.9305, for unshuffled data, in comparison, LDA produced 88% of accuracy. 

## SVM - Question 2 - with shuffling
## ------------------------------------------------------------------------
embryoSingleCells$tg <- factor( ifelse( x$Embryonic.day == "E3.25", "E3.25", "other") )

## ------------------------------------------------------------------------
embryoSingleCells$tg = sample(embryoSingleCells$tg) # randomise

## ------------------------------------------------------------------------
with( embryoSingleCells, table( tg ) )   # show counts in each class

## ------------------------------------------------------------------------
task <- makeClassifTask( id = "Hiiragi", data = embryoSingleCells, target = "tg" )

## ------------------------------------------------------------------------
fv = generateFilterValuesData(task, method = "anova.test")
plotFilterValues(fv)  # look at the distribution of the gene p-values

## ------------------------------------------------------------------------
filtered.task = filterFeatures(task, method = "anova.test", abs=30) 

ps = makeParamSet(
  makeDiscreteParam("C", values = 2^(-12:12)),
  makeDiscreteParam("sigma", values = 2^(-12:12))
)

ctrl = makeTuneControlGrid()
rdesc <- makeResampleDesc( method = "CV", stratify = TRUE, iters = 3) 
res = tuneParams("classif.ksvm", task = filtered.task, resampling = rdesc, par.set = ps, control = ctrl)
outerS = makeResampleDesc("CV", iters = 3)
rS = resample(lrn, filtered.task, resampling = outerS)
#[Tune] Result: C=4; sigma=0.25 : mmce.test.mean=0.148
#[Resample] Result: mmce.test.mean=0.257

#predictive accuracy for SVM = 0.743 with shuffled data.
#this accuracy is 8% less accurate than with unshuffled data.



