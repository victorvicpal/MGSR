## Multivariate Gaussian Subspatial Regression

This repository contains the required functions to apply the **MGSR**. 


----------
# Table of Contents
1. [Summary](#summary)
2. [Algorithm](#algorithm)
3. [Installation](#installation)
4. [Reference](#reference)
5. [Example](#example)


-------------
## Summary

Our proposal is a mixture between Factorial Techniques and Gaussian Processes.

Factorial techniques ([Biplot](https://www.wikiwand.com/en/Biplot), [PCA](https://www.wikiwand.com/en/Principal_component_analysis), [CA](https://www.wikiwand.com/en/Correspondence_analysis) ...) are useful to represent our observations in terms of unobserved variables called factors. 
All these techniques provide a set of coordinates linked to the observations, which display information of our analyzed variables. 
These kind of procedures are merely descriptive and have a low predictive power.

On the other hand, [Gaussian Processes](https://www.wikiwand.com/en/Gaussian_process) are statistical methods where observations occur in a continuous domain (mainly time or space). 
Furthermore, variables have a multivariate normal distribution. Gaussian Processes use similarity between points to predict the value in an unobserved point.

The **MGSR** applies a Gaussian Process Regression to the created subspatial of the Factorial Technique.
We make use of this subspace to simulate a continuous domain that permit the application of Gaussian Processes, such as Cokriging.

## Algorithm

1. Factorial Technique
2. [Cross-variogram](https://github.com/victorvicpal/MGSR/blob/master/crossvariogram.R)
    * [CV plot (cross variogram and LMC)](https://github.com/victorvicpal/MGSR/blob/master/plot.crossvariogram.R)
3. [Linear Model of Coregionalization](https://github.com/victorvicpal/MGSR/blob/master/lmc.R)
4. [Subspatial Grid](https://github.com/victorvicpal/MGSR/blob/master/grid.R)
5. [Cokriging](https://github.com/victorvicpal/MGSR/blob/master/cokrig.R)

> **Note:**

> - Cross-variograms usually follows a "Power Distribution" due to the small scale of Factorial Techniques.
> - Unlike geostatistical analyses, we don't have a real field where boundaries restrict our study. On the other hand, this aspect is more positive than negative because we can portray a more simple grid without losing information.

## Installation
```R
install.packages('devtools')
library(devtools)
install_github("victorvicpal/MGSR")
library(MGSR)
```

## Reference
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.264102.svg)](https://doi.org/10.5281/zenodo.264102)

## Example

###Iris Data
```R
data(iris)
```
####Data Exploration
```R
summary(iris)
apply(iris[which(iris$Species=='versicolor'),1:4],2,function(x,y) plot(density(x))) #density function
```
####"Versicolor" Train/subset
```R
Versicolor <- iris[which(iris$Species=='versicolor'),-5]

ind_subTrain <- sample(50,10)

subTrain <- Versicolor[ind_subTrain,]
Train <- Versicolor[-ind_subTrain,]
```
####Data Standardization
```R
means_train <- apply(Train,2,mean)
sd_train <- apply(Train,2,sd)
Train_st <- Train

for (i in 1:length(Train[1,]))
{Train_st[,i] <- (Train[,i]-means_train[i])/sd_train[i]}
```

####Principal Component Analysis
```R
PC_train <- princomp(Train_st)
biplot(PC_train)
```
####CrossVariogram
```R
CV_train <- crossvariogram(as.data.frame(PC_train$scores[,1:2]),as.data.frame(Train_st),10)
plot.crossvariogram(CV_train)
```

####lcm fitting
Tip: *Range value may vary. Check different values within the "Power" function.*
```R
RES_train <- lmc(CV_train,'Pow',1.7)
plot.crossvariogram(CV_train,RES_train)
```

####Grid
```R
xygrid <- GRID_MGSR(as.data.frame(PC_train$scores[,1:2]),0.05)
```

####Cokriging
```R
Z_train_st <- cokrig(RES_train,xygrid)
Z_train <- Z_train_st

for (i in 1:length(Train[1,]))
{Z_train[,i+2] <- Z_train_st[,i+2]*sd_train[i]+means_train[i]}
```

####Predicting Subset values
```R
ind_pred <- apply(dist2(subTrain[,1:4],Z_train[,3:6]),1,which.min)
residuales <- subTrain[,1:4]-Z_train[ind_pred,3:6]

par(mfrow=c(2,2))
apply(residuales,2,function(x) plot(density(x)))

for (i in 1:4)
{
  qqnorm(residuales[,i])
  qqline(residuales[,i])
}
```

Try to fit  a model with *virginica* and *setosa* species.
