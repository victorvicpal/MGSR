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
```
install.packages('devtools')
library(devtools)
install_github("victorvicpal/MGSR")
library(MGSR)
```

## Reference
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.264102.svg)](https://doi.org/10.5281/zenodo.264102)

## Example

###Slump Data
[Information](https://archive.ics.uci.edu/ml/datasets/Concrete+Slump+Test)
```
install.packages('RCurl')	#HTTP requests package
library('RCurl')
url1 <- getURL('https://archive.ics.uci.edu/ml/machine-learning-databases/concrete/slump/slump_test.data')
Slump_data <- read.csv(text = url1,colClasses=c("NULL",NA,"NULL","NULL",NA,"NULL",NA,NA,"NULL",NA,NA),col.names = c('no','cement','slag','fly_ash','water','sp','coarse','fine','slump','flow','Comp_str'))
```

We create Train/Test datasets
```
ind_test <- sample(length(Slump_data[,1]),20)
Test <- Slump_data[ind_test,]
Train <- Slump_data[-ind_test,]
```

####Descriptive analysis
[Slump test info](en.wikipedia.org/wiki/Concrete_slump_test)
```
summary(Train)
#density functions
par(mfrow=c(2,5))
apply(Train,2,function(x) plot(density(x)))
apply(Train,2,boxplot)

#Data Standardization
means_train <- apply(Train,2,mean)
sd_train <- apply(Train,2,sd)
Train_st <- Train

for (i in 1:length(Train[1,]))
{Train_st[,i] <- (Train[,i]-means_train[i])/sd_train[i]}
```

####Principal Component Analysis
```
PC_train <- princomp(Train_st)
biplot(PC_train)
```
####CrossVariogram
```
n_opt <- n_pairs_opt(as.data.frame(PC_train$scores[,1:2]),as.data.frame(Train_st),9,35)
plot(n_opt$dif_pairs)
CV_train <- crossvariogram(as.data.frame(PC_train$scores[,1:2]),as.data.frame(Train_st),16)
plot.crossvariogram(CV_train)
```

####lcm fitting
Range value may vary. Check different values within the "Power" function.
```
RES_train <- lmc(CV_train,'Pow',1.4)
plot.crossvariogram(CV_train,RES_train)
```

####Grid
```
xygrid <- GRID_MGSR(as.data.frame(PC_train$scores[,1:2]),0.05)
```

####Cokriging
```
Z_train_st <- cokrig(RES_train,xygrid)
Z_train <- Z_train_st

for (i in 1:length(Train[1,]))
{Z_train[,i+2] <- Z_train_st[,i+2]*sd_train[i]+means_train[i]}
```

####Predicting Test values
```
ind_pred <- apply(dist2(Test[,1:4],Z_train[,3:6]),1,which.min)
residuales <- Test[,5:6]-Z_train[ind_pred,7:8]
par(mfrow=c(1,2))
plot(density(Test$flow),col='blue',main='Flow')
lines(density(Z_train$flow[ind_pred]),col='red')
plot(density(Test$Comp_str),col='blue',main='Comp_Str')
lines(density(Z_train$Comp_str[ind_pred]),col='red')
par(mfrow=c(2,1))
qqnorm(residuales$flow)
qqline(residuales$flow)
qqnorm(residuales$Comp_str)
qqline(residuales$Comp_str)
```
