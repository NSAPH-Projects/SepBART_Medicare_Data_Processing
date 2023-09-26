library("ModBart")
library("SoftBart") # MFM version
library("mvtnorm")
library("Rcpp")


## Call the SepBART function
source('00-functions-block.R')

## Load the data
try(load(file = "FILE PATHNAME of data_processed.RData"))
data_prcd <- as.matrix(data_prcd)

## Stack the data (with noise) to make the sample size larger.
data_stack <- data_prcd[rep(1:nrow(data_prcd), 2), ]
n_stack <- dim(data_stack)[1]
data_noise <- mvtnorm::rmvnorm(n_stack, sigma = diag(apply(data_stack,2,var))/4)
data_expd <- data_stack + data_noise
rm(data_stack, data_noise)

## Using the expanded data (data_expd) define the outcome, covariates, and exposure variables.
Y <- data_expd[,1] # death rate
X <- data_expd[,2:17] # pct_female"         "pct_dual"           "mean_bmi"          
# "smoke_rate"         "hispanic"           "pct_blk"            "medhouseholdincome"
# "medianhousevalue"   "poverty"            "education"          "popdensity"        
# "pct_owner_occ
W <- data_expd[,18:25]

## Define two levels of exposures we would like to compare.
W0 <- apply(W, 2, quantile, 0.4)
W1 <- apply(W, 2, quantile, 0.6)

## Fit the model.
fit <- SepBART_blk(Y = Y, X = X, W = W,
        Xtest = NULL, W1 = W1, W0 = W0,
        nMCMC=1000,BurnIn_portion=0.2,stepsize=5)
