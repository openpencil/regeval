######0. This should be the git directory #####
workdir <- "~/regeval"
dir.create(workdir)
setwd(workdir)

###### I. libraries and locales ######
source('./regeval_packages.R')


##### II. Load in algorithms #####
source('./regeval_algorithms.R')

##### III. Get ingredients for simulating data #####

##### a. Load in data and pre-process ######
load(file="example_dataset.rda")
load(file="example_response.rda") 

##### b. Create the Design Matrix for the evaluation ######
createdesignmatrix <- function(countdata) {
  ## Remove rows and columns with zero rowSums and colSums ##
  ## Convert sequence counts to log probabilities ##
  logprobx <-  logprobclean(countdata)
  ## Address sum-to-1 redundancy ## 
  ## Select reference variable. (Equation 3) ##
  referencevarindex <- order(colSums(countdata), decreasing=T)[1]
  ## Drop reference variable from the design matrix. ##
  nonredundantmatrix <- logprobx[,-referencevarindex]
  ## Scale the design matrix to make mean 0 and variance 1 for all columns
  nonredundantxscaled <- apply(nonredundantmatrix, 2, scale)
  return(nonredundantxscaled)
}

designmatrix <- createdesignmatrix(exampledata)


##### c. Get signal to noise from original data #####
getsnr <- function(designmatrix, exampleresponse) {
  xscaled <- designmatrix
  y <- exampleresponse
  rownames(y) <- rownames(xscaled)
  ## scaling y to making mean 0 and variance 1 ##
  yscaled <- scale(y$response)
  
  ## This will generate a list of 3 items: ("glmnetmodel", "alpha", "lambda")
  cv_glmnet <- elasticnet_cv(x=xscaled,y=yscaled, standardize_cv=F, lambdaseq=NULL)
  
  ## Plugs the alpha from above to build an optimized model
  model_glmnet <- glmnet(x=xscaled, y=yscaled, alpha=cv_glmnet$alpha, standardize=F, family="gaussian")
  
  ## Predict at the lambda selected by the crossvalidation.
  predicted <- predict(model_glmnet, type="response", s=cv_glmnet$lambda, newx=xscaled, exact=F)
  
  ## Compute SNR
  # SNR is the l2 norm of the prediction divided by the l2 norm of the residual error term
  # SNR = \frac{||\bm{X}\beta||_{2}}{\sqrt{n} \sigma} 
  # where $||\bm{X}\beta||_{2}$ is the $l_{2}$ norm of the predicted response. 
  numerator <- sqrt(sum(predicted^2)) # L2 norm of (beta x X)
  denominator <- sqrt(sum((yscaled-predicted)^2)) # L2 norm of (y - beta X)
  snr <- numerator/denominator
  return(snr)
}

## Get signal to noise ratio for the original data ##
snrofdata <- getsnr(designmatrix, exampleresponse)


##### IV. Set up the models #####
# debug
# signaltonoise <- 4.6
# xscaled <- designmatrix
# cutoff <- 0.01
# signaltonoise <- snrofdata
# betavalue <- "unif"

##### Set up models and output for simulation #####
simulation <- function(xscaled, signaltonoise, cutoff, betavalue="unif"){  
  seed <- floor(runif(n=1, min=0, max=1000))
  p <- ncol(xscaled)
  n <- nrow(xscaled)
  
  ## cutoff for deciding the number of relevant variables ##
  relevant <- round(cutoff*p,0)
  if (relevant%%2 == 1){
    relevant <- relevant+1
  }
  ## One iteration of the permutation to generate betas across x and y ##
  ## sampling relevant index without replacement ##
  relindex <- sample(x=1:p, size=relevant, replace=F)
  print(sprintf("Chosen feature index:%s", relindex))
  
  ## Construct beta vector ##
  beta <- rep(x=0, length.out=p)
  
    ## Uniform distribution ##
  if (betavalue=="unif"){
    unifbetas1 <- runif(n=floor(relevant/2),min=-1, max=-0.5)
    unifbetas2 <- runif(n=floor(relevant/2),min=0.5, max=1)
    beta[relindex] <- c(unifbetas1, unifbetas2)
    colnames(xscaled) <- paste("v", 1:ncol(xscaled), sep="")
    colnames(xscaled) <- paste(colnames(xscaled), gsub("0\\.", "p", as.character(round(beta, 4))), sep="")
    colnames(xscaled) <- gsub("\\-p", "m", colnames(xscaled))
  
    ## Ones and minus ones ##
  } else if (betavalue=="ones"){
    ones <- relindex[1:round(relevant/2, 0)]
    minusones <- setdiff(relindex, ones)
    beta[ones] <- 1
    beta[minusones] <- -1
    colnames(xscaled) <- paste("v", 1:ncol(xscaled), sep="")
    colnames(xscaled)[ones] <- paste(colnames(xscaled)[ones], "one", sep="")
    colnames(xscaled)[minusones] <- paste(colnames(xscaled)[minusones], "mone", sep="") 
  }
  
  ## predicted y ##
  ## %*% is for matrix multiplication ##
  ## (Equation 2.) X.beta is being generated here. ##
  y_hat <- xscaled %*% beta
  
  ## Page 104, Buhlmann.##
  ## SNR is a ratio of the 2 SDs. 1st: l2norm_y / sqrt(n) * SD of error term from regression ##
  ## l2 norm of y_hat ##
  l2norm_y <- sqrt(sum(y_hat^2))
  
  ## Setting the error level or the noise.
  ## SNR = sqrt(power in the signal / power in the noise)
  ## (Equation 5.)
  ##  For a given SNR, we sampled the random error $\epsilon$ from a normal distribution $\mathcal{N}(0, \sigma^{2}) $ where:
  ## \sigma = \frac{||\bm{X}\beta||_{2}}{\sqrt{n} \cdot \text{SNR}}
  sigmanoise <- l2norm_y/(sqrt(n)*signaltonoise)
  ## rnorm(n,mean=0,sd=sigmanoise) is the error term ##
  y <- y_hat + rnorm(n,mean=0,sd=sigmanoise)
  y <- scale(y)
  out <- runmodels(xscaled, y, seed, iterations=10000, oracle=beta)
  return(out)
}
  

##### V. Set up input for running the simulation #####
setuplist <- list(setup=designmatrix)
betavalues <- c("ones", "unif")
cutoffrange <- c(0.02, 0.03, 0.04)
names(cutoffrange) <- paste("cutoff", cutoffrange*100, sep="")
snrout <- c(0.25, 16, 4.6)
simuloptions <- expand.grid(setup=names(setuplist), 
                            betavalue=betavalues, 
                            cutoff=cutoffrange, 
                            snr=snrout, 
                            simulnum=c(1:100))

### Test set-up for a few simulations
testsimul <- expand.grid(setup=names(setuplist), 
                            betavalue=betavalues, 
                            cutoff=cutoffrange[2:3], 
                            snr=snrout[2:3], 
                            simulnum=c(1:3))

### For a few number of simulations, activate this ###
# simuloptions <- testsimul

##### VI. Run simulation #####
runsimulation <- mclapply(1:nrow(simuloptions), function(rownum){
  simulrow <- simuloptions[rownum,]
  setup <- simulrow$setup
  betavalue <- simulrow$betavalue
  cutoff <- simulrow$cutoff
  snr <- simulrow$snr
  simulnum <- simulrow$simulnum
  simulout <- simulation(setuplist[[setup]], snr, cutoff, betavalue)
  if(length(simulout) > 0){
    saveRDS(simulout, sprintf("simul_%s_%0.2f_%0.2f_%s_%d.rds", setup, snr, cutoff, betavalue, simulnum))  
    }
  }, mc.cores=4)

# NB: Ignore these warnings.
# 20: In predict.lm.spike(object = train_path[[fold]][[feature]],  ... :
#                           Implicit intercept added to newdata.