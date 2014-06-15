workdir <- "~/simulclean/"
dir.create(workdir)
setwd(workdir)

###### I. libraries and locales ######
source('~/microbiome/regeval_packages.R')

###### extra libraries ######
## Follow instructions for installing the BoomSpikeSlab package from:
## https://sites.google.com/site/stevethebayesian/googlepageforstevenlscott/boom
# 1. Install gfortran from here: http://cran.r-project.org/bin/macosx/tools/
# gfortran-4.2.3.pkg
# 2. install.packages("BH")
# From https://sites.google.com/site/stevethebayesian/googlepageforstevenlscott/boom
# 3. Download Boom_0.1.tar.gz (or later versions) and BoomSpikeSlab_0.4.1.tar.gz (or later versions) and copy into ~/R_source
# 4. install.packages(c("~/R_source/Boom_0.1.tar.gz"), repos=NULL, type="source")
# 5. install.packages(c("~/R_source/BoomSpikeSlab_0.4.1.tar.gz"), repos=NULL, type="source")
# load library
library("BoomSpikeSlab")

##### II. Load in algorithms #####
source('~/microbiome/regeval_algorithms.R')

##### III. Get ingredients for simulating data #####

##### a. Load in data and pre-process ######
load(file="example_dataset.rda")
load(file="example_response.rda") 

# convert into log probabilties, remove columns with all zeros.
#debug
#datainput <- exampledata
logprobx <-  logprobclean(exampledata)

##### b. Visualize the covariance structure of data ######
displaycovariancestruc <- function(logprobx){
  par(pin=c(2,2))
  par(plt=c(2,2,2,2))
  cov_matrix <- cov(logprobx) 
  ggmelt <- melt(as.matrix(cov_matrix))
  colnames(ggmelt) <- c("rows", "columns", "correlation")
  c <- ggplot(ggmelt, aes(x=rows, y=columns))
  c <- c + geom_tile(aes(fill = correlation), colour = "white")
  #c <- c + scale_fill_gradient(low="white", high = "#9E0041")
  myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
  c <- c + scale_fill_gradientn(colours = rev(myPalette(100)))
  c <- c + theme_grey(base_size = 9) 
  c <- c + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  c <- c + theme(legend.position = "none", axis.ticks = element_blank())
  c <- c + theme(axis.text.x = element_blank())
  ggsave("covarmatrix_entiredata.pdf", plot=c, width=10, height=9, units="in", limitsize=F)
}
## Generate the covariance visualization to see how complex high dimensions can be!
displaycovariancestruc(logprobx)

##### c. Get signal to noise from original data #####
getsnr <- function(logprobx, exampleresponse) {
  x <- logprobx
  y <- exampleresponse
  rownames(y) <- rownames(x)
  
  # scaling to making mean 0 and variance 1
  xscaled <- apply(x, 2, scale)
  yscaled <- scale(y$response)
  
  ## This will generate a list of 3 items: ("glmnetmodel", "alpha", "lambda")
  cv_glmnet <- elasticnet_cv(x=xscaled,y=yscaled, standardize_cv=F, lambdaseq=NULL)
  
  ## Plugs the alpha from above to build an optimized model
  model_glmnet <- glmnet(x=xscaled, y=yscaled, alpha=cv_glmnet$alpha, standardize=F, family="gaussian")
  
  ## Predict at the lambda selected by the crossvalidation.
  predicted <- predict(model_glmnet, type="response", s=cv_glmnet$lambda, newx=x, exact=F)
  
  ## Compute SNR
  #SNR is the l2 norm of the prediction divided by the l2 norm of the residual error term
  # SNR = \frac{||\bm{X}\beta||_{2}}{\sqrt{n} \sigma} 
  # where $||\bm{X}\beta||_{2}$ is the $l_{2}$ norm of the predicted response. 
  numerator <- sqrt(sum(predicted^2)) # L2 norm of (beta x X)
  denominator <- sqrt(sum((yscaled-predicted)^2)) # L2 norm of (y - beta X)
  snr <- numerator/denominator
  return(snr)
}

snrofdata <- getsnr(logprobx, exampleresponse)

##### IV. Set up the models #####
# debug
# signaltonoise <- 4.6
# xscaled <- logprobx
# cutoff <- 0.01
# signaltonoise <- snrofdata
# betavalue <- "unif"

##### Set up models and output for simulation #####
simulation <- function(xscaled, signaltonoise, cutoff, betavalue="unif"){  
  seed <- floor(runif(n=1, min=0, max=1000))
  p <- ncol(xscaled)
  n <- nrow(xscaled)
  
  ## cutoff for deciding the number of relevant variables
  relevant <- round(cutoff*p,0)
  if (relevant%%2 == 1){
    relevant <- relevant+1
  }
  #### One iteration of the permutation to generate betas across x and y ###
  # sampling relevant index without replacement
  relindex <- sample(x=1:p, size=relevant, replace=F)
  print(sprintf("Chosen feature index:%s", relindex))
  
  # Construct beta vector
  beta <- rep(x=0, length.out=p)
  
    # Uniform distribution
  if (betavalue=="unif"){
    unifbetas1 <- runif(n=floor(relevant/2),min=-1, max=-0.5)
    unifbetas2 <- runif(n=floor(relevant/2),min=0.5, max=1)
    beta[relindex] <- c(unifbetas1, unifbetas2)
    colnames(xscaled) <- paste("v", 1:ncol(xscaled), sep="")
    colnames(xscaled) <- paste(colnames(xscaled), gsub("0\\.", "p", as.character(round(beta, 4))), sep="")
    colnames(xscaled) <- gsub("\\-p", "m", colnames(xscaled))
  
    # Ones and minus ones 
  } else if (betavalue=="ones"){
    ones <- relindex[1:round(relevant/2, 0)]
    minusones <- setdiff(relindex, ones)
    beta[ones] <- 1
    beta[minusones] <- -1
    colnames(xscaled) <- paste("v", 1:ncol(xscaled), sep="")
    colnames(xscaled)[ones] <- paste(colnames(xscaled)[ones], "one", sep="")
    colnames(xscaled)[minusones] <- paste(colnames(xscaled)[minusones], "mone", sep="") 
  }
  
  ### predicted y ####
  ## %*% is for matrix multiplication
  # Equation 2. X.beta is being generated here.
  y_hat <- xscaled %*% beta
  
  #### Page 104, Buhlmann.
  #### SNR is a ratio of the 2 SDs. 1st: l2norm_y / sqrt(n) * SD of error term from regression
  #### l2 norm of y_hat ####
  l2norm_y <- sqrt(sum(y_hat^2))
  
  # Setting the error level or the noise.
  # SNR = sqrt(power in the signal / power in the noise)
  # Equation 5.
  #  For a given SNR, we sampled the random error $\epsilon$ from a normal distribution $\mathcal{N}(0, \sigma^{2}) $ where:
  # \sigma = \frac{||\bm{X}\beta||_{2}}{\sqrt{n} \cdot \text{SNR}}
  sigmanoise <- l2norm_y/(sqrt(n)*signaltonoise)
  # rnorm(n,mean=0,sd=sigmanoise) is the error term
  y <- y_hat + rnorm(n,mean=0,sd=sigmanoise)
  y <- scale(y)
  out <- runmodels(xscaled, y, seed, iterations=10000, oracle=beta)
  return(out)
}
  

##### V. Set up input for running the simulation #####
setuplist <- list(setup=logprobx)
betavalues <- c("ones", "unif")
cutoffrange <- c(0.01, 0.02, 0.03, 0.04)
names(cutoffrange) <- paste("cutoff", cutoffrange*100, sep="")
snrout <- c(0.25, 16, 4.6)
simuloptions <- expand.grid(setup=names(setuplist), 
                            betavalue=betavalues, 
                            cutoff=cutoffrange, 
                            snr=snrout, 
                            simulnum=c(1:100))

### Test set-up for a few simulations
testsimul <- expand.grid(setup=names(setuplist), 
                            betavalue=betavalues[2], 
                            cutoff=cutoffrange[2:3], 
                            snr=snrout[2:3], 
                            simulnum=c(1:5))

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