##### Set number of cores in your machine here #####
numcore <- 2
# register cores for parallel processing
registerDoMC(cores=numcore)

#####I. Stability function#####
# mystability_glmnet: gives the matrix (or path matrix for stability selection)
# Equation 8
# For any given variable $j$, the probability that it belongs to $S(\lambda)$ 
# is estimated by the proportion of subsamples in which it was selected. 
# {eq:ssprob}
# P(j \in S(\lambda) )= \frac{1}{B} \sum_{b=1}^{B} 1_{j \in S_{\lambda}(b)}
# where $b$ is the index of any given subsample and $B$ is the total number 
# of subsamples. We refer to this probability as the variable inclusion probability. 
# The path matrix is a \lambda by j table containing probabilities of being selected.
# Note: Generously adapted from quadrupen!

##### Generate path matrix #####
#debug
# subsamples <- 104 (Has to be a multiple of the number of cores. Traditionally, 100, but this does not work if we have for instance, 8 cores.)
# alphaval <- best_cv_parameter
# randomize=T
# NOTE: this code works only for continuous outcomes!!! 
# x <- xscaled

mystability_glmnet <- function(x,
                               y,
                               lambdaseq,
                               subsamples  = 100,
                               alphaval=0,
                               mystandardize=T,
                               bolasso=F,
                               randomize=T,
                               setseed=T
                               ) {
  #Fixed variables
  #one way to improve the stability selection is to introduce randomness to the application of the l1 penalty for every parameter.
  # Weakness is the parameter that controls this variable specific randomness.
  weakness    = 0.5
  verbose     = TRUE
  p <- ncol(x)
  n <- nrow(x)
  if (bolasso == T) {
    sample.size <- n
    replace <- T
    } else {
    sample.size = floor(n/2)
    replace <- F
  }
  
  if (setseed == T){ 
    set.seed(101)
  }
  # Getting folds for stability selection
  folds       = replicate(subsamples, sample(1:nrow(x), sample.size, replace=replace), simplify=FALSE)
  mc.cores    = detectCores()
  penscale <- rep(1,p)
  nlambda1 <- 100
  glmnet.control(mnlam=nlambda1)
  
  ## Prepare blocks of sub samples to run jobs parallely
  blocs <- split(1:subsamples, 1:mc.cores)
  
  ## Efficiency trick for muliple cores
  # blocs are the subsamples and bloc.stability processes all the subsamples by splitting the subsamples across the cores.
  # one path matrix from each core
  # all path matrices will be averaged across the cores.
  #debug
  #subsets <- blocs[[1]]
  # s <- 2
  bloc.stability <- function(subsets) {
    # Matrix produces a sparse matrix
    select <- Matrix(0,nlambda1,p+1)
    subsamples.ok <- 0
    if (setseed == T){
      set.seed(101)
    }
    activebetas <- sapply(1:length(subsets), function(s){
      # print(s)
      if (randomize) { 
        current_penscale <- penscale / runif(p,weakness,1)
        # generating p random numbers for penalizing the betas
      } else { 
        current_penscale <- penscale 
      }
      x1=x[folds[[subsets[s]]], ]
      y1=y[folds[[subsets[s]]]]
      # family = "gaussian" for continuous outcomes
      model <- glmnet(x=x1, y=y1, family="gaussian",standardize=mystandardize, alpha=alphaval, 
                      nlambda=nlambda1, lambda=lambdaseq, penalty.factor=current_penscale)
      
      beta <- coef(model)
      # find which betas are non-zero
      # some corner cases lead to models failing. Possibly some singularity problem (determinants are zero)
      active <- apply(beta, 1, function(value) { value != 0 })
      return(active)
    }, simplify=F)
    
    
    # if glmnet stops before generating betas for 100 values of lambda, subsample.ok is not incremented and these short ones are tossed.
    # select is a 100 X p matrix. if a variable is active in every subsample, 
    # select will keep adding this variable across the active matrices.
    # this is aggregation across sub-samples.
    # length(subsamples.ok) = number of contrasts
    # length(subsets) --> number of splits of the 144 subsamples across the 8 cores = 144/8 = 18
    # subsamples.ok = vector of the number of splits that have a contrast  
    # print(subsamples.ok/length(subsets))    
    # dim(active) 100 X p: Each cell is boolean. Tells us whether the variable was selected at a particular value of lambda
    # if there are nlambda rows in the active matrix, aggregate the matrix else go to the next active matrix 
    # which has a complete set of nlambda1 rows
    for(s in 1:length(activebetas)){
      if (nrow(activebetas[[s]]) == nlambda1) {
        subsamples.ok <- subsamples.ok + 1
        select <- select + activebetas[[s]]
      }
    }

    if (subsamples.ok < 0.5*length(subsets)) {
      cat("\nWarning: more than 50% of the subsamples were discarded 
          in that core due to early stops of the fitting procedure. 
          You should consider largest 'min.ratio' or strongest 'lambda2'.")
    }
    # there are 8 blocs of 13 subsamples each
    # 
    #contribution of each core can be at most 1/number of cores.
    #subsamples.ok is the number of subsamples that are okay.
    #length(blocs) is number of subsets of n/2 subsamples.
    contribcore <- select/(subsamples.ok*length(blocs))
    return(contribcore)
  }
  
  ## Now launch all the jobs
  # 8 items in the list - 1 from each core. Each core gives us a path matrix
  prob.bloc <- mclapply(blocs, bloc.stability, mc.cores=mc.cores)
  
  ## Construct the probability path
  # path matrix = each row corresponds to value of lambda.
  # number in a cell = selection probability. % of subsamples where the variable was selected.
  # values in each path cell are frequency of any variable being selected in 104 (number of subsamples) trials for a given lambda.
  # this is adding the contributions across the 8 cores.
  path <- Matrix(0,nlambda1,p+1)
  for (b in 1:length(prob.bloc)) {
    path <- path + prob.bloc[[b]]
  }
  pathnew <- as.matrix(path)
  return(pathnew)
}

#####II.  Run elastic net double crossvalidation of alpha and lambda #####
# elasticnet_cv : generates optimal alpha and lambda
# The ENC procedure is a specific instance of the penalized regression setting.
# Equation 6
# {eq:glmnet}
# \hat{\beta} = \argmin_{\beta} \left\{ \sum_{i=1}^{n} |y_{i} -  X_{i}\beta|^{2} + \lambda(\alpha (||\beta||_{1} + (1-\alpha)\frac{1}{2}||\beta||^{2}_{2}) \right\}  
# where $\sum_{i=1}^{n} |y_i - X_i\beta|^2$ is the Mean Squared Error (MSE), 
# $\lambda$ is the tuning parameter that penalizes the MSE by the size of the regression coefficients.
# We are selecting the lambda and alpha by minimizing MSE 
# Reference in regeval_simulation.R: model_glmnet_cv <- elasticnet_cv(x=xscaled, y=y, standardize_cv=F, lambdaseq=lambdaseq_alph1)
#debug
# x <- xscaled
# y <- y
# standardize_cv <- F
# lambdaseq <- lambdaseq_alph1
# i <- 10

elasticnet_cv <- function(x, y, standardize_cv, lambdaseq) {
  # Set the number of folds to the number of samples
  numfolds <- dim(x)[1]
  # Grid for alpha crossvalidation
  alphas <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999)
  # In a k fold crossvalidation, fold ID gives the iteration in which the sample would be in the test set.  
  foldid <- sample(rep(seq(numfolds),length=dim(x)[1]))
  # Go through alpha grid
  # Run crossvalidation for lambda.
  # Each model for each alpha is run by a parallel core and put into a list called lassomodels
  lassomodels <- foreach(i = c(1:length(alphas))) %dopar% {
  # the function finds the best lambda for a given alpha
  # within each model there is cross-validation happening for lambda for each alpha.
  model <- try(cv.glmnet(x=x, y=y, family="gaussian",
                         nfolds=numfolds, 
                         type.measure="deviance", 
                         foldid=foldid,
                         standardize=standardize_cv, 
                         alpha=alphas[i],
                         parallel =T,
                         lambda=lambdaseq))
# For some of the folds, there will be <3 observations and for these option "grouped" will be set to FALSE
# Warning message:
# Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold 
  }
  
  # there are two lambdas per model
  # minimum lamda
  # lambda within 1 standard deviation of the error
  # find the best alpha
  best_alpha_index <- 0
  lowest_error <- 0
  for (i in c(1:length(alphas))) {
    # if a model fails, "try-error" will return true. "try-error" is an object that traps errors
    # inherits is a function that will be true if try-error has collected an error from the model
    # we want to avoid any errors recorded in "try-error" in the list of models we just generated
    # example too small a dataset
    if (!inherits(lassomodels[[i]], "try-error")) {
      # First we will find the index of the lambda corresponding to the lambda.min
      index <- which(lassomodels[[i]]$lambda.min == lassomodels[[i]]$lambda)
      # high lambda means more penalty.
      # lambdas are arranged from highest to lowest
      # alpha = 1 ==> lasso
      # alph = 0 ==> ridge    
      # cvm is the cross-validated error, in this case, deviance.
      # print(i)
      error <- lassomodels[[i]]$cvm[index]
      # print(error)
      if (best_alpha_index == 0 || error < lowest_error) {
        best_alpha_index <- i # picks an alpha from the grid of alphas
        lowest_error <- error # picks the lowest deviance from the grid
      }
    }
  }
  #print(best_alpha_index)
  out <- list(c(), c(), c())
  if (best_alpha_index != 0) {
    # print the lassomodel at the best_alpha_index
    lasso_model <- lassomodels[[best_alpha_index]]
    alphaval <- alphas[best_alpha_index]
    # Use lambda which gives the lowest cross validated error
    lambda <- lasso_model$lambda.min
    out <- list(lasso_model, alphaval, lambda)
  }
  names(out) <- c("glmnetmodel", "alpha", "lambda")
  return(out)
}


# Reference in regeval_simulation.R: model_stability <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F, 
#cv_parameter=updatedalpha, setseed=F, lambdaseq=lambdaseq)

#debug
# x <- as.matrix(xscaled)
# y <- as.matrix(y) 
# mystandardize <- F
# cv_parameter <- updatedalpha
# setseed <- F 
# lambdaseq <- lambdaseq
# bolasso <- F
# modality <- "stability"
# weightedlambda <- T

elastic_net_stability <- function(x, y, mystandardize, cv_parameter=-1,
                                       modality="stability", setseed=T, lambdaseq=NULL, 
                                       bolasso=F, weightedlambda=F) {
  # Removing zero columns silences the error on quadrupen's implementation of elastic net.
  zero_indicator_boolean_matrix <- (x==0)
  zero_columns <- which(colSums(zero_indicator_boolean_matrix) == nrow(x))
  if (length(zero_columns) >0) {
  x <- x[,-zero_columns]
  }
  if (cv_parameter == -1) {
      cv <- elasticnet_cv(x=x, y=y, seed=1001, family="gaussian", standardize_cv=mystandardize, lambdaseq=lambdaseq)
      best_cv_parameter <- cv$alpha
      best_l1_parameter <- cv$lambda
      } else {
      best_cv_parameter <- cv_parameter
      best_l1_parameter <- 0
    }
  print(sprintf("Best alpha:%0.2f",best_cv_parameter))
  
  if (setseed == T) {
    set.seed(101)
  }
  # mystability has an extra line: probabilities = path to expose the stability probabilities of ALL taxa.
  # probabilities is plotted using ggpplot in the source code of plot.stability
  # n is the number of samples
  # p is the number of x variables
  n <- dim(x)[1]
  p <- dim(x)[2]
  # floor of 5.3 is 5. Min of 5.3, 5 is 5
  num_variables <- min(floor(n/log(p)),p)
  
  #debug
  #alphaval <- best_cv_parameter
  # mystability uses either glmnet or quadrupen based on values of useglmnet
  path  <- mystability_glmnet(x,y, subsamples=104,
                              alphaval = best_cv_parameter,
                              mystandardize=mystandardize, setseed=setseed, 
                              lambdaseq=lambdaseq, bolasso=bolasso)
  
  
  if (modality == "stability") {
    # Incorporating Reviewer #2 suggestions
    if (weightedlambda == T) {
        weightedlambda <- rowMeans(path)
        weightedlambdanorm <- weightedlambda/sum(weightedlambda)
        # Rows are the models
        # Columns are the taxa
        # colSums - sum of the stability probabilities across a 100 models down each column. Each model is indexed by 1 lambda value.
        weightit <- sapply(1:ncol(path), function(pcol){
          weightedrow <- weightedlambdanorm*path[,pcol]
          return(weightedrow)
          ##### Quick plot for visualizing how weights are distributed across lambda grid######
          # ggframe <- data.frame(lambda=c(1:100), weights=weightedlambdanorm)
          # p <- ggplot(ggframe, aes(x=lambda, y=weights))
          # p <- p + geom_point()
          # p <- p + lightertheme
          # ggsave(filename = "weightedlambda.pdf", plot = p,width =4, height = 3, units = "in", limitsize = F )  
        })
        score <- colSums(weightit)
        
      } else {
        score <- colSums(path)
      }
    # columns of x are names of taxa.
    names(score) <- colnames(path)
    # Sort the stability scores
    score <- sort(score, decreasing=T)
    out <- list(best_l1_parameter, best_cv_parameter, score, x, y)
    names(out) <- c("lambda", "alpha", "stability", "xmatrix", "ymatrix")
  } else if (modality == "pfer") {
    select <- pfer_select(path)
    out <- list(best_l1_parameter, best_cv_parameter, select, x, y)
    names(out) <- c("lambda", "alpha", "pferselect", "xmatrix", "ymatrix")
  } else if (modality == "path") {
    # output path probs
    out <- list(best_l1_parameter, best_cv_parameter, path, x, y)
    names(out) <- c("lambda", "alpha", "path", "xmatrix", "ymatrix")
  }
  return(out)  
}


#######III. PFER algorithm of Meinshausen and Buhlmann ####### 
# pfer_select: gives names of selected variables.
# derived an upper bound on the expected number of false positives in any 
# given variable selection algorithm
# {eq:pfer}
# \E[FP] \le \frac{1}{2 \pi - 1 } \frac{q^{2}}{p}
# where $\E[FP]$ is the expected number of false positives, $q$ is the expected 
# number of selected influential variables and $\pi \in (\frac{1}{2}, 1)$ is a 
# tuning parameter. Based on Equation~\ref{eq:pfer}, they developed an algorithm 
# that computes the optimal set of variables which minimizes the number of false positives, 
# also known as per-family error rate (PFER). 
# Note: Generously adapted from quadrupen!
                  
pfer_select <- function(path, PFER=2, cutoff=0.6) {
# Try PFER = 1 vs 2 == E(FP) (these are the expected number of false positives)
# Choose cutoff in the range (0.6,0.9)
selection <- rep("unselected",ncol(path))
# estimate the average number of selected variables on the current path. 
# q = average number of selected variables at a given lambda (across 104 or however many subsamples)
# path matrix = each row corresponds to value of lambda.
# number in a cell = selection probability. % of subsamples where the variable was selected.
# values in each path cell are frequency of any variable being selected in 104 (number of subsamples) trials for a given lambda.
q <- rowSums(path >= cutoff)

# Pick the path controlling the PFER at the desired level
# qLim = maximum num of selected variables for a given PFER and a cutoff probability
# higher qLim means higher errors
# needs cutoff >= 0.5
p <- ncol(path)
qLim <- sqrt(PFER * (2 * cutoff-1) * p)
# Goal is to find the smallest index where q <= qLim 
# For each lambda, is the average number of selected variables <= maximum value ? (T/F)
# iq = [T T T T  F F F F F F F F] 
iq <- (q <= qLim) 
# Find the lowest index (value of lambda) at which the average # of vars selected exceeds max 
# which.min() = index of the first FALSE in the vector
# first FALSE is the lambda value where the q > qLim
# If this index = 1 (i.e. there are lots of selected vars even at the very beginning!), special cases:
# if iq is empty, set iq = 1
# if iq is not empty, set iq = 100 (max value of iq)
# which.min() - 1 = position of the last true variable 
# iq is going across lambda
# which.min()=1 means all are FALSE
#  => in this case, sum(iq) = # of lambda's where (q <= Lim)
# if that is zero
# it can mean:
## a) all entries are False i.e. q > qLim => then set iq to 1
## b) all entries are True i.e. q <= qLim => then seq iq to the last lambda index 100
iq <- ifelse(which.min(iq) != 1, which.min(iq)-1,
             ifelse(sum(iq) == 0, 1, length(iq)))
# At the above value of iq, make a decision on which variables have a selection probability > cutoff probability
#print(dim(path))
#print(iq)
#print(cutoff)
selection[path[iq, ] > cutoff] <- "selected"
# names of those variables
selected <- colnames(path)[which(selection == "selected")]
return(selected)
}

#####IV. BMA models #####
# expected_model_size <- bestmodelsize_bmac
# iterations <- 100
# y <- ydata
# x <- xdata
# defaultmode <- "logistic"
# iterations <- 10000
# expected_model_size <- 1
##### a. BMA models for regeval #####
runspike <- function(expected_model_size, y, x, iterations=10000, defaultmode="gaussian"){
  if (typeof(y) == "character"){
    y <- as.integer(as.factor(y))-1
  }
  if (defaultmode == "gaussian"){
    spikesmodel <- lm.spike(y ~ x, niter=iterations, expected.model.size=expected_model_size)
    spikesummary <- summary(spikesmodel, burn=round(0.1*iterations, 0))
    betas <- spikesummary$coefficients[,"mean"]
    inclprob <- spikesummary$coefficients[,"inc.prob"]
  } else if (defaultmode == "logistic"){
    x <- as.matrix(x)
    spikesmodel <- logit.spike(y ~ x, niter=iterations, expected.model.size=expected_model_size)
    spikesummary <- summary(spikesmodel, burn=round(0.1*iterations, 0))
    betas <- spikesummary[, "mean"]
    inclprob <- spikesummary[,"inc.prob"]
  }
  names(betas) <- gsub("xscaled|x", "", names(betas))
  #betindices <- grep("intercept", names(betas), invert=T, ignore.case=T)
  #betas <- betas[betindices]
  return(list(betas=betas, inclprob=inclprob, model=spikesmodel))
}

##### b. BMA models for bayesianmice #####
# defaultmode <- "logistic"
runbmac <- function(x, y, defaultmode="gaussian", modelsizearray=c(1, 3, 7, 10)){
  if (nrow(x) < 25) {
    folds<- 3 
  } else {
    folds<- 5
  }
  if (length(modelsizearray) > 0){
  bestmodelsize_bmac <- crossvalidate(x, y, folds=folds, variant="bmac",
                                      alphaval=-1, lambdaseq=c(), 
                                      bestlambdaindex=100, seed=101,  iterations=10000, defaultmode, modelsizearray)
  } else {
    bestmodelsize_bmac <-1
  }
  model_bmac <- runspike(expected_model_size=bestmodelsize_bmac, y=y, x=x, iterations=10000, defaultmode)
  return(list(ems=bestmodelsize_bmac, model=model_bmac))
}



#####V. Crossvalidate for the following models being evaluated ######
# Variable Selection based on Inclusion Probabilities at a single lambda (LSC, LRC)
# Bayesian model averaging with Spike-and-Slab regression (BMAC)
# debug
# folds <- 5
# variant <- "lrc"
# bestlambdaindex <- -1
# x <- xscaled
# y <- y
crossvalidate <- function(x, y, folds, variant, alphaval, lambdaseq, 
                          bestlambdaindex=-1, seed, iterations=10000, defaultmode="gaussian", modelsizearray=c(1, 3, 7, 10)) {
  set.seed(seed)
  
  ### Computing the folds
  foldid <- sample(rep(x=(1:folds),times=(nrow(x)/folds)), replace=F)
  xy <- cbind(y, x)
  
  ### Assigning test/train partitions
  assigntesttrain <- function(fold){
    testset <- xy[which(foldid==fold),]   
    trainset <- xy[which(foldid!=fold),] 
    out <- list(trainset, testset)
    names(out) <- c("train", "test")
    return(out)
  }
  testtrain <- lapply(1:folds, assigntesttrain)
  n <- nrow(x)
  p <- ncol(x)
  
  ## lrc is stability selection at a crossvalidated lambda with resampling
  ## this is technique behind the bolasso algorithm
  bolasso <- ifelse(grepl("^lrc", variant), T, F)
  
  # This is for models other than BMA, namely lrc, lsc
  ## Build the models
  if (grepl("bmac", variant)==F) {
    stability_path <- lapply(testtrain, function(elementlist){
      train_path <- mystability_glmnet(x=elementlist$train[,-1],
                                       y=elementlist$train[,1],
                                       lambdaseq=lambdaseq,
                                       subsamples = 104,
                                       alphaval=alphaval,
                                       mystandardize=F,
                                       bolasso=bolasso, #bootstrap rather than subsampling
                                       randomize=T,
                                       setseed=F)
    })
    #debug
    #fold_path <- stability_path[[2]] #lambdarow <- fold_path[5,] #bestlambdaindex <- 53 #proportion <- 0.01
    # Variable Selection based on Inclusion Probabilities at a single lambda (LRC, LSC)
    if (grepl("(lrc|lsc)", variant)) {
      gettopvar <- lapply(stability_path, function(fold_path){
        apply(fold_path, 1, function(lambdarow){
          # this is a n/log(p) based on data partition (subsampling or resampling)
          numcut <- floor(nrow(fold_path)/log(ncol(fold_path)))
          lambdasort <- sort(lambdarow, decreasing=T)
          topvar <- names(lambdasort[1:numcut])
        })
      })
    } 
    # datasubset <- testtrain[[1]] # gettopvars has five different sets of 100 lambda each # topvars <- gettopvar[[1]]
    ## Compute training set parameters
    ## For each of the 100 lambdas, there will be a MSE error.
    ## We have five models on five test tests and five errors
    ## Average these errors.
    ## Leave one out crossvalidation is unstable and non-deterministic with small sample-sizes. So 5-fold crossvalidation here.
    ## feature <- 70
    ## x <- 3
    maxrange <- 100
    trainreg <- sapply(1:maxrange, function(feature){
      sapply(1:folds, function(x){
        if(maxrange == 100){
          topvars <- gettopvar[[x]][,feature]
        }else{
          topvars <- gettopvar[[x]][[feature]]  
        }
        #construct regular regression with topvars for each lambda
        trainreg <- glmnet(x=testtrain[[x]]$train[, grep('ntercept', topvars, value=T, invert=T)], 
                           y=testtrain[[x]]$train[,1],
                           family="gaussian", alpha=alphaval, standardize=F, lambda=lambdaseq)
      },simplify=F)
    },simplify=F)
    #s=0 unregularized regression, lambda = 0. Select this model from trainreg glmnet 
    #exact = when there is no model corresponding to lambda=0 
    # inherentlambda <- 3 # x <- 4 #feature <- 1
    computemse <- sapply(1:maxrange, function(feature){
      mseval <- sapply(1:folds, function(fold){
        if (maxrange == 100){
          topvars <- gettopvar[[fold]][,feature]
        } else{
          topvars <- gettopvar[[fold]][[feature]]  
        }
        predicttestreg <- predict(object=trainreg[[feature]][[fold]], 
                                  newx=testtrain[[fold]]$test[,grep('ntercept', topvars, value=T, invert=T)], 
                                  s=0, type="response", exact=F) 
        coeff <- predict(object=trainreg[[feature]][[fold]], s=0, type="coefficients", exact=F) 
        
        y_test <- testtrain[[fold]]$test[,1]
        mse_test <- mean((y_test - predicttestreg)^2)
        return(mse_test)
      }, simplify=T)
      avgmse <- mean(mseval)
      return(avgmse)
    }, simplify=T)
  } else {
    # This is for BMA
    #debug #elementlist <- testtrain[[1]]
    ## Building the model over 4 different model sizes
    train_path <- mclapply(testtrain, function(elementlist){
      xscaled <- elementlist$train[,-1]
      y <-elementlist$train[,1]
      modelout <- mclapply(modelsizearray, function(num_vars){
        if (defaultmode == "gaussian") {
          spikesmodel <- lm.spike(y ~ xscaled, niter=iterations, expected.model.size=num_vars)
        } else if (defaultmode == "logistic"){
          xscaled <- as.matrix(xscaled)
          spikesmodel <- logit.spike(y ~ xscaled, niter=iterations, expected.model.size=num_vars)
        }
        return(spikesmodel)
      },mc.cores=3)
      return(modelout)
    }, mc.cores=4)
    
    #debug
    #feature <- 1
    #fold <- 1
    ## Taking the model parameters, we predict response on the test partition and compute mse. 
    # Average this MSE over the five folds
    # feature is the expected model size
    computemse <- sapply(1:length(modelsizearray), function(feature){
      mseval <- sapply(1:folds, function(fold){
        #print(sprintf("fold: %s", fold))
        #print(sprintf("feature:%s", feature))
        if (defaultmode == "gaussian"){
          #E(y|X)=\beta*X
          predicttestreg <- predict(object=train_path[[fold]][[feature]],
                                  newdata=as.matrix(testtrain[[fold]]$test[,-1]),
                                  type="response", burn=round(0.1*iterations, 0)) 
          medianpred <- apply(predicttestreg, 1, median)
          y_test <- testtrain[[fold]]$test[,1]
          mse_test <- mean((y_test - medianpred)^2)
        } else if (defaultmode == "logistic") {
          #P(y=1|X)
          #fold<-1
          #feature <- 6
          predicttestreg <- predict(object=train_path[[fold]][[feature]],
                                    newdata=as.matrix(testtrain[[fold]]$test[,-1]),
                                    type="prob", burn=round(0.1*iterations, 0)) 
          medianpred <- apply(predicttestreg, 1, function(x){
            class1 <- length(which(x>=0.50))
            class0 <- length(which(x<0.50))
            # For each sample this will return a 0 or 1
            voting <- ifelse(class1 >= class0, 1, 0)
            return(voting)
          })
          y_test <- testtrain[[fold]]$test[,1]
          # this is the error of classification of the categorical outcome
          mse_test <- length(which(y_test != medianpred))
        }
        return(mse_test)
      }, simplify=T)
      avgmse <- mean(mseval)
      return(avgmse)
    }, simplify=T)
  }
  #print(sprintf("MSE:%0.2f", computemse))
  
  ## Identify the parameters with the minimum MSE
  best_feature_index <- which.min(computemse)
  # for lambda, return  the  index
  best_feature <- ifelse(grepl("^bmac", variant), modelsizearray[best_feature_index], best_feature_index)
  return(best_feature)
}




########VI. Random utility functions ###########

#####a. Get lambda sequence #####
get.lambda1.l1 <- function(xty,nlambda1,min.ratio) {
  lmax <- max(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}

#####c. Normalize counts, smooth zeroes and convert to log proportions #####
normalize_smooth <- function(x) {
  x <- as.numeric(x)
  x <- x+1
  x <- log(x/sum(x))
}

logprobclean <- function(datainput){
  removezeros <- datainput[which(rowSums(datainput) != 0), which(colSums(datainput) != 0)]
  # add smoother. log zero is not defined, so add 1 and take the log. 1 is the smallest unit of reads.
  addsmoother <- removezeros+1
  #numrow <- 3
  probs <- t(apply(addsmoother, 1, function(x){
    y <- x/sum(x, na.rm=T)  
    return(y)
  }))
  logprobs <- log(probs)
  return(logprobs)
}

#####d. Calculate L1 and L2 norms #####
getl1norm <- function(somevalue){
  l1 <- sum(abs(somevalue))
  return(l1)
}
getl2norm <- function(somevalue){
  l2 <- sqrt(sum((somevalue)^2))
  return(l2)
}

#####e. Model components for all models ######

##### Borrowed from quadrupen!!! #####
standardize <- function(x,y,intercept,penscale,zero=.Machine$double.eps,
                        call.from.mv = FALSE) {
  
  n <- length(y)
  p <- ncol(x)
  ## ============================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  if (intercept) {
    xbar <- colMeans(x)
    ybar <- mean(y)
  } else {
    xbar <- rep(0,p)
    ybar <- 0
  }
  
  ## ============================================
  ## NORMALIZATION
  if (call.from.mv) { ## already scaled...
    normx <- rep(1,p)
  } else {
    normx <- sqrt(drop(colSums(x^2)- n*xbar^2))
    if (any(normx < zero)) {
      warning("A predictor has no signal: you should remove it.")
      normx[abs(normx) < zero] <- 1 ## dirty way to handle 0/0
    }
  }
  ## xbar is scaled to handle internaly the centering of X for
  ## sparsity purpose
  xbar <- xbar/normx
  normy <- sqrt(sum(y^2))
  
  ## normalizing the predictors...
  x <- sweep(x, 2L, normx, "/", check.margin = FALSE)
  ## and now normalize predictors according to penscale value
  if (any(penscale != 1)) {
    x <- sweep(x, 2L, penscale, "/", check.margin=FALSE)
    xbar <- xbar/penscale
  }
  ## Building the sparsely encoded design matrix
  if (inherits(x, "sparseMatrix")) {
    xs    <- as(x, "dgCMatrix")
    x     <- list(Xi = xs@i, Xj = xs@p, Xnp = diff(xs@p), Xx = xs@x)
    xty   <- drop(crossprod(y-ybar,scale(xs,xbar,FALSE)))
  } else {
    x     <- list(Xx = as.matrix(x))
    xty   <- drop(crossprod(y-ybar,scale(x$Xx,xbar,FALSE)))
  }
  
  return(list(xty=xty))
}

#####VII. Run all the models #######
#iterations <- 50
#oracle <- beta
runmodels <- function(xscaled, y, seed, iterations, oracle, weightedlambdamodels){
  # number of samples
  n <- nrow(xscaled)
  p <- ncol(xscaled)
  
  ###### TRUE BETA #####
  # "True" beta for evaluation gold-standard
  oraclebeta <- oracle
  names(oraclebeta) <- colnames(xscaled)
  truevars <- names(which(oraclebeta != 0))
  
  # standard lasso
  # get crossvalidated alpha and lambda
  # no scaling of alpha to start with when we generate the initial sequence of lambda
  
  # Generate the lambda sequence using quadrupen code
  # t(xscaled) %*% y or xscaled - mean(xscaled) * y - mean(y)
  # minimum value of lambda at which all betas are still zeroes
  betastillzero <- standardize(xscaled, y, intercept=T, penscale=rep(1,p))
  lambdaseq_alph1 <- get.lambda1.l1(xty=betastillzero$xty, nlambda1=100, min.ratio=0.01)
  
  ## all lambda values are on a log scale
  model_glmnet_cv <- elasticnet_cv(x=xscaled, y=y, standardize_cv=F, lambdaseq=lambdaseq_alph1)
  
  ## Updates after crossvalidation
  updatedalpha <- model_glmnet_cv$alpha
  lambdaindex_glmnet <- which(lambdaseq_alph1 == model_glmnet_cv$lambda)
  lambdaseq <- lambdaseq_alph1/updatedalpha
  
  ##################VARIANTS start here####################
  ###### ENC #####
  # Equation 6
  #eq:glmnet
  #\hat{\beta} = \argmin_{\beta} \left\{ \sum_{i=1}^{n} |y_{i} -  X_{i}\beta|^{2} + \lambda(\alpha (||\beta||_{1} + (1-\alpha)\frac{1}{2}||\beta||^{2}_{2}) \right\}  
  # where $\sum_{i=1}^{n} |y_i - X_i\beta|^2$ is the Mean Squared Error (MSE), 
  # $\lambda$ is the tuning parameter that penalizes the MSE by the size of the regression coefficients. 
  # $\alpha$ is a tuning parameter which balances the $l_{1}$ and $l_{2}$ penalties. 
  # An $\alpha$ of 1 promotes sparsity in the model while an $\alpha$ of 0 ensures that 
  # correlated variables are assigned similar regression coefficients. An optimal value of $\alpha$ finds a balance between the two penalties.
  # The following model has 100 values of lambda in a predefined sequence.
  model_enc <- glmnet(as.matrix(xscaled) , as.matrix(y), alpha=updatedalpha, 
                      standardize=F, family="gaussian", nlambda=100, lambda=lambdaseq)
  
  #extract coefficients from the model corresponding to cross-validated lambda.
  # exact=T re-estimates the coefficients at the specified value of lambda.
  enc_coefficients <- coef(object=model_enc, s=model_glmnet_cv$lambda, exact=F)
  betas_enc <- as.vector(enc_coefficients)
  # remove the beta of of the intercept term
  names(betas_enc) <- rownames(enc_coefficients)
  nonzerovars_enc <- names(betas_enc)[which(betas_enc!=0)]
  
  ###### SS #####
  # Variant: SS  | Subsampling for training sets
  # Equations 7, 8, 9
  # stability selection at the lasso CV param
  model_ss <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F, 
                                    cv_parameter=updatedalpha, setseed=F, lambdaseq=lambdaseq)
  model_ssw <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F, 
                                     cv_parameter=updatedalpha, setseed=F, 
                                     lambdaseq=lambdaseq, weightedlambda = T)
  ## controlling for nefarious conditions which generate NAs
  sumcheck <- function(modeldata) {
    sum_ss <- sum(modeldata$stability)
    booleancheck <-  (is.na(sum_ss) || is.nan(sum_ss) || is.infinite(sum_ss)) 
    return(booleancheck)
  }
  if (sumcheck(model_ss) == T){
    return(list())
  }
  if (sumcheck(model_ssw) == T){
    return(list())
  }
  ###### SR #####
  # Variant: SR | Resampling for training sets
  # Equations 7, 8, 9
  # stability selection at the lasso CV param
  # stability selection with resamples rather than subsamples
  model_sr <- elastic_net_stability(as.matrix(xscaled), as.matrix(y), mystandardize=F, 
                                    cv_parameter=updatedalpha, setseed=F, lambdaseq=lambdaseq, bolasso=T)
  model_srw <- elastic_net_stability(as.matrix(xscaled), as.matrix(y), mystandardize=F, 
                                    cv_parameter=updatedalpha, setseed=F, lambdaseq=lambdaseq, bolasso=T, weightedlambda = T)
  ## controlling for nefarious conditions which generate NAs
  if (sumcheck(model_sr) == T){
    return(list())
  }
  if (sumcheck(model_srw) == T){
    return(list())
  }

  ###### BMA #####
  # Equations 12, 13, 14, 15 
  # expected_model_size <- 1
  if (weightedlambdamodels == F){
    model_bma <- runspike(expected_model_size=1, y=y, x=xscaled, iterations)
    ###### BMAC #####
    # Equations 12, 13, 14, 15 
    # Tuning the expected.model.size. Expected model size is approximately the number of relevant variables we expect to be important to the response.
    # This is, at best, a guess. Default (and the most sparse option) is 1.
    # We tune the expected.model.size with crossvalidation.
    bestmodelsize_bmac <- crossvalidate(xscaled, y, folds=5, variant="bmac",
                                      alphaval=updatedalpha, lambdaseq=lambdaseq, 
                                      bestlambdaindex=lambdaindex_glmnet, seed=seed, iterations)
    model_bmac <- runspike(expected_model_size=bestmodelsize_bmac, y=y, x=xscaled, iterations)
    ###### PS #####
    # Variant: Stability selection minimizing per family error rate  (PS | PR)
    # {eq:pfer}
    # \E[FP] \le \frac{1}{2 \pi - 1 } \frac{q^{2}}{p}
    # where $\E[FP]$ is the expected number of false positives, 
    # $q$ is the expected number of selected influential variables and 
    # $\pi \in (\frac{1}{2}, 1)$ is a tuning parameter. 
    # Based on Equation~\ref{eq:pfer}, they developed an algorithm that computes the optimal set of variables which minimizes the number of false positives,  also known as per-family error rate (PFER).
    # Variant PS, based on subsamples
    model_ps <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F,
                                      cv_parameter=updatedalpha, modality="pfer", setseed=F, lambdaseq=lambdaseq)
    ###### PR #####
    # Variant: PR, based on resamples
    model_pr <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F,
                                    cv_parameter=updatedalpha, modality="pfer", setseed=F, lambdaseq=lambdaseq, bolasso=T)
    ###### LR #####
    # Variant: Variable LR Selection based on Inclusion Probabilities at a single lambda (LR) 
    # randomize: penalization lambda is decreased by random factors
    model_lr <- mystability_glmnet(xscaled, y, subsamples = 104, alphaval=updatedalpha, mystandardize=F,
                                 bolasso=T, randomize=T, setseed=F, lambdaseq=lambdaseq)
    # the lambda value for LR is determined by glmnet
    inclprob_lr <- model_lr[lambdaindex_glmnet,]
  
    # Ordered vars and probabilities according to the inclusion probabilities.
    orderedvars_lr <- names(inclprob_lr)[order(inclprob_lr, decreasing=T)]
    orderedprobs_lr <- inclprob_lr[orderedvars_lr]
    ###### LRC #####
    # Variant: Variable LR Selection based on Inclusion Probabilities at a single lambda (LRC) 
    bestlambdaindex_lrc <- crossvalidate(xscaled, y, folds=5, variant="lrc", 
                                       alphaval=updatedalpha, lambdaseq=lambdaseq, seed=seed)
    inclprob_lrc <- model_lr[bestlambdaindex_lrc,]
    # Ordered vars and probabilities according to the inclusion probabilities.
    orderedvars_lrc <- names(inclprob_lrc)[order(inclprob_lrc, decreasing=T)]
    orderedprobs_lrc <- inclprob_lrc[orderedvars_lrc]
    ###### LS #####
    # Variant: Variable LR Selection based on Inclusion Probabilities at a single lambda (LS) 
    model_ls <- elastic_net_stability(as.matrix(xscaled) , as.matrix(y), mystandardize=F,
                                      cv_parameter=updatedalpha, modality="path", setseed=F, lambdaseq=lambdaseq)
  
    # select the lambda from glmnet
    inclprob_ls <- model_ls$path[lambdaindex_glmnet,]
  
    # Ordered vars and probabilities according to the inclusion probabilities.
    orderedvars_ls <- names(inclprob_ls)[order(inclprob_ls, decreasing=T)]
    orderedprobs_ls <- inclprob_ls[orderedvars_ls]
  
    ###### LSC #####
    # Variant: Variable LR Selection based on Inclusion Probabilities at a single lambda (LSC) 
    bestlambdaindex_lsc <- crossvalidate(xscaled, y, folds=5, variant="lsc",
                                         alphaval=updatedalpha, lambdaseq=lambdaseq, seed=seed)
    inclprob_lsc <- model_ls$path[bestlambdaindex_lsc,]
    # Ordered vars and probabilities according to the inclusion probabilities.
    orderedvars_lsc <- names(inclprob_lsc)[order(inclprob_lsc, decreasing=T)]
    orderedprobs_lsc <- inclprob_lsc[orderedvars_lsc]
    ###### OUTPUT #####
    # The complete models
    models <- list(enc=model_enc,
                 ss=model_ss,
                 sr=model_sr,
                 ps=model_ps,
                 pr=model_pr,
                 lr=model_lr,
                 ls=model_ls,
                 bma=model_bma,
                 bmac=model_bmac)
    # Variables Selected: These methods do not give an inclusion probability list
    varsout <- list(truevars=truevars, #"Simulated Truth"
                    enc=nonzerovars_enc, #ENC
                    ps=model_ps$pferselect, #PS
                    pr=model_pr$pferselect #PR
                    )
    # Inclusion Probabilities
    problist <- list(lr=orderedprobs_lr, #LR
                     lrc=orderedprobs_lrc, #LRC
                     ls=orderedprobs_ls, #LS
                     lsc=orderedprobs_lsc, #LSC
                     ss=model_ss$stability, #SS
                     sr=model_sr$stability, #SR
                     bma=model_bma$inclprob, #BMA
                     bmac=model_bmac$inclprob #BMAC
                     )
  } else {
    models <- list(ss=model_ss,
                   ssw=model_ssw,
                   sr=model_sr,
                   srw=model_srw
                   )
    # Variables Selected: These methods do not give an inclusion probability list
    varsout <- list(truevars=truevars) #"Simulated Truth
    # Inclusion Probabilities
    problist <- list(ss=model_ss$stability,  #SS
                     ssw=model_ssw$stability,#SSW
                     sr=model_sr$stability,  #SR
                     srw=model_srw$stability #SRW
                     )
  }
  return(list(models=models, varsout=varsout, inclprobs=problist))
}

###### VIII. Non-Linear Models ######

## Random Forests Classifier ##
randomforestit <- function(x, y){
  set.seed(101)
  fy <- as.factor(y)
  #tunedforest <- tuneRF(x, fy, trace=T, plot=T, doBest=T, ntreeTry=5000)
  rfresults <- randomForest(x, fy, ntree=10000)
  impdf <- data.frame(importance(rfresults))
  impdf$var <- rownames(impdf)
  impsorted <- impdf[order(impdf$MeanDecreaseGini, decreasing=T),]
  #bootstrapval <- bootstrapit(x, y, iterations=100, numtree=1000)
  outlist <- list(medianval=impsorted, rfobject=rfresults)
  return(outlist)
}

