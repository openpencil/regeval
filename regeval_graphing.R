######0. This should be the git directory #####
workdir <- "~/regeval"
dir.create(workdir)
setwd(workdir)

###### Ia. libraries and locales ######
source('./regeval_packages.R')


###### Ib. Data that was used for the simulation ######
load(file="example_dataset.rda")
load(file="example_response.rda") 
numsamples <- nrow(exampledata)

###### II. Load common functions ######
evaluateselection <- function(selectedvars, truevars, allvars){
  # remove intercept term: not necessary
  allvars <- setdiff(allvars, "(Intercept)")
  
  nonselectedvars <- setdiff(allvars, selectedvars)
  falsevars <- setdiff(allvars, truevars)
  #false positive rate
  falsepositives <- length(setdiff(selectedvars, truevars))
  truenegatives <- length(intersect(falsevars, nonselectedvars))
  fpr <- falsepositives/(falsepositives+truenegatives)
  #true positive rate
  truepositives <- length(intersect(selectedvars, truevars))
  tpr <- truepositives/(length(truevars))
  # recall
  recall <- tpr
  # precision
  precision <- truepositives/length(selectedvars)
  if(!is.nan(precision) && !is.nan(recall) && precision == 0 && recall == 0){
    fscore <- 0
  } else {
    fscore <- 2*(precision*recall)/(precision+recall)  
  }
  eval <- c(fpr, tpr, fscore, precision, recall)
  names(eval) <- c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL")
  return(eval)
}

# Calculate the plotting points for the ROC curve
genthreshold <- function(multicol){  
# Uncomment for debugging!
#   multicol <- bindqx[model== "lr" & 
#                        snr == "16.00" & 
#                        cutoff == 0.02 & 
#                        simulnum == "1" &
#                        onesunif == "ones"]
  variable <- multicol[,colnames(multicol)=="variables",with=F] 
  qx <- multicol[,grep("^q", colnames(multicol)),with=F]
  rescale<- multicol[,which(colnames(multicol)=="rescale"),with=F] 
  lrvalue <- log(rescale$rescale + exp(qx[,1,with=F]))
  # For each quantile cutoff, find variables whose proportions fall above the cutoff
  # cut <- qx[1,15,with=F]
  metrics <- sapply(qx[1,], function(cut){
    madethecut <- variable$variable[which(lrvalue$q1 >= cut)]
    truevars <- grep("one|m|p", variable$variables, value=T)
    evalresult <- evaluateselection(truevars, selectedvars=madethecut, allvars=variable$variables)
  })
  metricsrow <- unlist(c(metrics["FPR",], metrics["TPR",], qx[1,]))
  names(metricsrow) <- c(paste("FPR", names(metricsrow[1:length(qx)]), sep="_"), 
                         paste("TPR", names(metricsrow[1:length(qx)]), sep="_"),
                         paste("tau", names(metricsrow[1:length(qx)]), sep="_"))
  
  calcauc <- sapply(2:ncol(metrics), function(colindex){
    # area of trapezium = 1/2 * base * sum of || sides
    # x=fpr, y=tpr
    auc <- 0.5 * abs((metrics["FPR",colindex] - metrics["FPR",colindex-1]) *
                       (metrics["TPR",colindex] + metrics["TPR",colindex-1]))
  })
  aucvalue <- sum(calcauc)
  return(data.table(t(metricsrow), auc=aucvalue))
}

# How do the inclusion probabilities correlate with the relative magnitudes of true betas?
extractbeta <- function(varname){
  #v212m3759 
  betaval <- ifelse(grepl("(p|m)", varname), as.numeric(gsub(".*[mp](\\d{1,}$)", "\\1", varname))*1e-4, 0)
  return(betaval)
}

rankcorrelate <- function(multicol){
  betavector <- extractbeta(multicol$variables)
  scor <- cor(multicol$rescale, betavector, method="spearman")
  scor <- ifelse(is.na(scor), 0, scor)
  return(list(scor=scor))
}

# Method of lagged differences(LD), scree ##
selectionvector <- function(multicol){
  #   # #Uncomment for debugging!
  #     multicol <- bindpre[model== "lr" & 
  #                          snr == "16.00" & 
  #                          cutoff == 0.02 & 
  #                          simulnum == "1" &
  #                          onesunif == "ones",
  #                          c("variables", "rescale"), with=F]
  multiorder <- multicol[order(multicol$rescale, decreasing=T),]
  # Find the first order lagged differences
  difforder <- -(diff(multiorder$rescale, 1))
  # Find the index corresponding to the maximum lagged difference
  maxdiffindex <- which(max(difforder) == difforder)[1]
  # Build a vector of 1s and 0s. 1s: Selected variables, 0s. Variables not selected.
  select <- ifelse(1:nrow(multiorder) <= maxdiffindex, 1, 0)
  names(select) <- multiorder$variable
  select <- select[multicol$variable]
  selectedvars <- multicol$variable[which(select==1)]
  truevars <- grep("one|m|p", multicol$variables, value=T)
  metrics <- evaluateselection(selectedvars, truevars, allvars=multicol$variables)
  return(data.frame(FPR=metrics[1], TPR=metrics[2], FSCORE=metrics[3], PRECISION=metrics[4], RECALL=metrics[5]))
}

#Sample the rows from many columns
mysamplemulti <- function(multivalue, minval, do=T){
  set.seed(203)
  if (do==T){
    sampleit <- sample(1:nrow(multivalue), minval, replace=F)
    final <- multivalue[sampleit,]
  } else {
    final <- value
  }
  return(final)
}

###### III. Read in results of simulations ######
allrds <- list.files(path=workdir, pattern=".*\\d{1,}\\.rds", full.names=T)

###### IV. Harvest all the simulation results ######
# debug
# rdsfile <- allrds[1]
# model <- "lr"
simulationresults <- lapply(allrds, function(rdsfile){
  print(sprintf("loading:%s",rdsfile))
  snr <- gsub(".*(16\\.00|25|4\\.60)_0.*rds$", "\\1", rdsfile)
  cutoff <- gsub(".*_(0.0\\d)_(\\d*).*.rds", "\\1", rdsfile)
  simulnum <- gsub(".*_(\\d*).*.rds", "\\1", rdsfile)
  onesunif <- gsub(".*(ones|unif).*rds$", "\\1", rdsfile)
  rdslist <- readRDS(rdsfile)
  rdsnames <- setdiff(names(rdslist$betasout), "truebeta")
  ## Gather up inclusion probabilities for all variables
  # model <- "bma"
  modelprob <- sapply(names(rdslist$inclprobs), function(model){
  # catch corner cases where simulations didn't produce any output
    if (length(rdslist$inclprobs[[model]])==0) {
      return(list())
    } else {
      modelip <- data.table(ldply(rdslist$inclprobs[[model]]))
      setnames(modelip, colnames(modelip), c("variables", "ipval"))
      modelip[,variables:=gsub("^xscaled", "", variables),]
      modelip[,model:=model,]
      return(modelip)
    }
  },simplify=F)
  
  # Extract selected vars for enc, ps and pr (methods with no inclusion probabilities)
  # Evaluate model performance by comparing selected vars against the true vars
  #db 
  #model <- "enc"
  truevars <- rdslist$varsout$truevars
  allvars <- rownames(rdslist$models$enc$beta)
  modeleval <- sapply(setdiff(names(rdslist$varsout), "truevars"), function(model){
    selectedvars <- rdslist$varsout[[model]]
    evaluation <- evaluateselection(selectedvars, truevars, allvars)
    return(evaluation)
  }, simplify=F)
  evalply <- data.table(ldply(modeleval))
  #db
  #modelp <- "bma"
  ipvalues <- rbindlist(modelprob)
  ipvals <- ipvalues[,`:=`(snr=snr, cutoff=cutoff, simulnum=simulnum, onesunif=onesunif), ]
  evals <- evalply[,`:=`(snr=snr, cutoff=cutoff, simulnum=simulnum, onesunif=onesunif), ]
  return(list(eval=evalply, ip=ipvals)) 
})

saveRDS(simulationresults, "simulationresults.rds")

###### V. Combine all the simulation results ######
# simulationresults <- readRDS("simulationresults.rds")
# type <- "ip"
combineandsample <- sapply(c("ip", "eval"), function(type){
  gatherup <- sapply(1:length(simulationresults), function(num){
    collectip <- simulationresults[[num]][[type]]
    },simplify=F)
  bindup <- rbindlist(gatherup)
  if (type == "eval"){
    setnames(bindup, colnames(bindup)[1], "model")
    lengthsubset <- bindup[,list(lenval=length(FPR)),by=c("model", "snr", "cutoff", "onesunif")]
  } else {
    # Remove the intercept
    bindup <- bindup[variables!="(Intercept)",]
    # For the stability models, IPs range from 0 to 100, because IPs are summed over all lambdas.
    # Divide these stability values by 100
    bindpre <- bindup[,rescale:=ifelse(grepl("ss|sr", model), ipval/100, ipval),]
    nonzero <- setdiff(bindpre$rescale, 0)
    grandmin <- min(nonzero)
    # Sequence of thresholds equally spaced on the logscale
    # 22 inclusion probability thresholds between 0.00 and 1.01 on the log-linear scale
    # Adding 0.01. 
    qx <- c(seq(log(grandmin), 0, length.out=21), 0.01)
    names(qx) <- paste("q", 1:length(qx), sep="")  
    qxmat <- matrix(rep(qx, each=nrow(bindpre)), ncol=length(qx))
    colnames(qxmat) <- names(qx)
    bindqx <- cbind(bindpre, qxmat)
    colinterest <- c("variables", grep("^q", colnames(bindqx), value=T), "rescale")
    rocdata <- bindqx[ ,genthreshold(.SD), by=c("model", "snr", "cutoff", "simulnum", "onesunif"), .SDcols=colinterest]  
    # Save this as an RDS just in case the simulation is massive and you need the results in a jiffy.
    saveRDS(rocdata, "rocdata.rds")
    # Reads the partial results 
    # rocdata <- readRDS("rocdata.rds")
    lengthsubset <- rocdata[,list(lenval=length(auc)),by=c("model", "snr", "cutoff", "onesunif")]
    corcol <- c("variables", "rescale")
    bindunif <- bindpre[onesunif=="unif"]
    bindcor <- bindunif[,list(scor=rankcorrelate(.SD)$scor), by=c("model", "snr", "cutoff", "simulnum"), .SDcols=corcol]
    # Save partial results 
    saveRDS(bindcor, "bindcor.RDS")
    # Reads the partial results
    # bindcor <- readRDS("bindcor.RDS")
    # bindpre consists of all simulations; unifs and ones
    selectld <- bindpre[,selectionvector(.SD), by=c("model", "snr", "cutoff", "simulnum", "onesunif"), .SDcols=corcol]
    # Rename rocdata as bindup because the eval branch has bindup
    bindup <- rocdata
  }
  # determine minimum number of simulations available for a set of conditions.
  minnumsample <- min(lengthsubset$lenval)
  sdcols <- setdiff(colnames(bindup), c("simulnum","model", "snr", "cutoff", "onesunif"))
  
  # From each set of conditions (model, snr, cutoff, onesunif), 
  # sample the minimum of "minnumsample" simulations without replacement 
  # This is to make sure that every set of conditions has the same number of simulations.
  sampledrows <- bindup[,mysamplemulti(.SD, minval=minnumsample, do=T), 
                        by=c("model", "snr", "cutoff", "onesunif"), .SDcols=sdcols]
  if (type == "eval"){
    return(sampledrows)
  } else{
    sampledcorr <- bindcor[,mysamplemulti(.SD, minval=minnumsample, do=T), by=c("model", "snr", "cutoff"), .SDcols=c("scor")]
    sampledld <- selectld[,mysamplemulti(.SD, minval=minnumsample, do=T), 
                                 by=c("model", "snr", "cutoff", "onesunif"), 
                                 .SDcols=c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL")]
    return(list(ip=sampledrows, cor=sampledcorr, selectld=sampledld))
  }
},simplify=F)

# Ignore this warning.
# Warning messages:
# In ifelse(grepl("(p|m)", varname), as.numeric(gsub(".*[mp](\\d{1,}$)",  ... :
# NAs introduced by coercion

# Save partial results
saveRDS(combineandsample, "combinesample.rds")

# Read in partial results in case you are picking up analysis after a break.
# combineandsample  <- readRDS("./combinesample.rds")


###### VI. Get ROC data corresponding to the median AUC ######
# The ROC curve for a binary classification problem plots 
# the true positive rate as a function of the false positive rate. 
# The points of the curve are obtained by sweeping the 
# inclusion probability threshold from 100% to 0%. 

inclprob <- combineandsample$ip$ip
# calculate mean tpr, fpr, auc and tau for each set of conditions
tprfprtauauc <- grep("tpr|fpr|tau|auc", colnames(inclprob), ignore.case=T, value=T)
meantprfpr <- inclprob[,lapply(.SD, mean), by=list(model, snr, cutoff, onesunif), .SDcols=tprfprtauauc]
tprfprtau <- grep("TPR|FPR|tau", colnames(inclprob), value=T)

meltrocdata <- function(datasub){
  essentialcols <- setdiff(colnames(datasub), tprfprtau)
  meltitbytau <- sapply(1:(length(tprfprtau)/3), function(num){
    tausuffix <- rep(sprintf("q%s", num), nrow(datasub))
    tpfptau <- c(sprintf("TPR_q%s", num), sprintf("FPR_q%s", num), sprintf("tau_q%s", num), essentialcols)
    gettripletcols  <- cbind(datasub[,tpfptau,with=F], tausuffix)
    setnames(gettripletcols, colnames(gettripletcols), c("TPR", "FPR","tau", essentialcols, "tauindex"))
    return(gettripletcols)
  },simplify=F)
}
rocmean <- rbindlist(meltrocdata(meantprfpr))

#########VII. General graphing code #############
## colours: colours for each variant ##
newcolours <- c("#d73027" , "#f46d43", "#e08214", "#66bd63", "#1a9850", "#1d91c0", "#225ea8","#c2a5cf", "#9970ab","#dd3497", "#ae017e")
names(newcolours) <- c("BMA", "BMAC", "ENC", "LR", "LRC", "LS", "LSC", "PR", "PS", "SR", "SS")

modifylabel <- function(var, value){
  value <- as.character(value)
  if (var=="cutoff") {
    value <- paste("P", value, sep=":")
  }
  else if (var=="snrname") {
    value[value=="SNR:04.6"] <- "SNR:4.6"
  }
  return(value)
} 

plotlayers <-function(plottingdata, ymetric, ggp){
  p <- ggp + facet_grid(snrname~cutoff,  labeller=modifylabel)
  p <- p + xlab("Approaches") + ylab(ymetric)
  p <- p + lightertheme
  return(p)
}

rocfin <- datasub
prepdata <- function(rocfin, oneorunif){
  if (oneorunif != "") {
    rocfin <- rocfin[onesunif==oneorunif,]
  }
  getsnrname <- function(snrval){
    newsnr <- ifelse(grepl("16", snrval), "SNR:16.0",
                     ifelse(grepl("25", snrval), "SNR:0.25",
                            ifelse(grepl("4.6", snrval), "SNR:04.6", "undefined")))
    return(newsnr)
  }
  rocfin[,snrname:=getsnrname(snr),]
  rocfin[,model:=toupper(model),]
  rocfin[,modelnum:=getmodelnum(rocfin),]
  return(rocfin)
}

getmodelnum <- function(plottingdata){
  index <- sprintf("%02d", match(plottingdata$model, sort(names(newcolours))))
  modelnumber <- paste(index, plottingdata$model, sep=":")
  return(modelnumber)
}

stripnum <- function(value){
  value <- gsub("\\d{2}:", "", value)
  return(value)
}

##########VIII. Plot ROC curves and AUC boxplots ##############

# what <- "mean"
# oneunif <- "unif"

## Make extra columns for labels for graphs ##
findata <- sapply(list(rocmean, inclprob), function(datasub){
  sapply(unique(datasub$onesunif), function(type){
    fin <- prepdata(datasub, type)
  },simplify=F)
},simplify=F)
names(findata) <- c("mean", "all")

## Plot ROC curves ##
plotroc <- function(what, oneunif){
  rocfin <- findata[[what]][[oneunif]]
  ggp <- ggplot(rocfin, aes(x=FPR, y=TPR, colour=model))
  ggp <- ggp + geom_abline(aes(intercept=0, slope=1), linetype="dashed")
  ggp <- ggp + geom_path()
  p <- plotlayers(plottingdata=rocfin, ymetric="", ggp=ggp)
  p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
  p <- p + labs(colour="Approaches")
  p <- p + scale_color_manual(values=newcolours)
  p <- p + scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  ggsave(sprintf("ROC_%s_%s.pdf", what, oneunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
}

# FINAL PLOTS
# Receiver Operating Characeristics (ROC curves): Figure 2a and Figure 2b 
plotroc(what="mean", oneunif="ones")
plotroc(what="mean", oneunif="unif")


#db
# oneunif <- "unif"
## Plot AUC boxplots ##
plotaucbox <- function(oneunif){
  rocfin <- findata$all[[oneunif]]
  ggp <- ggplot(rocfin, aes(x=modelnum, y=auc, fill=model))
  ggp <- ggp + geom_boxplot(outlier.size=1, fatten=0.5)
  p <- plotlayers(plottingdata=rocfin, ymetric="AUC", ggp=ggp)
  p <- p + labs(fill="Approaches")
  p <- p + scale_fill_manual(values=newcolours)
  p <- p + scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + theme(axis.text.x = element_text(size=10, angle = 290, hjust = 0, colour="black"))
  p <- p + scale_x_discrete(label=stripnum)
  ggsave(sprintf("AUC_boxplot_%s.pdf", oneunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
}
# FINAL PLOTS
# Variability in AUC: Figure 3a and Figure 3b
plotaucbox("unif")
plotaucbox("ones")

###### IX. Plot FSCORE across the approaches ######

## Get data for plotting from combineandsample ##
allbeta <-combineandsample$eval
allldselect <- combineandsample$ip$selectld
vardata <- list(noip=allbeta, yesip=allldselect)

vardata_ones <- sapply(vardata, function(datasub){
  pdata <- prepdata(rocfin=datasub, oneorunif="ones")
  return(pdata)
},simplify=F)

vardata_unifs <- sapply(vardata, function(datasub){
  pdata <- prepdata(rocfin=datasub, oneorunif="unif")
  return(pdata)
},simplify=F)

# debug
# metric <- "FSCORE"

## Vector for iterating over beta=+-1 and beta \in \mathcal{U}
ou <- c("ones", "unifs")

# debug
# onesorunif <- ou[2]
sapply(ou, function(onesorunif){
  x <- "yesip" # methods that give IP rankss
  y <- "noip"  # methods that do not give IP ranks
  print(sprintf("vardata_%s", onesorunif))
  datasub <- get(sprintf("vardata_%s", onesorunif)) 
  
  mergeddata <- rbind(datasub[[x]], datasub[[y]])
  
  # Revised calculation for Fscore
  mergeddata[,fscore_rev:=ifelse(PRECISION == 0 | RECALL == 0, 0, FSCORE),]
  
  # Rename the revised Fscore back
  mergeddata[,FSCORE:=fscore_rev,]
  
  # Melt data 
  mergemelt <- melt(mergeddata, measure.vars=c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL"))
  
  # Eliminate 0.25 SNR
  mergemelt <- mergemelt[snr!=25]
  
  # FINAL PLOTS
  # Variability in FSCORE: Figure 4a and Figure 4b
  # Variability in FPR: Figure 5a and Figure 5b
  sapply(c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL"), function(metric){
    allmodels <- mergemelt[variable==metric & !is.na(value)]
    p <- ggplot(allmodels, aes(x=modelnum, y=value, fill=model))
    p <- p + geom_boxplot(outlier.size=1, fatten=0.5)
    p <- p + facet_grid(snrname~cutoff, labeller=modifylabel)
    p <- p + scale_fill_manual(values=newcolours)
    p <- p + scale_x_discrete(labels=stripnum)
    p <- p + lightertheme
    p <- p + xlab("Approaches") + ylab(toupper(metric))
    p <- p + labs(fill="Approaches")
    p <- p + theme(axis.text.x = element_text(size=10, angle = 300, hjust = 0, colour = "black"))
    ggsave(sprintf("VARSELECT_%s_%s.pdf", metric, onesorunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
  },simplify=F)
 },simplify=F)


###### X. Plot rank correlation across the approaches ######

## Get data for plotting from combineandsample ##
allcorr <-  combineandsample$ip$cor

plotcorbox <- function(datasub){
  ggdata <- prepdata(datasub, oneorunif="")
  ggdata <- ggdata[snr != 25]
  ggp <- ggplot(ggdata, aes(x=modelnum, y=scor, fill=model))
  ggp <- ggp + geom_boxplot(outlier.size=1, fatten=0.5)
  p <- plotlayers(plottingdata=ggdata, ymetric="Spearman's Rank Correlation", ggp=ggp)
  p <- p + labs(fill="Approaches")
  p <- p + scale_fill_manual(values=newcolours)
  p <- p + theme(axis.text.x = element_text(size=10, angle = 290, hjust = 0, colour="black"))
  p <- p + theme(axis.title.y = element_text(vjust=1, colour="black"))
  p <- p + scale_x_discrete(label=stripnum)
  ggsave("Correlation_boxplot.pdf", plot=p, width=9, height=6.5, units="in", limitsize=F)
}

# FINAL PLOTS
# Variability in Rank Correlation: Figure 6
plotcorbox(combineandsample$ip$cor)

# You can safely ignore the warning message below.
# An explanation for this warning is here: http://stackoverflow.com/a/20688045
# Warning message:
#   In `[.data.table`(rocfin, , `:=`(snrname, getsnrname(snr)), ) :
#   Invalid .internal.selfref detected and fixed by taking a copy of the whole table so that := can add this new column by reference. At an earlier point, this data.table has been copied by R (or been created manually using structure() or similar). Avoid key<-, names<- and attr<- which in R currently (and oddly) may copy the whole data.table. Use set* syntax instead to avoid copying: ?set, ?setnames and ?setattr. Also, in R<=v3.0.2, list(DT1,DT2) copied the entire DT1 and DT2 (R's list() used to copy named objects); please upgrade to R>v3.0.2 if that is biting. If this message doesn't help, please report to datatable-help so the root cause can be fixed.