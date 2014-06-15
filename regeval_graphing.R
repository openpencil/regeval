######0. This should be the git directory #####
workdir <- "~/yourgitdirectory/regeval"
dir.create(workdir)
setwd(workdir)

###### I. libraries and locales ######
source('./regeval_packages.R')


###### Ib. Data that was used for the simulation ######
load(file="example_dataset.rda")
load(file="example_response.rda") 
numsamples <- nrow(exampledata)

###### II. Load common functions ######
evaluateselection <- function(selectedvars, truevars, allvars){
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

genrescale <- function(value){
  rvalue <- rescale(x=value, to=c(0,1))
  return(rvalue)
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
selectionvector <- function(multicol, method){
#   # #Uncomment for debugging!
#     multicol <- bindpre[model== "lr" & 
#                          snr == "16.00" & 
#                          cutoff == 0.02 & 
#                          simulnum == "1" &
#                          onesunif == "ones",
#                          c("variables", "rescale"), with=F]
  multiorder <- multicol[order(multicol$rescale, decreasing=T),]
  if (method == "ld") {
    # Find the first order lagged differences
    difforder <- -(diff(multiorder$rescale, 1))
    # Find the index corresponding to the maximum lagged difference
    maxdiffindex <- which(max(difforder) == difforder)[1]
    # Build a vector of 1s and 0s. 1s: Selected variables, 0s. Variables not selected.
    select <- ifelse(1:nrow(multiorder) <= maxdiffindex, 1, 0)
  } else if (method == "nbylogp") {
    p <- nrow(multiorder)
    nbylogp <- floor(numsamples/log(p))
    select  <- ifelse(1:nrow(multiorder) <= nbylogp, 1, 0)
  }
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

##### IV. Harvest simulation results #####
# debug
# rdsfile <- allrds[1]
# model <- "lr"
simulationresults <- mclapply(allrds, function(rdsfile){
  snr <- gsub(".*(16\\.00|25|4\\.60)_0.*rds$", "\\1", rdsfile)
  cutoff <- gsub(".*_(0.0\\d)_(\\d*).*.rds", "\\1", rdsfile)
  simulnum <- gsub(".*_(\\d*).*.rds", "\\1", rdsfile)
  onesunif <- gsub(".*(ones|unif).*rds$", "\\1", rdsfile)
  rdslist <- readRDS(rdsfile)
  rdsnames <- setdiff(names(rdslist$betasout), "truebeta")
  ## Gather up inclusion probabilities for all variables
  # model <- "bma"
  modelprob <- sapply(names(rdslist$inclprobs), function(model){
    modelip <- data.table(ldply(rdslist$inclprobs[[model]]))
    setnames(modelip, colnames(modelip), c("variables", "ipval"))
    modelip[,variables:=gsub("^xscaled", "", variables),]
    modelip[,model:=model,]
    return(modelip)
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
},mc.cores=8)
saveRDS(simulationresults, "simulationresults.rds")

##### Combine all simulation results #####
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
    # rescale the inclusion probabilities from 0 to 1 for each combination of model, snr, cutoff, simulnum and one/unif
    # for the stability models, IPs range from 0 to 100, because IPs are summed over all lambdas. However, for
    # all the other models that calculate the stability based on 1 lambda or are BMA, IPs are between 0 and 1.
    # Rescaling everything between 0 and 1 reconciles these two.
    bindpre <- bindup[ ,rescale:=genrescale(value=ipval), by=c("model", "snr", "cutoff", "simulnum", "onesunif")] 
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
    bindup <- bindqx[ ,genthreshold(.SD), by=c("model", "snr", "cutoff", "simulnum", "onesunif"), .SDcols=colinterest]  
    # Save this as an RDS just in case the simulation is massive and you need the results in a jiffy.
    saveRDS(bindup, "rocdata.rds")
    #Reads the saved file
    #bindup <- readRDS("rocdata.rds")
    lengthsubset <- bindup[,list(lenval=length(auc)),by=c("model", "snr", "cutoff", "onesunif")]
    corcol <- c("variables", "rescale")
    bindunif <- bindpre[onesunif=="unif"]
    bindcor <- bindunif[,list(scor=rankcorrelate(.SD)$scor), by=c("model", "snr", "cutoff", "simulnum"), .SDcols=corcol]  
    # bindpre consists of all simulations; unifs and ones
    selectld <- bindpre[,selectionvector(.SD, method="ld"), by=c("model", "snr", "cutoff", "simulnum", "onesunif"), .SDcols=corcol]
    selectnbylogp <- bindpre[,selectionvector(.SD, method="nbylogp"), by=c("model", "snr", "cutoff", "simulnum", "onesunif"), .SDcols=corcol]
  }
  # minimum number of simulations available for any set of conditions.
  minnumsample <- min(lengthsubset$lenval)
  sdcols <- setdiff(colnames(bindup), c("simulnum","model", "snr", "cutoff", "onesunif"))
  
  # From each set of conditions (model, snr, cutoff, onesunif), 
  # sample a minimum of "minnumsample" simulations without replacement 
  # This is so that every set of conditions has the same number of simulations.
  sampledrows <- bindup[,mysamplemulti(.SD, minval=minnumsample, do=T), 
                        by=c("model", "snr", "cutoff", "onesunif"), .SDcols=sdcols]
  if (type == "eval"){
    return(sampledrows)
  } else{
    sampledcorr <- bindcor[,mysamplemulti(.SD, minval=minnumsample, do=T), by=c("model", "snr", "cutoff"), .SDcols=c("scor")]
    sampledld <- selectld[,mysamplemulti(.SD, minval=minnumsample, do=T), 
                                 by=c("model", "snr", "cutoff", "onesunif"), 
                                 .SDcols=c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL")]
    samplednbylogp <- selectnbylogp[,mysamplemulti(.SD, minval=minnumsample, do=T), 
                                 by=c("model", "snr", "cutoff", "onesunif"), 
                                 .SDcols=c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL")]
    return(list(ip=sampledrows, cor=sampledcorr, selectld=sampledld, selectnlp=samplednbylogp))
  }
},simplify=F)

# Ignore this warning.
# Warning messages:
# In ifelse(grepl("(p|m)", varname), as.numeric(gsub(".*[mp](\\d{1,}$)",  ... :
#                                                           NAs introduced by coercion
saveRDS(combineandsample, "combinesample.rds")

############Read in sampled data##############
# just in case you are picking up analysis after a break.
# combineandsample  <- readRDS("./combinesample.rds")


####### V. Get ROC data corresponding to the median AUC #######
# The ROC curve for a binary classification problem plots 
# the true positive rate as a function of the false positive rate. 
# The points of the curve are obtained by sweeping the 
# classification threshold from the most positive classification 
# value to the most negative. 

inclprob <- combineandsample$ip$ip
# calculate mean tpr, fpr, auc and tau for each set of conditions
tprfprtauauc <- grep("tpr|fpr|tau|auc", colnames(inclprob), ignore.case=T, value=T)
meantprfpr <- inclprob[,lapply(.SD, mean), by=list(model, snr, cutoff, onesunif), .SDcols=tprfprtauauc]

# calculate median of auc from big dataset for each combination of things in by below:
aucmedian <- inclprob[,list(aucmed=round(median(auc), 2)), by=list(model, snr, cutoff, onesunif)]
auc_cols <- colnames(aucmedian)
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

#########VI. General graphing code #############
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

prepdata <- function(rocfin, oneorunif){
  if (oneorunif != "") {
    rocfin <- rocfin[onesunif==oneorunif,]
  }
  rocfin[,snrname:=ifelse(grepl("16", snr), "SNR:16.0",
                          ifelse(grepl("25", snr), "SNR:0.25",
                                 ifelse(grepl("4.6", snr), "SNR:04.6", "undefined")))]
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
##########VII. Plot ROC curve##############

#db
#what <- "mean"
#oneunif <- "unif"
##Make extra columns for labels for graphs##
findata <- sapply(list(rocmean, inclprob), function(datasub){
  sapply(unique(datasub$onesunif), function(type){
    fin <- prepdata(datasub, type)
  },simplify=F)
},simplify=F)
names(findata) <- c("mean", "all")

plotroc <- function(what, oneunif){
  rocfin <- findata[[what]][[oneunif]]
  ggp <- ggplot(rocfin, aes(x=FPR, y=TPR, colour=model))
  ggp <- ggp + geom_abline(aes(intercept=0, slope=1), linetype="dashed")
  ggp <- ggp + geom_path()
  p <- plotlayers(plottingdata=rocfin, ymetric="", ggp=ggp)
  # Adjust coordinates of the legend to occupy an empty space in the graph
  p <- p + theme(legend.justification=c(0,0), legend.position=c(0.85,-0.025))
  p <- p + theme(legend.background=element_rect(fill = "transparent"))
  p <- p + xlab("False Positive Rate") + ylab("True Positive Rate")
  p <- p + labs(colour="Approaches")
  p <- p + scale_color_manual(values=newcolours)
  p <- p + scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  ggsave(sprintf("ROC_%s_%s.pdf", what, oneunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
}

# FINAL PLOTS
# Reference: Receiver Operating Characeristics (ROC curves): Figure 2/3 
plotroc(what="mean", oneunif="ones")
plotroc(what="mean", oneunif="unif")


#db
#oneunif <- "unif"
plotauc <- function(oneunif){
  datasub <- prepdata(aucmedian, oneorunif=oneunif)
  ggp <- ggplot(datasub, aes(x=modelnum, y=aucmed, fill=model))
  ggp <- ggp + geom_bar(stat="identity")
  p <- plotlayers(plottingdata=datasub, ymetric="AUC", ggp=ggp)
  #p <- p + theme(legend.justification=c(0,0), legend.position=c(0,0))
  p <- p + theme(legend.background=element_rect(fill = "transparent"))
  p <- p + labs(fill="Approaches")
  p <- p + scale_fill_manual(values=newcolours)
  p <- p + scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + theme(axis.text.x = element_text(size=10, angle = 290, hjust = 0, colour="black"))
  p <- p + scale_x_discrete(label=stripnum)
  ggsave(sprintf("AUC_%s.pdf", oneunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
}

# FINAL PLOTS
# Reference: Area under the ROC curves (median AUC) : Figure 4/5 
plotauc(oneunif="unif")
plotauc(oneunif="ones")


plotaucbox <- function(oneunif){
  rocfin <- findata$all[[oneunif]]
  ggp <- ggplot(rocfin, aes(x=modelnum, y=auc, fill=model))
  ggp <- ggp + geom_boxplot()
  p <- plotlayers(plottingdata=rocfin, ymetric="AUC", ggp=ggp)
  #p <- p + theme(legend.justification=c(0,0), legend.position=c(0,0))
  p <- p + theme(legend.background=element_rect(fill = "transparent"))
  p <- p + labs(fill="Approaches")
  p <- p + scale_fill_manual(values=newcolours)
  p <- p + scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1), limits=c(0,1))
  p <- p + theme(axis.text.x = element_text(size=10, angle = 290, hjust = 0, colour="black"))
  p <- p + scale_x_discrete(label=stripnum)
  ggsave(sprintf("AUC_boxplot_%s.pdf", oneunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
}
# FINAL PLOTS
# Reference: Variability in AUC: Figure 6
plotaucbox("unif")
plotaucbox("ones")

#######VIII. Plots other than ROC#########
allbeta <-combineandsample$eval
allldselect <- combineandsample$ip$selectld
allnlpselect <- combineandsample$ip$selectnlp
allcorr <-  combineandsample$ip$cor

sdcols <- setdiff(colnames(allbeta), c("simulnum","model", "snr", "cutoff", "onesunif"))
mymedian <- function(x){median(x, na.rm=T)}
evalmedian <- allbeta[,lapply(.SD, mymedian), by=c("model", "snr", "cutoff", "onesunif"), .SDcols=sdcols]
evalmedian[,method:="MODEL",]
#method of lagged differences
selectldmedian <- allldselect[,lapply(.SD, mymedian), by=c("model", "snr", "cutoff", "onesunif"), .SDcols=sdcols]
#method: n/log(p)
selectnlpmedian <- allnlpselect[,lapply(.SD, mymedian), by=c("model", "snr", "cutoff", "onesunif"), .SDcols=sdcols]
selectldmedian[,method:="LD",]
selectnlpmedian[,method:="NLP",]
selectmedian <- rbind(selectldmedian, selectnlpmedian)
cormedian <- allcorr[,lapply(.SD, mymedian), by=c("model", "snr", "cutoff"), .SDcols="scor"]


## New columns for plotting ##
mmdata <- list(emed=evalmedian, smed=selectmedian)
mmdata_ones <- sapply(mmdata, function(datasub){
  pdata <- prepdata(rocfin=datasub, oneorunif="ones")
  return(pdata)
},simplify=F)
mmdata_unifs <- sapply(mmdata, function(datasub){
  pdata <- prepdata(rocfin=datasub, oneorunif="unif")
  return(pdata)
},simplify=F)


#db
#type <- "med"
#metric <- "FPR"
ou <- c("ones", "unifs")
#onesorunif <- ou[2]

# FINAL PLOTS
# Reference: FScore: Comparison of Variable Selection Strategies Figure 7-12
sapply(ou, function(onesorunif){
  x <- "smed" #var selection
  y <- "emed" #models without IP
  print(sprintf("mmdata_%s", onesorunif))
  datasub <- get(sprintf("mmdata_%s", onesorunif)) 
  mergeddata <- rbind(datasub[[x]], datasub[[y]])
  mergemelt <- melt(mergeddata, measure.vars=c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL"))
  #metric <- "FSCORE"
  twographs <- sapply(c("FPR", "TPR", "FSCORE", "PRECISION", "RECALL"), function(metric){
    ipdata <- mergemelt[variable==metric & model %in% unique(datasub$smed$model)]   
    p <- ggplot(ipdata, aes(x=modelnum, y=value, fill=method))
    p <- p + geom_bar(stat="identity", position="dodge")
    p <- p + facet_grid(snrname~cutoff, labeller=modifylabel)
    p <- p + lightertheme
    p <- p + xlab("Approaches") + ylab(toupper(metric))
    p <- p + scale_fill_discrete(name="Selection", labels=c("LD", "n/log(p)"))
    p <- p + scale_x_discrete(labels=stripnum)
    p <- p + theme(legend.text=element_text(face="italic"))
    p <- p + theme(axis.text.x = element_text(size=10, angle = 300, hjust = 0, colour = "black"))
    ggsave(sprintf("LDVSNLP_%s_%s.pdf",metric, onesorunif), width=9, height=6.5, units="in", limitsize=F)
    
    allmodels <- mergemelt[grepl("LD|MODEL", method) & variable==metric,]
    p <- ggplot(allmodels, aes(x=modelnum, y=value, fill=model))
    p <- p + geom_bar(stat="identity")
    p <- p + facet_grid(snrname~cutoff, labeller=modifylabel)
    p <- p + scale_fill_manual(values=newcolours)
    p <- p + scale_x_discrete(labels=stripnum)
    p <- p + lightertheme
    p <- p + xlab("Approaches") + ylab(toupper(metric))
    p <- p + labs(fill="Approaches")
    p <- p + theme(axis.text.x = element_text(size=10, angle = 300, hjust = 0, colour = "black"))
    ggsave(sprintf("LDVSMODEL_%s_%s.pdf", metric, onesorunif), plot=p, width=9, height=6.5, units="in", limitsize=F)
  },simplify=F)
  },simplify=F)  


# FINAL PLOTS
# Figure 13 Spearman’s Correlation: Variable Ranking, 
plotcorbar <- function(datasub){
  ggdata <- prepdata(datasub, oneorunif="")
  ggp <- ggplot(ggdata, aes(x=modelnum, y=scor, fill=model))
  ggp <- ggp + geom_bar(stat="identity")
  p <- plotlayers(plottingdata=ggdata, ymetric="Spearman's Rank Correlation", ggp=ggp)
  p <- p + scale_fill_manual(values=newcolours)
  p <- p + scale_x_discrete(labels=stripnum)
  p <- p + theme(axis.text.x = element_text(size=10, angle = 300, hjust = 0, colour = "black"))
  ggsave("SpearmanCorr_VarRanking.pdf", plot=p, width=9, height=6.5, units="in", limitsize=F) 
}
plotcorbar(cormedian)
