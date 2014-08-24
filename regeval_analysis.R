###### 0. This should be the git directory ######
workdir <- "~/regeval"
dir.create(workdir)
setwd(workdir)

###### I. libraries and locales ######
source('./regeval_packages.R')
source('./regeval_algorithms.R')
source('./regeval_corrplot.R')
source('./regeval_colored_dendrogram.R')
source('./regeval_colorlegend.R')

###### Ib. Load and process data  ######
load(file="example_dataset.rda")
load(file="example_response.rda") 

## Process y ## 
yscaled <- scale(exampleresponse$response)

## Function for processing x ##
processx <- function(countdata, refindex) {
  ## Remove rows and columns with zero rowSums and colSums ##
  ## Convert sequence counts to log probabilities ##
  logprobx <-  logprobclean(countdata)
  ## Address sum-to-1 redundancy ## 
  ## Select reference variable. (Equation 3) ##
  referencevarindex <- order(colSums(countdata), decreasing=T)[refindex]
  ## Drop reference variable from the design matrix. ##
  nonredundantmatrix <- logprobx[,-referencevarindex]
  ## Scale the design matrix to make mean 0 and variance 1 for all columns
  nonredundantxscaled <- apply(nonredundantmatrix, 2, scale)
  return(nonredundantxscaled)
}

###### II. Run analysis ######
modelresults <- sapply(c(1, 2), function(element){
  xscaled <- processx(exampledata, element)
  results <- runmodels(xscaled, yscaled, seed=101, iterations=10000,
                       oracle=rep(0, ncol(xscaled)))
  return(results)
}, simplify=F)

# Ignore these non-consequential memory-based warning messages coming from the doMC package.
# rsession(5276) malloc: *** error for object 0x7fe72d67dd40: pointer being freed was not allocated
# *** set a breakpoint in malloc_error_break to debug
# This is not a serious error. Ignore
# In predict.lm.spike(object = train_path[[fold]][[feature]],  ... :
# Implicit intercept added to newdata

names(modelresults) <- c("ref1", "ref2")

# Save model results
saveRDS(modelresults, "modelresults.RDS")


###### III. Average models across two reference variable baselines ######
meltit <- function(datasub, model){
  dmelt <- melt(datasub$inclprobs[[model]])
  dmelt$var <- rownames(dmelt)
  return(dmelt)
}

averagemodels <- function(){
  data1 <- modelresults$ref1
  data2 <- modelresults$ref2
  ipaverage <- sapply(names(data1$inclprobs), function(model){
    d1 <- meltit(data1, model)
    d2 <- meltit(data2, model)
    d1d2 <- merge(d1, d2, by="var", all=T)
    d1d2$mean <- rowMeans(d1d2[,c("value.x", "value.y")], na.rm=TRUE)
    meanvector <- d1d2$mean
    names(meanvector) <- d1d2$var
    means <- meanvector[grep("Intercept", names(meanvector), invert=T)]
    return(means)
  }, simplify=F)
  varintersect <- sapply(names(data1$varsout), function(model){
    unionvars <- union(data1$varsout[[model]], data2$varsout[[model]])
    unions <- grep("Intercept", unionvars, invert=T, value=T)
    return(unions)
  }, simplify=F)
  return(list(ip=ipaverage, var=varintersect))
}

## Get the average model results ##  
averageresults <- averagemodels()

###### IV. Combine non-ip and ip method findings  ######

compileresults <- function(){
  dfip <- averageresults$ip
  dfvar <- averageresults$var
  # Arrange IP values in descending order
  ipvals <- sapply(names(dfip), function(model){
    ip <- dfip[[model]]
    names(ip) <- gsub("^x", "", names(ip))
    ipdf <- data.frame(ip)
    if (grepl("(ss|sr)", model)){
      ipdf$ip <- ipdf$ip/100
    }
    ipdf <- ipdf[order(ipdf$ip, decreasing=T),,drop=F]
    ipdf$num <- 1:nrow(ipdf)
    ipdf$vars <- rownames(ipdf)
    return(data.table(ipdf))
  }, simplify=F)
  allip <- ldply(ipvals)
  colnames(allip)[1] <- "model" 
  allip$model <- toupper(allip$model)
  # Package up the non-IP methods
  varselect <- sapply(setdiff(names(dfvar), "truevars"), function(model){
    dv <- data.frame(vars=dfvar[[model]], ip=1, num=1)
    return(dv)
  }, simplify=F)
  varvalues <- ldply(varselect)
  setnames(varvalues, ".id", "model")
  varvalues$model <- toupper(varvalues$model)
  combomodels <- ldply(list(ip=allip, var=varvalues))
  modelcast <- cast(data=combomodels, formula=model~vars, fill=0, value="ip")
  rownames(modelcast) <- modelcast$model
  allmodels <- modelcast[,-1]
  return(list(modelcast=allmodels, modelips=allip))
}
allmodelscompiled <- compileresults()

###### V. Cluster models together based on their inclusion probability profiles ######

## Generate colours for models ##
rgb <- c("#800026", "#860026", "#8D0026", "#940026", "#9A0026", "#A10026", 
         "#A80026", "#AE0026", "#B50026", "#BC0026", "#C00225", "#C40523", 
         "#C80822", "#CD0B21", "#D10D20", "#D5101F", "#D9131E", "#DD161D", 
         "#E1191C", "#E51E1D", "#E7231E", "#EA2920", "#ED2F21", "#F03523", 
         "#F23A24", "#F54026", "#F84627", "#FA4B29", "#FC522B", "#FC592D", 
         "#FC602F", "#FC6731", "#FC6D33", "#FC7435", "#FC7B37", "#FC8239", 
         "#FC893B", "#FD8F3C", "#FD933E", "#FD9740", "#FD9B42", "#FD9F43", 
         "#FDA345", "#FDA747", "#FDAB49", "#FDAF4A", "#FEB34D", "#FEB752", 
         "#FEBC56", "#FEC05B", "#FEC460", "#FEC864", "#FECD69", "#FED16D", 
         "#FED572", "#FED977", "#FEDB7B", "#FEDD80", "#FEE084", "#FEE289", 
         "#FEE48E", "#FEE692", "#FEE897", "#FEEB9B", "#FFEDA0", "#FFEFA5", 
         "#FFF1AA", "#FFF3AF", "#FFF5B3", "#FFF7B8", "#FFF9BD", "#FFFBC2", 
         "#FFFDC7", "#FFFFCC",
         "#FFFFFF",
         "#FFFFD9", "#FDFED5", "#FCFDD2", "#FAFDCF", "#F9FCCC", "#F7FCC8", 
         "#F6FBC5", "#F4FBC2", "#F3FABF", "#F1F9BB", "#F0F9B8", "#EFF8B5", 
         "#EDF8B2", "#EBF7B1", "#E8F6B1", "#E4F4B1", "#E1F3B1", "#DEF2B2", 
         "#DBF1B2", "#D8EFB2", "#D5EEB2", "#D2EDB3", "#CFECB3", "#CCEBB3", 
         "#C9E9B3", "#C5E8B4", "#BFE6B4", "#B9E3B5", "#B4E1B5", "#AEDFB6", 
         "#A8DDB6", "#A2DAB7", "#9CD8B8", "#96D6B8", "#91D4B9", "#8BD1B9", 
         "#85CFBA", "#7FCDBA", "#7ACBBB", "#75C9BC", "#70C7BD", "#6BC5BD", 
         "#66C3BE", "#61C2BF", "#5CC0C0", "#57BEC0", "#52BCC1", "#4DBAC2", 
         "#48B8C2", "#43B6C3", "#3FB4C3", "#3CB1C3", "#39AEC3", "#36ABC2", 
         "#33A8C2", "#30A5C2", "#2EA2C1", "#2B9FC1", "#289CC1", "#2599C0", 
         "#2296C0", "#1F93C0", "#1D90BF", "#1D8CBD", "#1D88BB", "#1E84B9", 
         "#1E7FB7", "#1F7BB6", "#1F77B4", "#1F73B2", "#206FB0", "#206BAE", 
         "#2167AC", "#2163AA")
gencolours <- colorRampPalette(rgb)
twocolour <- c("#FFFFFF", "#ADDD8E")
gencolor <- colorRampPalette(twocolour)

newcolours <- c("#d73027", "#f46d43", "#e08214", "#66bd63", 
                "#136c39", "#1d91c0", "#19467e", "#c2a5cf", 
                "#9970ab", "#dd3497", "#ae017e")
names(newcolours) <- c("BMA", "BMAC", "ENC", "LR", 
                       "LRC", "LS", "LSC", "PR", 
                       "PS", "SR", "SS")

rankdistribution <- function(){
  datasub <- allmodelscompiled$modelips
  p <- ggplot(datasub, aes(x=num, y=ip, colour=model))
  p <- p + geom_line()
  p <- p + scale_colour_manual(values=newcolours, name="Approaches")
  p <- p + theme(legend.justification=c(0,0), legend.position=c(0.7,0.05))
  p <- p + theme(legend.background=element_rect(fill = "transparent"))
  p <- p + xlab("Rank of the genera") + ylab("Inclusion probability")
  p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  p <- p + scale_x_continuous(breaks=c(1, seq(25, max(datasub$num), by=25)), 
                                labels=c(1, seq(25, max(datasub$num), by=25)))
  p <- p + lightertheme
  p <- p + ggtitle("Inclusion Probability Profiles")
  ggsave("IP_rank_distribution.pdf", plot=p, width=6.7, height=3.7, units="in", limitsize=F)
}

# FINAL PLOTS
# Comparison of inclusion probability profiles across all the eleven approaches: Figure 9a and Figure 9b (Top panels)
rankdistribution()

corrplotit <- function(){
  rankcor <- 1/(1 + as.matrix(dist(allmodelscompiled$modelcast, "euclidean")))
  colnames(rankcor) <- rownames(allmodelscompiled$modelcast)
  rownames(rankcor)  <- rownames(allmodelscompiled$modelcast)
  pdf("ALLMODELS_corrcomparison.pdf", width=12, height=8)  
  mycorrplot(method="color", corr=rankcor, col=gencolours(592),
             cl.pos="b",cl.length=5, cl.cex=1,
             tl.col="black", cex.axis=0.5, tl.cex=1,
             type="lower", diag=T, bg="transparent",
             addgrid.col="#f4f4f5", addCoefasPercent=T, addCoef.col="black",
             mar=c(0,0,0,0), order="hclust", hclust.method="average")
  dev.off()
  clustorder <- hclust(dist(allmodelscompiled$modelcast, method="euclidean"), method="average")
  clustorder$labels <- rownames(allmodelscompiled$modelcast)
  pdf("ALLMODELS_clusterdendrogram.pdf", width=7.7, height=3)  
  A2Rplot(clustorder, k=11, boxes=F, main="",
          col.up="gray50", col.down=newcolours, 
          lwd.up=2, 
          lty.up="solid")
  dev.off()
}

# FINAL PLOTS
# Comparison of inclusion probability profiles across all the eleven approaches: Figure 9a and Figure 9b (Bottom panels)
corrplotit()

##### VI. Compare variable selections #####
# Variable selection using the Method of lagged differences(LD) ###
# df <- bactresults$ip$bma

selectvars <- function(ipvector){
  dforder <- ipvector[order(ipvector, decreasing=T)]
  # Find the first order lagged differences
  difforder <- -(diff(dforder, 1))
  # Find the index corresponding to the maximum lagged difference
  maxdiffindex <- which(max(difforder) == difforder)[1]
  # Build a vector of 1s and 0s. 1s: Selected variables, 0s. Variables not selected.
  select <- ifelse(1:length(ipvector) <= maxdiffindex, 1, 0)
  names(select) <- names(dforder)
  selectedvars <- names(dforder[which(select==1)])
  return(selectedvars)
}

# Isolate only the selected variables
compileselected <- function(){
  dfip <- averageresults$ip
  dfvar <- averageresults$var
  # Select variables by method of lagged differences and grab their inclusion probabilities
  ipselect <- sapply(names(dfip), function(model){
    dm <- dfip[[model]]
    sv <- selectvars(dm)
    withip <- dm[sv]
    names(withip) <- gsub("^x|xscaled", "", names(withip))
    if(model %in% c("ss", "sr")){
      withip <- withip/100
    }
    withmelt <- melt(withip)
    ips <- data.frame(var=rownames(withmelt), ip=withmelt)
    return(ips)
  },simplify=F)
  ipvalues <- ldply(ipselect)
  # For non-ip assigning methods, get all the variables and assign them an IP of 1.0
  varselect <- sapply(setdiff(names(dfvar), "truevars"), function(model){
    dv <- data.frame(var=dfvar[[model]], value=1)
    return(dv)
  }, simplify=F)
  varvalues <- ldply(varselect)
  # Combine ip and non-ip methods
  allvars <- rbindlist(list(ipvalues, varvalues))
  datacast <- cast(allvars, formula=.id~var, value="ip", fill=0)
  rownames(datacast) <- datacast$.id
  return(datacast)
}

## Generate selected variable matrix ##
selectedvarmatrix <- compileselected()

labelit <- function(feedin){
  feedout <- gsub("^.*\\d{2}_", "", feedin)
  return(feedout)
}

comparevarselections <- function(){
  # Matrix of only the selected variables
  ipcast  <- selectedvarmatrix
  # Find maximum number of selected variables for graphing dimensions
  numvar <- ncol(ipcast)
  # Matrix of all the variables
  allcast <- allmodelscompiled$modelcast
  
  csums <- sapply(2:ncol(ipcast), function(numcol){
    sum(as.numeric((ipcast[,numcol])), na.rm=T)
  })
  ggorder <- ipcast[,c(1, (order(csums, decreasing=T)+1))]
  rownames(ggorder) <- ggorder$.id
  ggorder <- ggorder[,-1]
  # To make variables with the highest IP appear on top
  ggorder <- cbind(ggorder[,rev(colnames(ggorder))], Approach=rep(-1, nrow(ggorder)))
  ### Start code for heatmap model comparison ####
  rownames(ggorder) <- paste(sprintf("%02d", seq(1, nrow(ggorder), 1)), rownames(ggorder), sep="_")
  colnames(ggorder) <- paste(sprintf("%02d", seq(ncol(ggorder), 1, -1)), colnames(ggorder), sep="_")
  ggorder$model <- rownames(ggorder)
  ggframe <- melt(data.frame(ggorder), id.vars="model")
  ggframe$modelvar <- paste(toupper(gsub("\\d{2}_", "", ggframe$model)), 
                            gsub("X\\d{2}_", "", ggframe$variable), sep="_")
  allcast$model <- rownames(allcast)
  allmelt <- melt(data.frame(allcast), id.vars="model")
  allmelt$modelvar <- paste(allmelt$model, allmelt$variable, sep="_")
  merged <- merge(allmelt, ggframe, by="modelvar", all.y=T)
  merged$percent <- ifelse(merged$value.y== -1, "", round(merged$value.x*100, 0))
  merged$colorval <- ifelse(merged$value.y== -1, 0, merged$value.y)
  mergeddt <- data.table(merged)
  # one slice of the data for retrieving all the model names
  mergeddt <- mergeddt[grepl("X02", variable.y)]
  mergeddt[,modellab:=toupper(gsub("^\\d+_", "", model.x))]
  p <- ggplot(merged, aes(x=model.y, y=variable.y))
  p <- p + geom_tile(aes(fill = colorval), colour = "white")
  p <- p + scale_fill_gradientn(colours = gencolor(2), guide="none")
  p <- p + scale_y_discrete(labels=labelit)
  p <- p + scale_x_discrete(labels=labelit, breaks=NULL)
  p <- p + theme(axis.text.x = element_blank(), 
                 axis.text.y = element_text(colour = "black", vjust=0.1, size=12))
  p <- p + geom_text(aes(label=percent))
  p <- p + xlab("") + ylab("")
  p <- p + geom_text(data = mergeddt,
                     aes(label = modellab), y = ncol(ggorder)-1, size = 3.5, fontface="bold")
  p <- p + ggtitle("Influential variables selected across the 11 approaches")
  calcheight <- 0.28 * (numvar+1)
  ggsave("ALLMODELS_selectedvarscomparison.pdf", plot=p, width=8, height=calcheight, units="in", limitsize=F)
  }

# FINAL PLOTS
# Summary of influential variables selected across all models: Figure 7 and Figure 8 ##
comparevarselections()
