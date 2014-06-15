######0. This should be the git directory #####
workdir <- "~/yourgitdirectory/regeval"
dir.create(workdir)
setwd(workdir)

###### I. libraries and locales ######
source('./regeval_packages.R')
source('./regeval_algorithms.R')
source("./regeval_corrplot.R")
library("BoomSpikeSlab")

###### Ib. Load and process data  ######
load(file="example_dataset.rda")
load(file="example_response.rda") 

logprobx <-  logprobclean(exampledata)
yscaled <- scale(exampleresponse$response)
xscaled <- apply(logprobx, 2, scale)

#####II. Run analysis #####
modelresults <- runmodels(xscaled, yscaled, seed=101, iterations=10000, oracle=rep(0, ncol(xscaled)))
saveRDS(modelresults, "modelresults.RDS")

#####III. Generate LD variable selections #####
# Variable selection using the Method of lagged differences(LD) ###
selectionvector <- function(df){
  dforder <- df[order(df$ele, decreasing=T),]
  # Find the first order lagged differences
  difforder <- -(diff(dforder$ele, 1))
  # Find the index corresponding to the maximum lagged difference
  maxdiffindex <- which(max(difforder) == difforder)[1]
  # Build a vector of 1s and 0s. 1s: Selected variables, 0s. Variables not selected.
  select <- ifelse(1:nrow(df) <= maxdiffindex, 1, 0)
  names(select) <- dforder$var
  selectedvars <- dforder$var[which(select==1)]
  return(selectedvars)
}

# Just in case you are picking off where you left earlier.
# modelresults <- readRDS("modelresults.RDS")
#db
inclprobs <- sapply(modelresults$inclprobs, function(ele){
  names(ele) <- gsub("^x", "", names(ele))
  eledf <- data.frame(ele)
  eledf$var <- rownames(eledf)
  eledf$ele <- rescale(eledf$ele, to=c(0,1))
  selectedvars <- selectionvector(eledf)
  selectip <- data.frame(vars=selectedvars, ip=eledf[selectedvars, "ele"])
  return(selectip)
}, simplify=F)


#ele <- "enc"
nonip <- sapply(names(modelresults$varsout), function(ele){
  if (grepl("true", ele)){
    return(c())
  } else {
    datasub <- modelresults$varsout[[ele]]
    df <- data.frame(vars=datasub, ip=1)
  }
  return(df)
}, simplify=F)

allip <- ldply(inclprobs)
allnonip <- ldply(nonip)
iplist <-  rbind(allip, allnonip)
ipcast <- cast(iplist, formula=.id~vars, value="ip", fill=0)



########IV. Graphing ########

#####a. Colours######
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

########b. Construct heatmap and correlation plots################
labelit <- function(feedin){feedout <- toupper(gsub("^.*\\d{2}_", "", feedin))}

makeggplot <- function(ipcast){
  csums <- sapply(2:ncol(ipcast), function(numcol){
    sum(as.numeric((ipcast[,numcol])), na.rm=T)
  })
  ggorder <- ipcast[,c(1, (order(csums, decreasing=T)+1))]
  rownames(ggorder) <- ggorder$.id
  ggorder <- ggorder[,-1]
  ### Start code for heatmap model comparison ####
  rownames(ggorder) <- paste(sprintf("%02d", seq(1, nrow(ggorder), 1)), rownames(ggorder), sep="_")
  colnames(ggorder) <- paste(sprintf("%02d", seq(1, ncol(ggorder), 1)), colnames(ggorder), sep="_")
  ggorder$model <- rownames(ggorder)
  ggframe <- melt(data.frame(ggorder), id.vars="model")
  ggframe$percent <- round(ggframe$value*100, 0)
  
  #FINAL PLOT
  #Reference: Figure 14, Figure 16
  p <- ggplot(ggframe, aes(x=model, y=variable))
  p <- p + geom_tile(aes(fill = value), colour = "white")
  p <- p + scale_fill_gradientn(colours = gencolor(2), guide="none")
  p <- p + theme_grey(base_size = 12) 
  p <- p + scale_y_discrete(labels=labelit)
  p <- p + scale_x_discrete(labels=labelit)
  p <- p + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))
  p <- p + geom_text(aes(label=percent))
  p <- p + xlab("") + ylab("")
  ggsave("All_Models_Comparisons.pdf", plot=p, width=6, height=8, units="in", limitsize=F)

  #### Start correlation plot
  #FINAL PLOT
  #Reference: Figure 15, Figure 17
  xinterest <- xscaled[,unique(tolower(labelit(ggframe$variable)))]
  corx <- cor(x=xinterest, method="spearman")
  pdf("Selected_Var_Correlations.pdf", width=15, height=12)
  mycorrplot(method="color", corr=corx, col=gencolours(10000),
           cl.pos="b",cl.length=5, cl.cex=1,
           tl.col="black", cex.axis=0.5, tl.cex=1,
           type="lower", diag=F, bg="transparent",
           addgrid.col="#f4f4f5", addCoefasPercent=T, addCoef.col="black",
           mar=c(0,0,0,0), order="FPC") 
  dev.off()
}
makeggplot(ipcast)



