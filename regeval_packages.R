##### I. Install necessary packages #####
installpackages<- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  else if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  } 
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) {# alternative to library
      require(need_to_attach[i], character.only = TRUE)
    }
  }
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}

installpackages(c("ggplot2", "plyr", "Hmisc", "reshape", "reshape2", 
                  "scales", "gridExtra", "glmnet","doMC", "quadrupen",
                  "RColorBrewer", "mice", "mixOmics", "corrplot", 
                  "igraph", "data.table", "dirmult", "mht", "digest"))

##### II. Set locales #####
Sys.setlocale("LC_CTYPE","en_US.UTF-8")
Sys.setlocale("LC_TIME","en_US.UTF-8")
Sys.setlocale("LC_COLLATE","en_US.UTF-8")
Sys.setlocale("LC_MONETARY","en_US.UTF-8")
Sys.setlocale("LC_MESSAGES","en_US.UTF-8")


##### III. Favourite ggplot theme #####
lightertheme <- theme( panel.background = element_rect(fill = "#f5f5f4", colour = "#fbfbfc"),
                       panel.border = element_rect(colour = "grey80", linetype="solid", fill=NA),
                       panel.grid.major = element_line(colour = "grey90", size=0.08),
                       panel.grid.minor = element_line(colour = "grey90", size=0.08),
                       panel.margin = unit(0.3, "lines"),
                       plot.background=element_rect(fill="transparent"),
                       plot.margin=unit(c(0,0,0,0),"mm"),
                       text = element_text(family="Helvetica", colour="black", size=12),
                       plot.title = element_text(size = 12),
                       strip.background=element_rect(fill="transparent", color="#ffffff", size=1),
                       strip.text.x=element_text(size=12, colour="black"),
                       strip.text.y=element_text(size=12, colour="black"),
                       axis.text.x=element_text(size=12, colour="black"),
                       axis.text.y=element_text(size=12, colour="black"),
                       legend.background=element_rect(fill = "transparent"),
                       legend.background=element_blank(),
                       legend.key=element_blank())