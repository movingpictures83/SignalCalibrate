## ----para, echo = FALSE, results='hide'---------------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(dev="png",fig.show="hold",
               fig.width=4,fig.height=4.5,fig.align="center",
               message=FALSE,collapse=TRUE)

## ----library------------------------------------------------------------------
library(rnaseqcomp)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {

## ----data---------------------------------------------------------------------
# load the dataset in this package
quant = list(rsem=as.matrix(read.csv(paste(pfix, parameters["rsem", 2], sep="/"), check.names=FALSE)), fluxcapacitor=as.matrix(read.csv(paste(pfix, parameters["fluxcapacitor", 2], sep="/"), check.names=FALSE)))
meta = read.csv(paste(pfix, parameters["meta", 2], sep="/"))
samp = read.csv(paste(pfix, parameters["samp", 2], sep="/"))

## ----meta---------------------------------------------------------------------
condInfo <- factor(samp$condition)
repInfo <- factor(samp$replicate)
evaluationFeature <- rep(TRUE, nrow(meta))
calibrationFeature <- meta$house & meta$chr == 'chr1'
unitReference <- 1

## ----filter-------------------------------------------------------------------
dat <- signalCalibrate(quant, condInfo, repInfo, evaluationFeature,
     calibrationFeature, unitReference, 
     calibrationFeature2 = calibrationFeature)
#class(dat)
#show(dat)

saveRDS(dat, outputfile)

## ----sd-----------------------------------------------------------------------
#plotSD(dat,ylim=c(0,1.4))

## ----nonexpplot---------------------------------------------------------------
#plotNE(dat,xlim=c(0.5,1))

## ----tx2----------------------------------------------------------------------
#plot2TX(dat,genes=meta$gene,ylim=c(0,0.6))

## ----diffroc------------------------------------------------------------------
#plotROC(dat,meta$positive,meta$fcsign,ylim=c(0,0.8))

## ----difffc-------------------------------------------------------------------
#meta$fcsign[meta$fcstatus == "off.on"] <- NA
#plotFC(dat,meta$positive,meta$fcsign,ylim=c(0,1.2))

}
