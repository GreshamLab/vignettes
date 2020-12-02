#gating.R

#list.files("/Users/david/Projects/data/flowdata", recursive = T)
library(ncdfFlow)
library(flowCore)
library(flowViz)
library(ggcyto)
library(ggforce)
library(tidyverse)
library(ggridges)
library(readxl)
library(flowWorkspace)
library(openCyto)
library(flowStats)


devtools::install_github("DillonHammill/CytoRSuiteData")

dir = '/Users/david/Projects/data/flowdata/Accuri'
files <- sort(list.files(path=dir,pattern = ".fcs", full.names=TRUE))
flowData <- read.ncdfFlowSet(files=files, pattern=".fcs", alter.names = TRUE)

#needs to differ depending on sample sheet file typ
#sample.sheet <- read_excel(paste(path=dir,"GAP1_multicheck_20May2019.xlsx", sep="/"))
SampleSheet <- read_csv("~/Projects/Data/flowdata/Accuri/SampleSheet.csv")

flowData

sampleNames(flowData) <- SampleSheet$Strain

pData(flowData)$Well <- SampleSheet$Well
pData(flowData)$Strain <- SampleSheet$Strain
pData(flowData)$Genotype <- SampleSheet$Genotype
pData(flowData)$Ploidy <- SampleSheet$Ploidy
pData(flowData)$Media <- SampleSheet$Media
pData(flowData)


plot(flowData[[1]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,1e6),smooth=F)
debris.gate <- locator(100, type='o', col='red')
gm.2 <- matrix(,length(debris.gate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- debris.gate$x
gm.2[,2] <- debris.gate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#Look at the gating on the controls
ggcyto(flowData[[1]], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + geom_gate(pg.nondebris)


