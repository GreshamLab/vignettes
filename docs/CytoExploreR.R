#CytoExploreR.R
#Author David Gresham
#Started 12/02/2020

#This is my attempt to use CytoExploreR to perform a complete analysis of our data.

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

library(devtools)
library(tidyverse)

#library(BiocManager)
#install(“cytolib”, “flowCore”, “flowWorkspace”, “openCyto”)

#Then you can install CytoExploreR from GitHub:

#devtools::install_github("DillonHammill/CytoExploreRData")
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(CytoExploreRData)

#Setup experiment
gating_set <- cyto_setup(path="/Users/david/Projects/data/flowdata/Accuri", select="fcs")

#To interactively edit the experiment details
cyto_details_edit(gating_set)

#Transform the data using cyto_transformer()
transformed <- cyto_transformer_logicle(gating_set)
FSC_SSC_transformed <- cyto_transformer_logicle(gating_set,
                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H"))
trans <- cyto_transformer_combine(transformed,
                                  FSC_SSC_transformed)
trans

transformed_accuri <- cyto_transform(gating_set,
                                    trans = trans)

#Transform the data using cyto_transform()
#the cyto_transform(0 function does not transform the FSC and SSC channels)
gating_set_transformed <- cyto_transform(gating_set,
                                         type="log")

#gate cells
cyto_gate_draw(transformed_accuri,
               alias = "Cells2",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "Accuri_test_data.csv",
               type = "boundary")

cyto_gate_draw(transformed_accuri,
               parent = "Cells2",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-A"),
               gatingTemplate = "Accuri_test_data.csv"
               )

cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = "One_Copy",
               channels = c("FSC-A","FL1-A"),
               gatingTemplate = "Accuri_test_data.csv"
               )


#if you need to redraw a gate you need to use cyto_gate_edit()
cyto_gate_edit(transformed_accuri,
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "Accuri_test_data.csv",
               type = "boundary")

