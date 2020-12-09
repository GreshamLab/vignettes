#CytoExploreR.R
#Author David Gresham
#Started 12/02/2020

#This code uses CytoExploreR to perform an analysis of our CNV reporter data.

#load library requirements
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


#Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#install(“cytolib”, “flowCore”, “flowWorkspace”, “openCyto”)

#Then you can install CytoExploreR from GitHub:

#devtools::install_github("DillonHammill/CytoExploreRData")
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(CytoExploreRData)

#Setup experiment using cyto_setup()
#This function will read .fcs files to a cutoset which are then added to a GatingSet
#in generating the cytoset a experiment details csv files is created.  Additional columns can be added by right clicking on a mac and inserting columns
#and a experiment markers csv file is created
accuri_gating_set <- cyto_setup(path="/Users/david/Projects/data/flowdata/Accuri", select="fcs")

#To interactively edit the experiment details
cyto_details_edit(accuri_gating_set)

#Transform the data using a logicle transformation.
#The first step applies the transformation using cyto_transformer_logicle()
#By default, this function does not apply the transformation to FSC and SSC, so we do that in a second step.
#we then need to combine the transformation using cyto_transformer_combined()
#finally we apply the transformation to all the data using cyto_transform
accuri_transformed <- cyto_transformer_logicle(accuri_gating_set)
accuri_FSC_SSC_transformed <- cyto_transformer_logicle(accuri_gating_set,
                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H"))

accuri_combined_transformed <- cyto_transformer_combine(accuri_transformed,
                                  accuri_FSC_SSC_transformed)

transformed_accuri <- cyto_transform(accuri_gating_set,
                                    trans = accuri_combined_transformed)

#####Gating cells.
#To gate cells we use the interactive function of Cytoexplorer
#The details of the gating are maintained in a .cvs file, which must be called
#The gating is done in a hierarchical manner, so that there is a parent and a child for each gate.
#By calling the entire gating set all cells from all fcs files are plotted
#individual files can be called as well

#The first gate defines cells based on forward scatter and side scatter
cyto_gate_draw(transformed_accuri,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv",
              )

#the next gate defines the singlets based on forward scatter height and width
cyto_gate_draw(transformed_accuri,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv"
               )

#####################

#To generate the gates based on individuals samples, we need to use the appropriate individual samples

#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
cyto_gate_draw(transformed_accuri[1],
               parent = "Single_cells",
               alias = "No GFP",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating2.csv"
               )

cyto_gate_edit(transformed_accuri,
               parent = "Single_cells",
               alias = "No GFP",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv"
)

cyto_gate_draw(transformed_accuri[2],
               parent = "Single_cells",
               alias = "One copy GFP",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv"
)

cyto_gate_draw(transformed_accuri[3],
               parent = "Single_cells",
               alias = "Two copy GFP",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv"
)

###in order to visualize an existing gate overlaied on cells

#first we need to extract the cells to be overlaied
negative <- cyto_extract(transformed_accuri, "Single_cells")[[1]]

#the overlay argument plots the cells as gray cells
cyto_gate_draw(transformed_accuri[3],
               parent = "Single_cells",
               alias = c("Neg","Two"),
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv",
               overlay=negative
)

#if you need to redraw a gate you need to use cyto_gate_edit()
cyto_gate_edit(transformed_accuri,
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "Accuri_test_data.csv",
               type = "boundary")


###Drawing a tree of the gated samples
cyto_plot_gating_tree(transformed_accuri[[3]],
                      stat="freq")


#Draw the gating scheme
cyto_plot_gating_scheme(transformed_accuri[[3]],
                        back_gate = TRUE,
                        gate_track = TRUE)


cyto_stats_compute(transformed_accuri[2],
                   alias = c("One copy GFP"),
                   stet - "median",
                   channels = "GFP")

