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
#1
install("cytolib", "flowCore", "flowWorkspace", "openCyto")

#Then you can install CytoExploreR from GitHub:

devtools::install_github("DillonHammill/CytoExploreRData")
devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
#library(CytoExploreRData)

#Setup experiment using cyto_setup()
#This function will read .fcs files to a cutoset which are then added to a GatingSet
#in generating the cytoset a experiment details csv files is created.  Additional columns can be added by right clicking on a mac and inserting columns
#and a experiment markers csv file is created
accuri_gating_set <- cyto_setup(path="/Users/david/Projects/data/flowdata/Accuri", select="fcs")

#To interactively edit the experiment details
cyto_details_edit(accuri_gating_set)


cyto_plot_explore(accuri_gating_set,
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)


cyto_plot_profile(accuri_gating_set,
                  parent = "root",
                  channels = c("FSC-A","GFP")
                  )

cyto_plot_gating_tree(accuri_gating_set)

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
#For this gate we use all the cells
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

cyto_plot_gating_scheme(transformed_accuri)
cyto_plot_gating_tree(transformed_accuri, stat="freq")
cyto_plot_gating_scheme(transformed_accuri[3], stat="freq")
#####################

#To generate the gates based on individuals samples, we need to use the appropriate individual samples

#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = "No GFP",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = "DGY1"),  #strain used to define no GFP signal
               gatingTemplate = "Accuri_gating.csv"
               )

cyto_plot_gating_scheme(transformed_accuri[1:2], stat="freq")


##To overlay the No GFP cells when gating the single copy GFP
negative <- cyto_extract(transformed_accuri, "Single_cells")[[1]]


#Define the one copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = "One copy GFP",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = "DGY500"),  #strain used to define one copy of GFP
               gatingTemplate = "Accuri_gating.csv",
               overlay=negative  #will plot the negative cells as gray plots on the same plot
)

cyto_plot_gating_scheme(transformed_accuri[1:3], stat="freq")

##To overlay the one GFP copy cells when gating the two copy GFP
one_copy <- cyto_extract(transformed_accuri, "Single_cells")[[2]]

#Define the two copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = "Two copy GFP",
               channels = c("FSC-A","FL1-A"),
               select = list(Strain = "DGY1315"),  #strain used to define two copies of GFP
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv",
               overlay=one_copy  #will plot the one copy cells as gray points on the same plot for reference
)

cyto_plot_gating_scheme(transformed_accuri[4], stat="freq")

##To overlay the one GFP copy cells when gating the two copy GFP
two_copy <- cyto_extract(transformed_accuri, "Single_cells")[[3]]

#Define the three copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = "More than two copies GFP",
               channels = c("FSC-A","FL1-A"),
             #  select = list(name = "DGY2158"),
               axes_limits = "data",
               gatingTemplate = "Accuri_gating.csv",
               overlay=two_copy  #will plot the two copy cells as gray points on the same plot for reference
)


cyto_plot_gating_scheme(transformed_accuri, stat="freq")

cyto_plot_gating_tree(transformed_accuri[[7]],
                      stat="freq")

cyto_plot_gating_scheme(transformed_accuri[1],
                        back_gate = TRUE,
                        gate_track = TRUE)

stats <- cyto_stats_compute(transformed_accuri,
                   parent = "Single_cells",
                   alias = c("No GFP", "One copy GFP", "Two copy GFP", "More than two copies GFP"),
                   stat="freq")

View(stats)

#############Drawing multiple gates simultaneously.

#It is possible to draw multiple gates by defining them


cyto_gate_draw(transformed_accuri,
               parent = "Single_cells",
               alias = c("Neg", "One", "Two","More"), #defines gate names (4 total in this case)
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = c("DGY1","DGY500","DGY1315","DGY2158")),  #control strains used to different copy numbers of GFP
               gatingTemplate = "Accuri_gating.csv",
               overlay = c(negative, one_copy, two_copy),  #the corresponding data for each control which has been extracted
               point_col = c("black", "green", "red", "blue")
               )


#############Additional

#if you need to redraw a gate you need to use cyto_gate_edit()
cyto_gate_edit(transformed_accuri,
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "Accuri_test_data.csv",
               type = "boundary")

###Drawing a tree of the gated samples
cyto_plot_gating_tree(transformed_accuri,
                      stat="freq")


#Draw the gating scheme for a single sample
cyto_plot_gating_scheme(transformed_accuri[[3]],
                        back_gate = TRUE,
                        gate_track = TRUE)


cyto_stats_compute(transformed_accuri,
                   alias = c("One copy GFP"),
                   stat = "median",
                   channels = "GFP")


#############Visualizing data


cyto_plot_explore(transformed_accuri,
                  channels_x = "FSC-A",
                  channels_y = "GFP"
)


cyto_plot(transformed_accuri,
          parent="Cells",
          channels = "FSC-A",
          alias = "",
          xlim = c(10000, 2000000),
          density_fill = rep("blue",19),
          density_stack = 0.7)

cyto_plot(transformed_accuri[1:19],
          parent="Cells",
          channels = "GFP",
        #  alias = c("One copy GFP"),
          alias = "",
          xlim = c(10000, 2000000),
          density_fill = rep("green",19),
          density_stack = 0.7)

cyto_plot(transformed_accuri[2],
          parent="Cells",
          alias = c("Neg", "One", "Two", "More"),
          channels = c("FSC-A","GFP"),
          xlim = c(100000, 3000000),
          ylim = c(10000, 3000000),
      )

