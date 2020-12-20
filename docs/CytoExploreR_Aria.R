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

##Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#install("cytolib", "flowCore", "flowWorkspace", "openCyto")

##Install CytoExploreR from GitHub:
#devtools::install_github("DillonHammill/CytoExploreRData")
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
#library(CytoExploreRData)

#Setup experiment using cyto_setup()
#This function will read .fcs files to a cutoset which are then added to a GatingSet
#in generating the cytoset a experiment details csv files is created.  Additional columns can be added by right clicking on a mac and inserting columns
#and a experiment markers csv file is created
aria_gating_set <- cyto_setup(path="/Users/david/Projects/data/flowdata/Aria", select="fcs")

#To interactively edit the experiment details
#cyto_details_edit(aria_gating_set)

#To tqke a look at the data
cyto_plot_explore(aria_gating_set,
                  channels_x = "FSC-A",
                  channels_y = c("FSC-A","GFP","mCherry-A")
)

cyto_plot_profile(aria_gating_set,
                  parent = "root",
                  channels = c("FSC-A","GFP","mCherry-A")
                  )

cyto_plot_gating_tree(aria_gating_set)

#Transform the data using a logicle transformation.
#The first step applies the transformation using cyto_transformer_logicle()
#By default, this function does not apply the transformation to FSC and SSC, so we do that in a second step.
#we then need to combine the transformation using cyto_transformer_combined()
#finally we apply the transformation to all the data using cyto_transform
aria_transformed <- cyto_transformer_logicle(aria_gating_set)
aria_FSC_SSC_transformed <- cyto_transformer_logicle(aria_gating_set,
                                              channels = c("FSC-A", "FSC-W", "SSC-A"))

aria_combined_transformed <- cyto_transformer_combine(aria_transformed,
                                  aria_FSC_SSC_transformed)

transformed_aria <- cyto_transform(aria_gating_set,
                                    trans = aria_combined_transformed)


#To tqke a look at the transformed data
cyto_plot_explore(transformed_aria,
                  density_modal = TRUE,
                  axes_limits = "data",
              #    display = 100000,
                  point_col_alpha = 0.5,
                   point_col = "black",
                  channels_x = "FSC-A",
                  channels_y = c("GFP")
)


cyto_plot_profile(transformed_aria,
                  parent = "root",
                  channels = c("FSC-A","GFP","mCherry-A"), #add as many channels as needed.
                  legend = "fill"
)

cyto_plot_profile(transformed_aria,
                  parent = "root"
)

#####Gating cells.
#To gate cells we use the interactive function of Cytoexplorer
#The details of the gating are recorded in a .cvs file, which must be specified in the function
#The gating is done in a hierarchical manner, so that there is a parent and a child for each gate.

#Cytoexplorer will merge the entire set of fcs files and plot them when gating is performed.  As a result it is not possible to distinguish the different samples.
#However, there are two ways that individual samples can be used for the purpose of gating as described below.

#Below I describe three different ways to gate the data using the gating set.

#1.  Gating using the entire set of experimental data.
#This is the default behavior when the gating set is called.

#The first gate defines cells based on forward scatter and side scatter
#For this gate we use all the cells
cyto_gate_draw(transformed_aria,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
              )

#the next gate defines the singlets based on forward scatter height and width
cyto_gate_draw(transformed_aria,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-W"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv"
               )

#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "Negative",
               channels = c("GFP","mCherry-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv"
               )

#Define the one copy GFP gate
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "One_copy",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#Define the two copy GFP gate
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "Two_copy_gfp",
               channels = c("FSC-A","GFP"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#Define the one copy mCherry gate
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "One_copy_mcherry",
               channels = c("FSC-A","mCherry-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#Define the two copy GFP gate
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "Two_copy_mcherry",
               channels = c("FSC-A","mCherry-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)


#Define the one copy GFP one copy mCherry gate
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "multi_copy",
               channels = c("GFP","mCherry-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#To visualize the effect of the gating on each sample
cyto_plot_gating_scheme(transformed_aria[6])

#To visualize the gating tree
cyto_plot_gating_tree(transformed_aria)

#To visualize the freq of cells in each gate for a single sample
cyto_plot_gating_tree(transformed_aria[[1]], stat="freq")

stats_freq <- cyto_stats_compute(transformed_aria,
                             parent = "Single_cells",
                  #           alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                             alias = "",
                             stat="freq"
                           #  save_as = "stats_freq1.csv")
View(stats_freq)

stats_count1 <- cyto_stats_compute(transformed_aria,
                                 parent = "Single_cells",
                                 alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                 stat="count",
                                 save_as = "stats_count1.csv")

View(stats_freq1)
View(stats_count1)


#####################
#2. Gating strategy uses individual samples
#In this case we visualize individual samples, which are used to define the gates

#Remove all gates to start clean with no gates applied
cyto_gate_remove(transformed_aria,
                 gatingTemplate = "aria_gating.csv",
                 alias = "Cells")  #will remove Cells gate and all descendant gates

#First we gate using all the cells to define the cells and singlets as in the first strategy
cyto_gate_draw(transformed_aria,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#the next gate defines the singlets based on forward scatter height and width
cyto_gate_draw(transformed_aria,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv"
)

####Now we use individual samples to define the gates, which are applied to all samples in the gatingset

#Identify the negative control/no GFP sample and use it to gate non-fluorescent cells
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "Negative",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = "DGY1"),  #strain used to define no GFP signal
               gatingTemplate = "aria_gating.csv"
               )

#cyto_plot_gating_scheme(transformed_aria2[1:2], stat="freq")

##To overlay the negative cells when gating the single copy GFP, we need to extract the data from the relevant sample, which is the first sample in this case
negative <- cyto_extract(transformed_aria, "Single_cells")[[1]]


#Define the one copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "One_copy",
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = "DGY500"),  #strain used to define one copy of GFP
               gatingTemplate = "aria_gating.csv",
               overlay=negative  #will plot the negative cells as gray plots on the same plot
                )

#cyto_plot_gating_scheme(transformed_aria[2], stat="freq")

##To overlay the one GFP copy cells when gating the two copy GFP we need to extract the relevant data
one_copy <- cyto_extract(transformed_aria, "Single_cells")[[2]]

#Define the two copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "Two_copy",
               channels = c("FSC-A","FL1-A"),
               select = list(Strain = "DGY1315"),  #control strain used to define two copies of GFP
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
               overlay=one_copy  #will plot the one copy cells as gray points on the same plot for reference
)

cyto_plot_gating_scheme(transformed_aria[3], stat="freq")

##To overlay the two GFP copy cells when gating the more than two copies we need to extract the data
two_copy <- cyto_extract(transformed_aria, "Single_cells")[[3]]

#Define the three copy GFP gate using the relevant control sample
cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = "multi_copy",
               channels = c("FSC-A","FL1-A"),
             #  select = list(name = "DGY2158"),  #if available use a known control
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
               overlay=two_copy  #will plot the two copy cells as gray points on the same plot for reference
              )

#cyto_plot_gating_scheme(transformed_aria2, stat="freq")

cyto_plot_gating_tree(transformed_aria[[7]],
                      stat="freq")

cyto_plot_gating_scheme(transformed_aria[7],
                        back_gate = TRUE,
                        gate_track = TRUE)

stats_freq2 <- cyto_stats_compute(transformed_aria,
                                  parent = "Single_cells",
                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq2.csv")

stats_count2 <- cyto_stats_compute(transformed_aria,
                                   parent = "Single_cells",
                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                   stat="count",
                                   save_as = "stats_count2.csv")

View(stats_freq2)
View(stats_count2)

#############Drawing multiple gates simultaneously.
#It is possible to draw multiple gates on the same plot by defining the samples, extracting the data and plotting them as overlays

#Remove all gates to start clean
cyto_gate_remove(transformed_aria,
                 gatingTemplate = "aria_gating.csv",
                 alias = "Cells")

######First we have to gate the cells and singlets using the entire dataset
cyto_gate_draw(transformed_aria,  #entire gating set is plotted for gating purposes.  It is downsamples to 25,000 events for
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv",
)

#the next gate defines the singlets based on forward scatter height and width
cyto_gate_draw(transformed_aria,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "aria_gating.csv"
)

##To overlay the negative cells we need to extract the data from the relevant sample, which is the first sample in this case
negative <- cyto_extract(transformed_aria, "Single_cells")[[1]] #DGY1
##To overlay the one GFP copy cells  we need to extract the relevant data
one_copy <- cyto_extract(transformed_aria, "Single_cells")[[2]] #DGY500
##To overlay the two GFP copy cells  we need to extract the data
two_copy <- cyto_extract(transformed_aria, "Single_cells")[[3]] #DGY1315

cyto_gate_draw(transformed_aria,
               parent = "Single_cells",
               alias = c("Negative", "One_copy", "Two_copy","multi_copy"), #defines gate names (4 total in this case)
               channels = c("FSC-A","FL1-A"),
               axes_limits = "data",
               select = list(Strain = c("DGY1","DGY500","DGY1315","DGY2158")),  #control strains used to different copy numbers of GFP.  May not be necessary to do this.
               gatingTemplate = "aria_gating_three.csv",
               overlay = c(negative, one_copy, two_copy),  #the corresponding data for each control which has been extracted
               point_col = c("black", "green", "red", "blue")
               )

cyto_plot_gating_tree(transformed_aria[[7]],
                      stat="freq")

cyto_plot_gating_scheme(transformed_aria[7],
                        back_gate = TRUE,
                        gate_track = TRUE)

stats3 <- cyto_stats_compute(transformed_aria,
                             parent = "Single_cells",
                             alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                             stat="freq")

stats_freq3 <- cyto_stats_compute(transformed_aria,
                                  parent = "Single_cells",
                                  alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq3.csv")

stats_count3 <- cyto_stats_compute(transformed_aria,
                                   parent = "Single_cells",
                                   alias = c("Negative", "One_copy", "Two_copy", "multi_copy"),
                                   stat="count",
                                   save_as = "stats_count3.csv")

View(stats)
View(stats_freq3)
View(stats_count3)

#############Additional examples of using functions

#if you need to redraw a gate you need to use cyto_gate_edit()
cyto_gate_edit(transformed_aria,
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               gatingTemplate = "aria_test_data.csv",
               type = "boundary")

###Drawing a tree of the gated samples
cyto_plot_gating_tree(transformed_aria,
                      stat="freq")

#Draw the gating scheme for a single sample
cyto_plot_gating_scheme(transformed_aria[[3]],
                        back_gate = TRUE,
                        gate_track = TRUE)

cyto_stats_compute(transformed_aria,
                   alias = c("One copy GFP"),
                   stat = "median",
                   channels = "GFP")


#############Visualizing data
#Function for visualizing the data

cyto_plot_explore(transformed_aria,
                  channels_x = "FSC-A",
                  channels_y = "GFP"
)


cyto_plot(transformed_aria,
          parent="Cells",
          channels = "FSC-A",
     #    alias = "",
          xlim = c(10000, 2000000),
          density_fill = rep("blue",19),
          density_stack = 0.7)

cyto_plot(transformed_aria[1:19],
          parent="Cells",
          channels = "GFP",
        #  label_text = "Strain",
        legend = "line",
        xlim = c(10000, 2000000),
          density_fill = rep("green",19),
          density_stack = 0.7)

cyto_plot(transformed_aria[2],
          parent="Cells",
          alias = "",
          channels = c("FSC-A","GFP"),
          xlim = c(100000, 3000000),
          ylim = c(10000, 3000000),
      )

