# ================ Description ================ 
# Visual Analysis as part of my master's thesis
# Author
# T. Rau
#
#### MAKE SURE TO OPEN CYTOSCAPE DESKTOP BEFORE RUNNING THE SCRIPT ####
#
# ================ Choose Running Parameters ================ 
# Choose, which dataset to run. 
# i) batch, set simulation to one out of "b", "batch"
# ii) chemostat, set simulation to one out of "c", "cstr", "chemostat"
simulation <- "c"
saveGraphs <- TRUE
# ================= SET WORKING DIRECTORY ================
# if working directory is not "./Masterarbeit" set it there with 
wDir = getwd()
# wDir <- setwd()
# ================ Libraries ================ 
library("ggplot2")
library("rlang")
library("network")
library("igraph")
library("BBmisc")
library("viridisLite")
library("viridis")
library("RColorBrewer")
library("scales")
library("RCy3")
library("methods")
library("RBGL")
library("glue")
library("sjmisc")
library("tidyverse")
library("cowplot")

# ================ Sources ================ 
#source("../functions.R")
source("functions.R")
# ================ Functions ================ 
# ================ Import ===================
# chemostat: SHIUMIx Chemostat: simulatedTrajectory_23-Nov-2019 17_22_42
# batch: SHIUMIx Batch aus dem paper: simulatedTrajectory_20-Dec-2019 11_57_36_Case_3_0.01_init_parallel
chemostat <- "simulatedTrajectory_23-Nov-2019 17_22_42"
batch <- "simulatedTrajectory_20-Dec-2019 11_57_36_Case_3_0.01_init_parallel"

## Set Import Directory/Folder
folderName = batch
if(!isBatch(simulation)){
  folderName = chemostat
}
importDir = paste(wDir, "data",folderName, sep = '/')
#importDir <- repath()
nodes <- read.csv(paste(importDir, "nodesList.txt", sep = '/'), header=T, as.is=T)
edges <- read.csv(paste(importDir, "edgeList.txt", sep = '/'), header=T, as.is=T)
trajectory <- read.csv(paste(importDir, "trajectory.txt", sep = '/'), header=T, as.is=T)

## remove all artificially created inflow/effluent nodes and edges
### first save them for later use (chemostat)
nodesIE <- nodes[which(nodes$type == "ie"), ]
edgesIE <- edges[which(edges$from == "ie") | which(edges$to == "ie"), ]
names(edgesIE) <- gsub("sense", "excrection", names(edgesIE), fixed = TRUE)
edgesIE$excretion <- rep(FALSE, nrow(edgesIE))
trajIE  <- trajectory[which(names(trajectory) %in% c("time", edgesIE$reacID))]  
### remove IE ...
modifiedDFs <- excludeNodesAndEdges(nodesDF = nodes, edgesDF = edges)
nodes <- modifiedDFs$nodes
edges <- modifiedDFs$edges
rmNodesIE <- modifiedDFs$rmNodeIds
rmEdgesIE <- modifiedDFs$rmEdgeReacIds

# ================ Set-Meta-Data ================ 
# FIXME: Should be automatically imported
batchTimeSlice <- c(0, 0.044, 0.090, 0.138, 0.250, 0.400, 0.750, 0.998)
batchPosSliceLabMoTraj <- c(0.1064, 0.106, 0.1064, 0.106, 0.1064, 0.1064, 0.1064, 0.1064)
batchPosSliceLabMoGrowth <- c(0.0435, 0.0405, 0.0435, 0.0405, 0.0435, 0.0435, 0.0435, 0.0435)
batchPosSliceLabCompTraj <- c(0.26, 0.28, 0.26, 0.28, 0.28, 0.28, 0.28, 0.28) # without h[e]
#batchPosSliceLabCompTraj <- c(0.39, 0.41, 0.39, 0.41, 0.41, 0.41, 0.41, 0.41) # with h[e]
batchPosSliceLabCompSelecGrowth <- c(3.5, 3.3, 3.5, 3.3, 3.5, 3.5, 3.5, 3.5) # wihtout h[e]
#batchPosSliceLabCompSelecGrowth <- c(5.3, 5.1, 5.3, 5.1, 5.3, 5.3, 5.3, 5.3) # with h[e]

chemTimeSlice <- c(0, 1.7247370, 3.7211585, 23.37732, 125.0131, 250.1034, 399.9201)
chemPosSliceLabMoTraj <- c(0.09, 0.10, 0.11, 0.12, 0.13, 0.13, 0.13) 
#chemPosSliceLabMoTraj <- c(c(0.105, 0.095, rep(0.105, 5))) # with xlim=30
chemPosSliceLabMoGrowth <- c(rep(0.027, 7))
chemPosSliceLabCompTraj <- c(rep(5.5,7))
chemPosSliceLabCompSelecGrowth <- c(rep(1.75,7))

batchMeta = list(title="SIHUMIx Batch Community", type="Batch", id = "batch", slice = batchTimeSlice, 
                 posSlLabMoTraj = batchPosSliceLabMoTraj, posSlLabMoGrowth = batchPosSliceLabMoGrowth,
                 posSlLabCompTraj = batchPosSliceLabCompTraj, posSlLabCompSelGrowth = batchPosSliceLabCompSelecGrowth)
chemMeta = list(title="SIHUMIx Chemostat Community", type="Chemostat", id = "chem", slice = chemTimeSlice,
                posSlLabMoTraj = chemPosSliceLabMoTraj, posSlLabMoGrowth = chemPosSliceLabMoGrowth,
                posSlLabCompTraj = chemPosSliceLabCompTraj, posSlLabCompSelGrowth = chemPosSliceLabCompSelecGrowth)

# choose right metaData according to input "simulation"
meta <- batchMeta
if(!isBatch(simulation)){
  meta <- chemMeta
}

# ================ Set consistent color palettes ================ 
# By naming the color vector, the ggplot2 function "scale_color_manual" will try to match names to colors
# --> CONSISTENCY over whole project
bNames <- nodes$name[which(nodes$type == "m")]
batchSortBNames <- bNames[c(7,5,2,1,4,8,6,3)]
batchSortingOrder <- c(7,5,2,1,4,8,6,3)
############TESTING COLORS
## Rainbow Palette: Default of ggplot2
colRain <- rainbow(8)
names(colRain) <- bNames
grid::grid.raster(colRain, interpolate = FALSE)
## Viridis Palette: Default "viridis" overall good performance
colVir <- viridis(8)
names(colVir) <- bNames
grid::grid.raster(colVir, interpolate = FALSE)
## RColorBrewer Palette: very versatile
colRBrewer <- brewer.pal(8, "Dark2")
names(colRBrewer) <- c(bNames[1:3], bNames[7], bNames[4:6], bNames[8])
grid::grid.raster(colRBrewer, interpolate = FALSE)
### colorbrewer2: from https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=8
### --> generate your own colors interactively
colCb <- c("#8dd3c7", "#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5")
grid::grid.raster(colCb, interpolate = FALSE)
names(colCb) <- c(bNames[1:3], bNames[7], bNames[4:6], bNames[8])

## Choose consistent color palette, which will then consistently used all over the script
## Defines which color palette will be used in ggplot's "scale_color_manual(value= *)"
## --> Choose one out of "rainbow", "viridis", "RColorBrewer", "CB2"
## c(colRain, colVir, colRBrewer, colCb)
moColor <- colRBrewer

# ================ Plotting ================ 
# ~==== 2.) trajectory ====
dirTrajectories <- paste(wDir, "plots", sep = '/')
## ~==== 2.i) biomass over time ====
trajBiomass <- selectTrajValuesByType(type = "m", traj = trajectory, edgesDf = edges)
biomassIDs <- names(trajBiomass)
biomassNames <- sapply(biomassIDs, getBiomassName, nodesDF = nodes)
names(trajBiomass) <- biomassNames
trajBiomass$time <- trajectory[,1]
trajBiomass <- trajBiomass[,c("time", biomassNames)]

# mit tibble dataframe
trajBiomass.t <- trajBiomass %>% pivot_longer(cols = colnames(trajBiomass)[-1], 
                                              names_to = "Microorganisms")
# Manually assign colors by using ColorPalettes from above
### trajectory
filename = paste0("3_2_", meta$id,"_mo_trajectory_sliced.jpeg")
gBiomass <- ggplot(trajBiomass.t, aes(x = time, 
                                      y = value, 
                                      col = Microorganisms, 
                                      group = Microorganisms)) + geom_line() +
  geom_vline(xintercept = meta$slice, linetype="dashed", 
             color = "black", size=0.5) +
  # y=meta$posSlLabMoTraj
  # meta$slice[c(1, 4:7)]
  annotate("label", x = meta$slice, y = meta$posSlLabMoTraj, label=round(meta$slice, 2)
           , size = 3.6, fill = "white", label.size = NA) + 
  # scale_y_log10() +
  # scale_x_log10() +
  # adding breaks = c("Anaero...4662", "Lac...DM1", "Esc...655") will reorder labels
  scale_color_manual(values = moColor, name="" ) +
  labs(title = meta$title, subtitle = "simulated trajectory of microorganisms") +
  theme_bw() +
  xlab("Time [h]") +
  ylab("Biomass Concentration [gDW/L]")+
  #xlim(NA, 50) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    #plot.caption = element_text(lineheight = 0.7),
    #legend.title = element_text(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size =14),
    axis.text = element_text(size =11)
  )
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories,
                      width = unit(10, "inches"), height = unit(7.07, "inches"))}

### growth rate
moDeriv <- derivTraj(trajBiomass)
moDerivRel <- moDeriv$rel
moDerivRel.t <- moDerivRel %>% pivot_longer(cols = colnames(trajBiomass)[-1], 
                                            names_to = "Microorganisms")
filename <- paste0("3_2_", meta$id,"_mo_growth.jpeg")
gBiomassGrowth <- ggplot(moDerivRel.t, aes(x = time, 
                                          y = value, 
                                          col = Microorganisms, 
                                          group = Microorganisms)) + geom_line() +
  geom_vline(xintercept = meta$slice, linetype="dashed", 
             color = "black", size=0.5) +
  # y=meta$posSlLabMoTraj
  # meta$slice[c(1, 4:7)]
  # y = c(c(0.026,0.0275,0.026,0.0275), meta$posSlLabMoGrowth[5:7])
  annotate("label", x = meta$slice, y = meta$posSlLabMoGrowth, label=round(meta$slice, 2)
           , size = 3.6, fill = "white", label.size = NA) + 
  #scale_y_log10() +
  scale_color_manual(values = moColor, name="") + 
  labs(title = meta$title, subtitle = "growth rate of microorganisms") +
  theme_bw() +
  xlab("Time [h]") +
  ylab("Biomass Growth [gDW/(L*h)]")+
  #xlim(NA, 30) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    #plot.caption = element_text(lineheight = 0.7),
    #legend.title = element_text(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size =14),
    axis.text = element_text(size =11)
  )
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories,
                      width = unit(10, "inches"), height = unit(7.07, "inches"))}

### Combine graphs gBiomass & gBiomassGrowth
prow <- plot_grid(gBiomass + theme(legend.position = "none"),
                  gBiomassGrowth + theme(legend.position = "none"),
                  align = 'vh',
                  labels = c("A", "B"),
                  hjust = -1,
                  nrow = 1)
legend_b <- get_legend(gBiomass + theme(legend.position="bottom",
                                        #legend.box.background = element_rect(color = "black"),
                                        legend.spacing.x = unit(0.2, "cm") ))
filename <- paste0("3_2_", meta$id,"_mo_traj&growth.jpeg")
gBiomassComb <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .2))
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories, 
                      width = unit(14.1, "inches"), height = unit(7.07, "inches"))}

## ~==== 2.ii) compounds over time ====
trajComp <- selectTrajValuesByType(type = "c", traj = trajectory, edgesDf = edges)
compShortNames <- sapply(names(trajComp), getShortName, nodesDF=nodes)
# compLongNames <- sapply(compShortNames, getLongName) # ergibt keinen Sinn hier longNames zu benutzen
names(trajComp) <- compShortNames
trajComp$time <- trajectory[,1]
trajComp <- trajComp[,c("time", compShortNames)]
trajComp.t <- trajComp %>% pivot_longer(cols = colnames(trajComp)[-1], 
                                        names_to = "Compounds")
# Automatically assign colors by using "scale_color_*" functions
filename = paste0("3_2_", meta$id,"_comp_trajectory.jpeg")
ggplot(trajComp.t, aes(x = time, 
                       y = value, 
                       col = Compounds, 
                       group = Compounds)) + geom_line() +
  # scale_y_log10() +
  scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "simulated trajectory of compounds") +
  xlab("Time [h]") +
  ylab("Compound Concentration [mmol/L]") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "none",
    legend.text = element_text(size = 12),
    axis.title = element_text(size =14),
    axis.text = element_text(size =11)
  )
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories)}


## ~==== 2.iii) Explore relevance of compounds ====
### Get compounds with highest relevance, hence highest overall sum
sumCompConc <- sumTrajValuesByType(traj = trajectory, type = "c", edgesDF = edges)
names(sumCompConc) <- sapply(names(sumCompConc), function(x, nodesDF = nodes){
  nodesDF$name[which(nodesDF$id == x)]})
slCompAbund <- sumCompConc[1:9]

### Play with derivations of compound concentration
compDerivations <- derivTraj(traj = trajComp)
compDerAbs <- compDerivations$abs
compDerRel <- compDerivations$rel
#### all absolute changes are turned positive and then summed up
#### in order go get a sense of compound importance
sumCompDerAbs <- apply(compDerAbs, 2, function(col){sum(abs(col))}) 
sumCompDerAbs <- sumCompDerAbs[order(sumCompDerAbs, decreasing = TRUE)]
slCompChanges <- sumCompDerAbs[1:9]
slCompMVP <- union(names(slCompAbund), names(slCompChanges))
comp2omit <- c("h2o[e]", "time", "h[e]")  # further candidates: "h[e]"
slCompMVP <- omitComps(slCompMVP, comp2omit)

#### ***MVP**** Just the compounds trajectory from MVP shortlist
trajCompMVP.t <- trajComp.t %>% filter(Compounds %in% slCompMVP)
filename = paste0("3_2_", meta$id,"_comp_trajectory_MVP.jpeg")
gCompTraj <- ggplot(trajCompMVP.t, aes(x = time, 
                                       y = value, 
                                       col = Compounds, 
                                       group = Compounds)) + geom_line() +
  geom_vline(xintercept = meta$slice, linetype="dashed", 
             color = "black", size=0.5) +
  annotate("label", x = meta$slice, y=meta$posSlLabCompTraj, label=round(meta$slice, 2), 
           size = 3.4, fill = "white", label.size = NA) + 
  scale_y_log10() +
  #scale_color_viridis(discrete = TRUE, name = "") + 
  labs(title = meta$title, subtitle = "simulated trajectory of selected compounds") +
  theme_bw() +
  xlab("Time [h]") +
  ylab("Compound Concentration [mmol/L]") +
  #xlim(NA, 30) +
  #scale_x_log10() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "none",
    legend.text = element_text(size = 12),
    axis.title = element_text(size =14),
    axis.text = element_text(size =11)
  )
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories,
                      width = unit(10, "inches"), height = unit(7.07, "inches"))}






#### ABSOLUTE Compound consumption / production over time
# delta(c)_(n-1) = c(t_n) - c(t_(n-1))
compDerAbs.t <- compDerAbs %>% pivot_longer(cols = colnames(compDerAbs)[-1], 
                                            names_to = "Compounds")
compDerAbs.t %>% filter(value >=1) # nur Acetat, Formate und Protone steigen am Ende sprunghaftig an
# Automatically assign colors by using "scale_color_*" functions
filename = paste0("3_2_", meta$id,"_comp_derivative_abs.jpeg")
ggplot(compDerAbs.t, aes(x = time, y = value, col = Compounds, 
                         group = Compounds)) + geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "Absolute compound consumption / production over time") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "none"
  ) +
  xlab("Time [h]") +
  ylab("Compound concentration change between steps [mmol/L]")
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories)}

#### RELATIVE Compound consumption / production over time
# delta(c) = (c(k) - c(k-1)) / (t - t_t-1)
compDerRel.t <- compDerRel %>% pivot_longer(cols = colnames(compDerRel)[-1], 
                                            names_to = "Compounds")
compDerRel.t %>% filter(value >=1) # nur Acetat, Formate und Protone steigen am Ende sprunghaftig an
# Automatically assign colors by using "scale_color_*" functions
filename = paste0("3_2_", meta$id,"_comp_derivative_rel.jpeg")
ggplot(compDerRel.t, aes(x = time, y = value, col = Compounds, 
                         group = Compounds)) + geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "Relative compound consumption / production over time") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "none"
  ) +
  xlab("Time [h]") +
  ylab("Compound Concentration change [mmol/(L*h)]")
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories)}

#### ***MVP*** RELATIVE Compound consumption / production over time
compDerRelMVP.t <- compDerRel.t %>% filter(Compounds %in% slCompMVP)
filename = paste0("3_2_", meta$id,"_comp_derivative_rel_MVP.jpeg")
gCompDerRel <- ggplot(compDerRelMVP.t, aes(x = time, 
                                           y = value, 
                                           col = Compounds, 
                                           group = Compounds)) + geom_line() +
  geom_vline(xintercept = meta$slice, linetype="dashed", 
             color = "black", size=0.5) +
  annotate("label", x = meta$slice, y=meta$posSlLabCompSelGrowth, label=round(meta$slice, 2), 
           size = 3.4, fill = "white", label.size = NA) + 
  scale_x_log10() +
  #scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "relative consumption / production of selected compounds") +
  theme_bw() +
  xlab("Time [h]") +
  ylab("Compound Concentration change [mmol/(L*h)]")+
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "none",
    legend.text = element_text(size = 12),
    axis.title = element_text(size =14),
    axis.text = element_text(size =11)
  )
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories,
                      width = unit(10, "inches"), height = unit(7.07, "inches"))}

### Combine graphs gCompTraj & gCompDerRel
prow <- plot_grid(gCompTraj + theme(legend.position = "none"),
                  gCompDerRel + theme(legend.position = "none"),
                  align = 'vh',
                  labels = c("C", "D"),
                  hjust = -1,
                  nrow = 1)
legend_b <- get_legend(gCompTraj + theme(legend.position="bottom",
                                         #legend.box.background = element_rect(color = "black"),
                                         legend.spacing.x = unit(0.4, "cm")) +
                         guides(col = guide_legend(nrow=2, byrow=TRUE)))
filename <- paste0("3_2_", meta$id,"_comp_traj&change.jpeg")
gCompComb <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .2))
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories, 
                      width = unit(12.9, "inches"), height = unit(7.07, "inches"))}

##### Sum of relative consump/prod over time
integRelDeriv <- integrateConsumpRates(compDerRel)
integRelDeriv.t <- integRelDeriv %>% pivot_longer(cols = colnames(integRelDeriv)[-1], 
                                                  names_to = "Change")
filename = paste0("3_2_", meta$id,"_comp_integratedDeriv_rel.jpeg")
ggplot(integRelDeriv.t, aes(x = time, y = value, col = Change, 
                            group = Change)) + geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "Total relative compound consumption / production over time", caption = "simulated by microbialSim") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    plot.caption = element_text(lineheight = 0.7),
    legend.position = "right"
  ) +
  xlab("Time [h]") +
  ylab("Compound Concentration change [mmol/L*h]")
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories)}

#### RELATIVE CHANGE in Compound consumption / production over time
# 
compSecondDerivation <- derivTraj(traj = compDerRel)
compSecondDerivationRel <- compSecondDerivation$rel
compSecondDerivationRel.t <- compSecondDerivationRel %>% pivot_longer(cols = colnames(compSecondDerivationRel)[-1], 
                                                                      names_to = "Compounds")
compSecondDerivationRel.t %>% filter(value >=1) # nur Acetat, Formate und Protone steigen am Ende sprunghaftig an
# Automatically assign colors by using "scale_color_*" functions
filename = paste0("3_2_", meta$id,"_comp_secondDerivative_rel.jpeg")
ggplot(compSecondDerivationRel.t, aes(x = time, y = value, col = Compounds, 
                                      group = Compounds)) + geom_line() +
  scale_color_viridis(discrete = TRUE) + 
  labs(title = meta$title, subtitle = "Change in Relative compound consumption / production over time", caption = "simulated by microbialSim") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    plot.caption = element_text(lineheight = 0.7),
    legend.position = "none"
  ) +
  xlab("Time [h]") +
  ylab("Compound change in production/consumption [mmol/L*h^2]")
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories)}


## ~==== 2.iiii) Flux analysis over time ====
### each flux in trajectory data needs to be scanned for uptake "moments"
### iff there are some, inverse edges will be created in edges and traj and
### values be copied along
inversedFlows <- getInverseFlowPeriods(traj = trajectory)
resultInversedFlowsIncorp <- incorpInverseEdges(traj = trajectory, invFp =  inversedFlows)
modEdges <- resultInversedFlowsIncorp$edges
modTraj <- resultInversedFlowsIncorp$trajectory

itrajFlux <- selectTrajValuesByType("r", traj = modTraj, edgesDf = modEdges, nodesDf = nodes)
trajFluxInteg <- integrateFluxes(traj = modTraj, edgeDF = modEdges)
trajFluxInteg.t <- trajFluxInteg %>% pivot_longer(cols = colnames(trajFluxInteg)[-1], 
                                                  names_to = "statistics")
# Manually assign colors by using ColorPalettes from above
### trajectory
filename = paste0("3_2_", meta$id,"_flux_trajectory_integrated_stats_log.jpeg")
ggplot(trajFluxInteg.t, aes(x = time, 
                            y = value, 
                            col = statistics, 
                            group = statistics)) + geom_line() +
  geom_vline(xintercept = meta$slice, linetype="dashed", 
             color = "black", size=0.5) +
  annotate("text", x = meta$slice, y=rep(max(trajFluxInteg.t$value), 8), label=round(meta$slice, 2), size = 2) + 
  scale_y_log10() +
  #scale_color_viridis(discrete = TRUE) + 
  scale_color_manual(values = colExcretion, name="") + 
  labs(title = meta$title, subtitle = "Excretion/Uptake flux statistics over time", caption = "simulated by microbialSim") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    plot.caption = element_text(lineheight = 0.7)
  ) +
  theme_bw() +
  xlab("Time [h]") +
  ylab("flux [mmol/gDW*h]")
if(saveGraphs){ggsave(filename = filename, path = dirTrajectories,
                      width = unit(10, "inches"), height = unit(7.07, "inches"))}



## ~==== 3.) Full visualization of metabolic network ====
### FIXME: Might not be useful, because in this representation,
### edges in both directions are included, which cannot be there together in real
### preparation to load into cytoscape
edgesCyto <- replEdgeNamesOfEdgesDF(modEdges) # change header of edges
cytoscapePing()
cytoscapeVersionInfo()
copyVisualStyle('Directed', 'Custom')
## ~==== 3.i) Full visualization of metabolic network ====
fullCytoNet <- createNetworkFromDataFrames(nodes, edgesCyto, 
                                           title="Full model without temporal data", 
                                           collection=meta$title)
applyCustomCytoStyle(edgesCyto)



## ~==== 4) Full visualization of metabolic network at different timepoints ====
### FIXME: Include optional data selection paragraph
### Data Selection
selection <- excludeNodesAndEdges(nodes, modEdges, nodeIDs = c("c_62", "c_64")) # water and proton

NonCarbonCompounds <- c("", "")
selectionCarbon <- excludeNodesAndEdges(nodes, edges, nodeIDs = NonCarbonCompounds)

### Slicing the trajectory and create list with edge and node sets at different timepoints
iSlices <- extractFullSlices(traj = modTraj, edgesDF = modEdges, nodesDF = nodes, 
                             timeSlices = meta$slice, output = "i", addMue = TRUE, rmZero = TRUE)
iSlicesSelect <- extractFullSlices(traj = modTraj, edgesDF = selection$edges, nodesDF = selection$nodes, 
                                   timeSlices = meta$slice, output = "i", addMue = TRUE, rmZero = TRUE)


### violin plot of fluxes at timeSLices
nRowAllFlux <- 0
for(idx in 1:length(iSlices)){
  nRowAllFlux <- nRowAllFlux + nrow(iSlices[[idx]][["edges"]])
}
trajSlicedFlux <- data.frame(timepoints = vector(mode = "numeric", length = nRowAllFlux),
                             #reacID = vector(mode = "character", length = nRowAllFlux),
                             excretion = vector(mode = "logical", length = nRowAllFlux),
                             weight = vector(mode = "numeric", length = nRowAllFlux),
                             stringsAsFactors = FALSE)
start <- 1
newStart <- 1
for(idx in 1:length(iSlices)){
  start <- newStart
  currentEdges <- iSlices[[idx]][["edges"]]
  currentLength <- nrow(currentEdges)
  newStart <- start + currentLength
  currentTime <- names(iSlices)[idx]
  print(currentLength)
  print(paste0(start, " to ", newStart, ": ", newStart-start))
  trajSlicedFlux$timepoints[start:(newStart-1)] <- rep(currentTime, currentLength)
  #trajSlicedFlux$reacID[start:(newStart-1)] <- currentEdges$reacID
  trajSlicedFlux$excretion[start:(newStart-1)] <- currentEdges$excretion
  trajSlicedFlux$weight[start:(newStart-1)] <- currentEdges$weight
}
trajSlicedFlux$timepoints <- as.factor(trajSlicedFlux$timepoints)
trajSlicedFlux$excretion <- as.factor(trajSlicedFlux$excretion)
#trajSlicedFluxPos <- trajSlicedFlux %>% filter(flux>0)
ggplot(trajSlicedFlux, aes(x=timepoints, y=weight, fill=excretion)) +
  scale_y_log10() +
  geom_violin(trim=TRUE, scale="count", width=1, position=position_dodge(width=0.7)) +
  #geom_boxplot(trim=TRUE) +
  theme_classic()

for (idx in 1:length(iSlices)){
  timePoint <- names(iSlices)[idx]
  summary(trajSlicedFlux[which(trajSlicedFlux$timepoints == timePoint), -which(names(trajSlicedFlux) == "reacID")])
}

### Cumulative Flux at eacht time Slices
cumDfTrajSlicedFlux <- data.frame(timepoints = vector(mode = "numeric", length = length(meta$slice)),
                                  all = vector(mode = "logical", length = length(meta$slice)),
                                  excretion = vector(mode = "logical", length = length(meta$slice)),
                                  uptake = vector(mode = "logical", length = length(meta$slice)),
                                  stringsAsFactors = FALSE)

cumTrajSlicedFlux <- lapply(vector(mode = 'list', length(meta$slice)), function(x){
  x <- vector(mode='list', 1)})
names(cumTrajSlicedFlux) <- meta$slice

for(idx in 1:length(meta$slice)){
  currentTraj <- trajSlicedFlux %>% filter(timepoints == meta$slice[idx])
  currentFlux <- currentTraj %>% .$weight
  currentExcretionFlux <- currentTraj %>% filter(excretion == TRUE) %>% .$weight
  currentUptakeFlux <- currentTraj %>% filter(excretion == FALSE) %>% .$weight
  
  sumExcretion <- sum(currentExcretionFlux)
  sumUptake <- sum(currentUptakeFlux)
  sumAll <- sumExcretion + sumUptake
  
  dezilEx <- quantile(currentExcretionFlux, seq(0,10)/10)
  dezilUp <- quantile(currentExcretionFlux, seq(0,10)/10)
  dezilAll <- quantile(currentFlux, seq(0,10)/10)
  
  meanEx <- mean(currentExcretionFlux)
  meanUp <- mean(currentUptakeFlux)
  meanAll <- mean(currentFlux)
  
  medianEx <- dezilEx[6]
  medianUp <- dezilUp[6]
  medianAll <- dezilAll[6]
  
  exList <- list(sum=sumExcretion, mean=meanEx, median=medianEx, dezil=dezilEx)
  upList <- list(sum=sumUptake, mean=meanUp, median=medianUp, dezil=dezilUp)
  allList <- list(sum=sumAll, mean=meanAll, median=medianAll, dezil=dezilAll)
  cumTrajSlicedFlux[[idx]] <- list(excretion = exList, uptake = upList, all = allList)
  
  cumDfTrajSlicedFlux$timepoints[idx] <- meta$slice[idx]
  cumDfTrajSlicedFlux$all[idx] <- sumAll
  cumDfTrajSlicedFlux$excretion[idx] <- sumExcretion
  cumDfTrajSlicedFlux$uptake[idx] <- sumUptake
}



## ~==== 4.i) Network characterization with igraph ====
### Exemplary take Network at t = 0
### List of all parameters which will be extracted/calc for every slice
### FIXME: Needs to be updated, because initialization of iSliceStats refers to the length of this vector

### header of nodes needs to be changed in order that igraph recognizes ids as names
iNodes <- nodes
names(iNodes) <- replNodeNamesToI(names(iNodes))

### Calculate Stats for each FULL time slice
extractionParaFull <- c("flux", "excretion", "dia", "diaPath", "wDia", "wDiaPath", "pageRank", "wPageRank", 
                        "vertBetween", "edgeBetween", "degreeCentrIn", "degreeCentrOut", "degreeCentrTot")
### on full iSlices
iSlicesStats <- lapply(vector(mode = 'list', length(meta$slice)), function(x){
  x <- vector(mode='list', length(extractionParaFull))})
names(iSlicesStats) <- meta$slice 
iSlicesStats <- iGraphStats(iSlices, iSlicesStats, iNodes)

### on iSlicesSelect
iSlicesStatsSelect <- lapply(vector(mode = 'list', length(meta$slice)), function(x){
  x <- vector(mode='list', length(extractionParaFull))})
names(iSlicesStatsSelect) <- meta$slice 
iSlicesStatsSelect <- iGraphStats(iSlicesSelect, iSlicesStatsSelect, iNodes)
warnings <- warnings()


## ~==== 5) Linearization at given timepoints ====
start_time <- Sys.time()
iSlicesLin <- linearizeSlices(iSlices)
end_time <- Sys.time()
end_time - start_time
### Check frequency of removed edges
freqRmEdges <- freqRmEdges(iSlicesLin = iSlicesLin)

start_time <- Sys.time()
iSlicesLinSelect <- linearizeSlices(iSlicesSelect)
end_time <- Sys.time()
end_time - start_time
freqRmEdgesSelect <- freqRmEdges(iSlicesLin = iSlicesSelect)


### Calculate Stats for each linearized time slice
extractionParaLin <- c(extractionParaFull, c("maxFluxPath", "dist"))

iSlicesLinStats <- lapply(vector(mode = 'list', length(meta$slice)), function(x){
  x <- vector(mode='list', length(extractionParaLin))})
names(iSlicesLinStats) <- meta$slice 
iSlicesLinStats <- iGraphStats(iSlicesLin, iSlicesLinStats, iNodes, linearized = TRUE)

### on iSlicesLinSelect
iSlicesLinStatsSelect <- lapply(vector(mode = 'list', length(meta$slice)), function(x){
  x <- vector(mode='list', length(extractionParaLin))})
names(iSlicesLinStatsSelect) <- meta$slice 
iSlicesLinStatsSelect <- iGraphStats(iSlicesLinSelect, iSlicesLinStatsSelect, iNodes, linearized = TRUE)

### Calculate Comparison Stats
sliceComparison <- compareSlices(iSlicesStats, iSlicesLinStats)
sliceComparisonSelect <- compareSlices(iSlicesStatsSelect, iSlicesLinStatsSelect)



### Load into cytoScape
cytoscapePing()
cytoscapeVersionInfo()
copyVisualStyle('Directed', 'Custom')

#### FULL SLICES
slicedCytoNets <- vector(mode = "list", length(meta$slice))
names(slicedCytoNets) <- meta$slice
# length(meta$slice)
for (i in c(1,3,4,5,7)){
  cytoEdges <- iSlices[[i]]$edges
  cytoNodes <- iSlices[[i]]$nodes
  names(cytoEdges) <- replEdgeNamesToCyto(names(cytoEdges))
  
  # Round values in order to be displayed clearly in cytoscape
  cytoEdges$flux <- round(cytoEdges$flux, 3)
  cytoNodes$value <- round(cytoNodes$value, 3)
  cytoNodes$mue <- round(cytoNodes$mue, 3)
  
  slicedCytoNets[i] <- createNetworkFromDataFrames(cytoNodes, cytoEdges, 
                                                   title=meta$slice[i], 
                                                   collection=paste0(meta$title, ": FULL sliced model"))
  applyCustomCytoStyle(cytoEdges)
}


#### LIN SLICES
linSlicedCytoNets <- vector(mode = "list", length(meta$slice))
names(linSlicedCytoNets) <- meta$slice

for (i in 1:length(meta$slice)){
  cytoEdges <- iSlicesLinSelect[[i]]$edges
  cytoNodes <- iSlices[[i]]$nodes
  names(cytoEdges) <- replEdgeNamesToCyto(names(cytoEdges))
  
  # Round values in order to be displayed clearly in cytoscape
  cytoEdges$flux <- round(cytoEdges$flux, 3)
  cytoNodes$value <- round(cytoNodes$value, 3)
  cytoNodes$mue <- round(cytoNodes$mue, 3)
  
  linSlicedCytoNets[i] <- createNetworkFromDataFrames(cytoNodes, cytoEdges, 
                                                      title=meta$slice[i], 
                                                      collection=paste0(meta$slice[i], ": NEW Linearized sliced model"))
  #applyCustomCytoStyle(cytoEdges)
}



#############################################################################
#### SELECTEDSELECTEDSELECTEDSELECTEDSELECTEDSELECTEDSELECTEDSELECTEDSELECTED
#############################################################################
#### FULL SLICES
cytoscapePing()
cytoscapeVersionInfo()
copyVisualStyle('Directed', 'Custom')
slicedCytoNetsSelect <- vector(mode = "list", length(meta$slice))
names(slicedCytoNetsSelect) <- meta$slice

for (i in 1:length(meta$slice)){
  cytoEdges <- iSlicesSelect[[i]]$edges
  cytoNodes <- iSlicesSelect[[i]]$nodes
  names(cytoEdges) <- replEdgeNamesToCyto(names(cytoEdges))
  
  # Round values in order to be displayed clearly in cytoscape
  cytoEdges$flux <- round(cytoEdges$flux, 3)
  cytoNodes$value <- round(cytoNodes$value, 3)
  cytoNodes$mue <- round(cytoNodes$mue, 3)
  
  slicedCytoNetsSelect[i] <- createNetworkFromDataFrames(cytoNodes, cytoEdges, 
                                                         title=meta$slice[i], 
                                                         collection=paste0(meta$title, ": Full sliced selected model"))
  applyCustomCytoStyle(cytoEdges)
}


#### LIN SLICES
linSlicedCytoNetsSelect <- vector(mode = "list", length(meta$slice))
names(linSlicedCytoNetsSelect) <- meta$slice

for (i in 1:length(meta$slice)){
  cytoEdges <- iSlicesLinSelect[[i]]$edges
  cytoNodes <- iSlicesSelect[[i]]$nodes
  names(cytoEdges) <- replEdgeNamesToCyto(names(cytoEdges))
  
  # Round values in order to be displayed clearly in cytoscape
  cytoEdges$flux <- round(cytoEdges$flux, 3)
  cytoNodes$value <- round(cytoNodes$value, 3)
  cytoNodes$mue <- round(cytoNodes$mue, 3)
  
  linSlicedCytoNetsSelect[i] <- createNetworkFromDataFrames(cytoNodes, cytoEdges, 
                                                            title=meta$slice[i], 
                                                            collection=paste0(meta$title, ": Linearized sliced selected model"))
  applyCustomCytoStyle(cytoEdges)
}




