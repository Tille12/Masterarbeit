# Masterarbeit


+++++++++++++++++++++++RUNNING THE VISUAL ANALYSIS++++++++++++++++++++++++++
=== INSTALL PACKAGES AND SOFTWARE ===
* install R version 3.6.3 or higher
* make sure to install the following packages:
* reagarding versions, "installedRpackages.txt" contains all packages with version, which were installed while the script was running
"ggplot2"
"rlang"
"network"
"igraph"
"BBmisc"
"viridisLite"
"viridis"
"RColorBrewer"
"scales"
"RCy3"
"methods"
"glue"
"sjmisc"
"tidyverse"
"cowplot"

* install Cytoscape Version 3.8.0 and open before running the R-Script
=== INSTALL PACKAGES AND SOFTWARE ===


=== INPUT PARAMETER ===
* open "plots.R" in an editor or Rstudio
* Within the script you can choose the input data (line 7)
* ... and whether or not Graphs should automativelly be saved in "./plots" (line 8)
* if the working directory (obtain it by getwd()) is not the "Masterarbeit" directory, change it accordingly (line 17)
* !!!SAVE AFTER MAKING YOUR CHOICE!!!
=== INPUT PARAMETER ===

=== RUNNING  "plots.R" ===
"plots.R" contains the R script to perform the analysis. It sources functions.R, which contains the major share of functions of the analysis.
It also contains some parts, which were not discusses in the thesis, and part of the exploratory work.
By Running the script line per line it is easier to follow along the analysis. But after providing the input parameter, it could also be run as entirely.
=!=!=!=! Start Cytoscape already before running the script !=!=!=!=!=

----------> run plots.R

=== RUNNING  "plots.R" ===
+++++++++++++++++++++++RUNNING THE VISUAL ANALYSIS++++++++++++++++++++++++++



