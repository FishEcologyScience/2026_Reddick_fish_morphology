## --------------------------------------------------------------#
## Script name: script00-00_user_interface.R
##
## Purpose of script:
##    A central location to order and source project scripts
##    Script naming:
##      - scriptXX-YY format (XX = class 00-99, YY = script 01-99)
##      - higher numbers depend on lower numbers
##      - letters (a,b,c) indicate no dependency between same-numbered scripts
##
## Author:
##
## Date Created:
##
## --------------------------------------------------------------#
## Modification Notes:
##
## --------------------------------------------------------------#


### Core Data Processing
#----------------------------#
source("02_scripts/script00-01_load_packages.R") # Load packages and helper functions first 
source("02_scripts/script01-01_import_format_singlefile.R") # Import and clean morphometric data 
source("02_scripts/script01-02_morphology_plots.R") # Run loop - per species to generate summary plots







