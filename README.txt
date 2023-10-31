This is a repository of files and R code to analyse the Huon Peninsula coral fossil data collected by John Pandolfi. This project examines temporal diversity trends in Holocene coral reefs.

Code used to generate results, figures and tables for DOI. Scripts makes extensive use of the code folding functionality in RStudio. Alt + O is the default shortcut on Windows and Linux versions of RStudio to collapse all folds.

full_analyses.R shows the complete data analysis pathway (from the output of data_acquisition_processing.R), as well as main and supplementary analyses. This script requires respository sub-folders etc as per the analysis script above.

####################
## Notable files: ##
####################

./huon_analyses.R - R script to process and analyse data and produce all figures and tables.

./raw.datafiles/huon_intercept_data.csv - this is the raw coral cliff intercept data, processed for analysis.
./raw.datafiles/huon_intercept_dataMETADATA.txt - metadata for raw coral cliff intercept data.

./raw.datafiles/huon_transect_data.csv - this is the transect-level data, processed for analysis.
./raw.datafiles/huon_transect_dataMETADATA.txt - metadata for transect-level data.