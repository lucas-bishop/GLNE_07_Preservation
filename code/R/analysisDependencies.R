# Checking working directory ----------------------------------------------

# Checking to make sure the working directory is set to the correct location for downstream analyses

# If using RStudio, double-click the Rproj file to start a new session with the proper working directory
# OR
# If using the command line, verify location and set manually using setwd() if needed
if (!any(grepl("\\.Rproj", list.files()))) {
  stop("Please set the working directory to the parent folder containing the R project file (.Rproj).")
}



# Checking package installations ------------------------------------------

# Checking to see if all required packages are installed and letting user know if any are misisng

# List of packages used in this analysis
analysis_packages <- c("tidyverse", "RColorBrewer", "ggforce")

# Checks to see which required packages are missing
packages_to_install <- analysis_packages[!(analysis_packages %in% installed.packages()[,"Package"])]

# Errors out script and tells user which packages are missing
if (length(packages_to_install) > 0) {
  stop(paste("Please install the following packages required for downstream analyses:",
             paste(packages_to_install, collapse = ", ")))
}


