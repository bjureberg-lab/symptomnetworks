
# ------------------------------------------------------------------------------
# Loading functions from pruning_source.R (does not need to be changed)
# ------------------------------------------------------------------------------

source("pruning_source.R")

# ------------------------------------------------------------------------------
# Choosing networks to prune and plot, and how much pruning should be done
# ------------------------------------------------------------------------------

# Load data
##' Change "path" below to a string of your own path (or keep as is, if your folder
##' with data files is called "data" and is in the same directory as this R script.)
data_files <- list.files(path = "data/", pattern = "*.csv", full.names = TRUE)

# Create patient IDs based on file names
##' Does not need to be changed (used to get patient IDs/file names at top of page in pdf)
patientIDs <- str_extract(data_files, "(?<=data/).*(?=.csv)")

# Prune and plot
##' Change "save_as" below to a string of your own file name, and change "gamma"
##' to adjust level of pruning. Note that 0 <= "gamma" <= 1 where 0 yields no pruning,
##' 1 yields the same number of edges as the number of nodes. The default is
##' gamma = 0.95 as in the thesis (and this will affect both pruned networks,
##' to ensure that they have the same number of edges).
##'
##' Tip: If you want to produce a pdf for only the first patient in your folder,
##' you can change the first argument below to data_files[1]. If you want to do
##' this for the first, say, 3 patients in the folder, you can change the argument
##' to data_files[1:3] etc. The default is to use all files in the given folder.
pruned_to_pdf(data_files, save_as = "pruned.pdf", gamma = 0.95)
