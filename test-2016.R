# Dependencies: import-2016.R
source("lib/function_library.R")

# Assumes 2016 data is extracted to inbound/data-2016/
# Download from https://drive.google.com/file/d/0B8TUkApaUlsGekFSblJWa25NM1E/edit
conf = list(
  inbound_dir = "inbound/data-2016",

  # Subdirectory to save temporary csv files.
  data_dir = "data",

  # Subdirectory to save output logs.
  output_dir = "output",

  # Set to T for extra output during execution.
  verbose = T,

  # Set auto-install to T for code to install any missing packages.
  auto_install = F,

  # Use up to this many cores if available.
  max_cores = 4
)

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = conf$auto_install, verbose = conf$verbose)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = conf$verbose)

input_file = paste0(conf$data_dir, "/import-2016.RData")

if (!file.exists(input_file)) {
  stop(paste("Can't find", input_file, ". Make sure to run import-2016.R first."))
}

load(input_file)

# First test on a single file.
names(files[[1]])

data = files[[1]]$data

# Setup parallelization? Use up to 4 cores.
num_cores = RhpcBLASctl::get_num_cores()

if (conf$verbose) {
  cat("Cores detected:", num_cores, "\n")
}

use_cores = min(num_cores, conf$max_cores)
options("mc.cores" = use_cores)

if (conf$verbose) {
  # Check how many parallel workers we are using:
  cat("Cores used:", getOption("mc.cores"), "\n")
}

q_lib = c("SL.mean", "SL.glmnet")
g_lib = c("SL.mean", "SL.glmnet")

set.seed(1, "L'Ecuyer-CMRG")

# Restrict to 250 observations to be similar to contest.
obs = sample(nrow(data), 250, replace = F)
A = data$z[obs]
Y = data$y[obs]
W = data_x[obs, ]

# Takes a minute or so to run using a simple library.
results = estimate_att(A = A,
                       Y = Y,
                       W = W,
                       SL.library = q_lib,
                       g.SL.library = g_lib,
                       pooled_outcome = T,
                       parallel = T,
                       verbose = conf$verbose)

# Parallel (4 cores): 31 seconds with stratified or pooled outcome regression.
# Non-parallel: 76 seconds with stratified, 78 seconds with pooled outcome regression. 
results$time

# Extract out unit estimates to make displaying more convenient.
unit_estimates = results$unit_estimates
results$unit_estimates = NULL

# Display everything except the unit estimates.
results

summary(unit_estimates)

# TODO: test command-line version like in analyze-2016.R
