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

  # Can be "simple" or "complex".
  #sl_lib_type = "simple",
  sl_lib_type = "complex",

  # Set to T for extra output during execution.
  verbose = T,
  
  parallel = F,

  # Set auto-install to T for code to install any missing packages.
  auto_install = T,

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

if (conf$parallel) {
  options("mc.cores" = use_cores)
}

if (conf$verbose) {
  # Check how many parallel workers we are using:
  cat("Cores used:", getOption("mc.cores"), "\n")
}

if (conf$sl_lib_type == "simple") {
  q_lib = c("SL.mean", "SL.glmnet")
  g_lib = c("SL.mean", "SL.glmnet")
} else {
  q_lib = c(list(# speedglm doesn't work :/ just use plain ol' glm.
    c("SL.glm", "All", "screen.corRank4", "screen.corRank8", "prescreen.nosq")#,
    #c("SL.mgcv", "All", "prescreen.nosq"),
    #c("sg.gbm.2500", "prescreen.nocat"),
    #"SL.xgboost",
    #"SL.xgboost_threads_4"
    # Effect modification learners can't be used with g, only Q.
  ),
  # create.Learner() grids.
  sl_glmnet_em15$names,
  sl_xgb$names,
  sl_ksvm$names, 
  list(
    #"SL.randomForest_fast",
    "SL.ranger_fast",
    c("SL.glmnet_fast", "All", "screen.corRank4", "screen.corRank8"),
    c("SL.nnet", "All", "screen.corRank4", "screen.corRank8"),
    c("SL.earth", "prescreen.nosq"),
    # Works only if parallel = F. Do not use with mcSuperlearner!
    "SL.bartMachine",
    "SL.mean"))
  
  # Need a separate g lib that does not include effect modification learners.
  g_lib = c(list(c("SL.glm", "All", "screen.corRank4", "screen.corRank8", "prescreen.nosq"),
                 #c("SL.mgcv", "All", "prescreen.nosq"),
                 #c("sg.gbm.2500", "prescreen.nocat"),
                 #"SL.xgboost",
                 #"SL.xgboost_threads_4",
                 "SL.ranger_fast"#,
  ), # create.Learner() grids.
  sl_xgb$names,
  sl_ksvm$names, 
  list(
    c("SL.glmnet_fast", "All", "screen.corRank4", "screen.corRank8"),
    c("SL.nnet", "All", "screen.corRank4", "screen.corRank8"),
    c("SL.earth", "prescreen.nosq"),
    # Works only if parallel = F. Do not use with mcSuperlearner!
    "SL.bartMachine",
    "SL.mean"))
}

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
                       parallel = F,
                       verbose = conf$verbose)

# Parallel (4 cores): X seconds with stratified or pooled outcome regression.
# Non-parallel: X seconds with stratified, X seconds with pooled outcome regression.
results$time

# Extract out unit estimates to make displaying more convenient.
unit_estimates = results$unit_estimates
results$unit_estimates = NULL

# Display everything except the unit estimates.
results

summary(unit_estimates)

# TODO: test command-line version like in analyze-2016.R
