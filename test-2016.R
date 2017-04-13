# Dependencies: import-2016.R
source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F,
                  verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

# Assumes 2016 data is extracted to inbound/data-2016/
# Download from https://drive.google.com/file/d/0B8TUkApaUlsGekFSblJWa25NM1E/edit
conf = list(
  inbound_dir = "inbound/data-2016",
  # Subdirectory to save temporary csv files.
  data_dir = "data",
  # Subdirectory to save output logs.
  output_dir = "output"
)

input_file = paste0(conf$data_dir, "/import-2016.RData")

if (!file.exists(input_file)) {
  stop(paste("Can't find", input_file, ". Make sure to run import-2016.R first."))
}

load(input_file)

# First test on a single file.
names(files[[1]])

data = files[[1]]$data

q_lib = c("SL.mean", "SL.glmnet")
g_lib = c("SL.mean", "SL.glmnet")

set.seed(1, "L'Ecuyer-CMRG")

# Takes a minute or so to run.
results = estimate_att(A = data$z,
                       Y = data$y,
                       W = data_x,
                       SL.library = q_lib,
                       g.SL.library = g_lib,
                       verbose = T)

results