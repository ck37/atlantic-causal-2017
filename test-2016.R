source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F,
                  verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

# Assumes 2016 data is extracted to inbound/data-2016/
# Download from https://drive.google.com/file/d/0B8TUkApaUlsGekFSblJWa25NM1E/edit
inbound_dir = "inbound/data-2016"

# Import X data.
data_x = read.csv(paste0(inbound_dir, "/x.csv"))
str(data_x)

# Loop over 20 outcome / treatment files and combine into separate input datasets.
file_names = list.files(inbound_dir, pattern = "^zy_.+\\.csv", full.names = F)
file_names

# Import files and save in a list.
files = list()
for (file_name in file_names) {
  files[[file_name]] = read.csv(paste0(inbound_dir, "/", file_name))
}

names(files)

# First test on a single file.
names(files[[1]])

file = files[[1]]

q_lib = c("SL.mean", "SL.glmnet")
g_lib = c("SL.mean", "SL.glmnet")

set.seed(1, "L'Ecuyer-CMRG")

# Takes a minute or so to run.
results = estimate_att(A = file$z,
                       Y = file$y,
                       W = data_x,
                       SL.library = q_lib,
                       g.SL.library = g_lib,
                       verbose = T)

results

# TODO:  Run algorithm on each file.
for (file in files) {
}