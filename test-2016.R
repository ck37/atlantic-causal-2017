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

analysis = list()

# Run algorithm on each file.
for (filename in names(files)) {
  cat("Analyzing", filename, "file", which(filename == names(files)), "of",
      length(files), "\n")

  file = files[[filename]]

  csv_filename = paste0(conf$data_dir, "/data-x", filename)

  # Save data to a temporary csv file for use in the later system2 call.
  # Include covariate data.
  data = cbind(file, data_x)
  write.csv(data, file = csv_filename)

  out1_filename = paste0(conf$data_dir, "/", filename, "-out1.csv")
  out2_filename = paste0(conf$data_dir, "/", filename, "-out2.csv")
  output_filename = paste0(conf$output_dir, "/", filename, ".out")

  # Run targeted_learner.R as a shell script to generate output.
  time_start = proc.time()

  # This will take 10+ minutes per dataset most likely.
  result = system2("./targeted_learning.R",
          args = c(csv_filename,
                   out1_filename,
                   out2_filename),
          stdout = output_filename, stderr = output_filename)
          # Enable this line below (and disable line above) to print output
          # directly to R console.
          #stdout = "", stderr = "")

  time_elapsed = proc.time() - time_start

  # result = 0 means that the command executed successfully.
  if (!result) {
    cat("Succeeded.\n")
  } else {
    cat("Failed with error code:", result, "\n")
  }

  # Save results
  results = list(
    data = data,
    time = time_elapsed,
    out1_filename = out1_filename,
    out2_filename = out2_filename,
    output_filename = output_filename
  )

  analysis[[filename]]  = results
}

save(analysis, file = paste0(conf$data_dir, "/test-2016.RData"))
