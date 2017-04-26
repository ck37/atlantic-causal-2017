source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F, verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

# Assumes 2017 data is extracted to inbound/pre_data
# Download file from ACIC 2017 website.
conf = list(
  inbound_dir = "inbound/pre_data",
  # Subdirectory to save temporary csv files.
  data_dir = "data",
  # Subdirectory to save output logs.
  output_dir = "output"
)

if (!dir.exists(conf$inbound_dir)) {
  stop(paste("Make sure to download 2017 data to:", conf$inbound_dir))
}

# Import X data.
data_x = read.csv(paste0(conf$inbound_dir, "/x.csv"), header = F)
str(data_x)

# Import X_Z data.
data_x_z = read.csv(paste0(conf$inbound_dir, "/X_subset_z.csv"), header = F)
str(data_x_z)

# Import X_Y data.
data_x_y = read.csv(paste0(conf$inbound_dir, "/X_subset_y.csv"), header = F)
str(data_x_y)

# Loop over 64 outcome / treatment files and combine into separate input datasets.
# Z files are treatment indicator.
# Y files are the outcome variable (continuous).
file_names = list.files(conf$inbound_dir,
                        pattern = "\\.[zy]\\.csv$",
                        full.names = F)
file_names

# Import files and save in a list.
files = list()
for (file_name in file_names) {
  files[[file_name]] = list(
    data = read.csv(paste0(conf$inbound_dir, "/", file_name))
  )
}

names(files)

# Save current results
save(files, file_names, data_x, data_x_y, data_x_z,
     file = paste0(conf$data_dir, "/import-2017.RData"))

