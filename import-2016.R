source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F, verbose = T)

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

if (!dir.exists(conf$inbound_dir)) {
  stop(paste("Make sure to download 2016 data to:", conf$inbound_dir))
}

# Import X data.
data_x = read.csv(paste0(conf$inbound_dir, "/x.csv"))
str(data_x)

# Loop over 20 outcome / treatment files and combine into separate input datasets.
file_names = list.files(conf$inbound_dir, pattern = "^zy_.+\\.csv", full.names = F)
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
save(files, file_names, data_x,
     file = paste0(conf$data_dir, "/import-2016.RData"))

