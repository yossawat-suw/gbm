library(styler)
# Define the directory containing your R files
dir_path <- "." # nolint

# Get a list of all R files in the directory
r_files <- list.files(path = dir_path, pattern = "\\.(R|Rmd)$", full.names = TRUE)

# Apply styler to each file
lapply(r_files, style_file)
