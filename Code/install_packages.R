# Installs every package the pipeline uses (run once). Tested with R 4.5.1.
pk <- c("readxl","writexl","readr","dplyr","tidyr","tibble","stringr","purrr","lubridate","ggplot2","scales","patchwork","sf","rnaturalearth","rnaturalearthdata",
        "survey","survival","mice","glmnet","ranger","FactoMineR","broom","nnet","magick")
miss <- setdiff(pk, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) try(install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev"))
cat("all packages available:", all(pk %in% rownames(installed.packages())), "\n")
