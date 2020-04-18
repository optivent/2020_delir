
library(tidyverse)
library(here)

data <- readRDS("input/data.rds")

alldata <- data$alldata
commondata <- data$commondata

write_rds(alldata, path = paste0(here("input"),"/alldata.rds"), compress = "none")
