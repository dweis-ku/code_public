# Converts HDExaminer exported "Peptide Pool Results" csv file into tidy format
# Saves tidy output as CSV and Rdata
# authored by David Weis, University of Kansas
# extensively modified to accommodate changes in parsing following readr 2.0
# release

library(tidyverse)
library(stringr)
library(janitor)

# load user parameters
source("scripts/user_parameters.R")
out_file_str <- str_remove(hx_data_file,"input/")
out_file_str <- str_c(str_sub(out_file_str, 1, -5),"_tidy")
out_file_csv <- str_c("output/", out_file_str, ".csv")
out_file_R <- str_c("output/", out_file_str, ".Rdata")

peptide_ID_file <- str_c(peptide_ID_file)
peptide_list <- read_csv(str_c(peptide_ID_file)) %>%
  select(pep_ID, pep_sequence)
hx_data_file <- str_c(hx_data_file)
# the first row contains HX labeling times and lots of white space
# capture the HX labeling times
hx_times <- read_csv(hx_data_file, n_max = 1, col_names = FALSE, 
                     col_types = cols(.default = "c"))
hx_times <- hx_times %>% unlist() %>% unname()
hx_times <- hx_times[!is.na(hx_times)]
# read the rest of the file dropping some useless columns
hx_data <- read_csv(hx_data_file,
                    skip = 1,
                    col_select = !contains(c("right", "protein", "_rt", 
                                             "search", "max", "percent")),                   
                    name_repair = make_clean_names
                    )
hx_data <- hx_data %>% rename_with(~ str_remove(.,"number_"))
# first instance of each variable does not have an index, this is annoying
hx_data <- hx_data %>% rename(d_1 = d)
if ("conf" %in% colnames(hx_data)) {
  hx_data <- hx_data %>% rename(conf_1 = conf)
}
if ("score" %in% colnames(hx_data)) {
  hx_data <- hx_data %>% rename(score_1 = score)
}

# DEV: this version assumes that the time unit is a single character at the end
# of the string 
# DEV: and that all times are in seconds

# determine if each variable is an HX time replicate of the preceding one
replicate_cnt <- vector("character", length(hx_times))
j <- 1
replicate_cnt[1] <- j
for (i in 2:length(hx_times)) {
  if (hx_times[i] == hx_times[i - 1]) {
    j <- j + 1
  } else {
    j <- 1
  }
  replicate_cnt[i] <- j
}

# dictionary uses the import_key to associate each measurement with its
# replicate number and HX labeling time
hx_dictionary <- tibble(
  # "import_key" = as.character(0:(length(hx_times) - 1)),
  "import_key" = as.character(1:length(hx_times)),
  "hx_time" = if_else(hx_times != "Full-D", str_sub(hx_times, 1, -2), "Inf"),
  "time_unit" = if_else(hx_times != "Full-D", str_sub(hx_times, -1), "s"),
  "replicate_cnt" = replicate_cnt
)
hx_dictionary$hx_time <- as.double(hx_dictionary$hx_time)
hx_dictionary$import_key <- as.integer(hx_dictionary$import_key)

# THIS NEEDS TO BE REPLACE WITH PIVOT_LONGER BUT IT IS NOT EASY BECAUSE 
# TRIPLETS OF COLUMNS NEED TO BE TIED TOGETHER

# gather and spread using the default variable names
# DEV: gather and spread are now deprecated :)
hx_data_tidy <- hx_data %>%
  gather(key, value, -state, -start, -end, -sequence, -charge) %>%
  separate(key, c("stat", "import_key"), sep = "_") %>%
  mutate(import_key = as.integer(import_key)) %>%
  spread(stat, value)

# bring in the HX times and replicate numbers from the dictionary
hx_data_tidy <- hx_data_tidy %>%
  left_join(hx_dictionary) %>%
  select(-import_key)

# some imported HX-MS data do not include confidence or score
# unreported scores are set to -1 to distinguish from NA
# DEV: problems will arise if the imported data is mixed with some results
# having scores or confidence, but others not
hx_data_tidy <- hx_data_tidy %>% mutate(
 conf = if (!("conf" %in% colnames(hx_data_tidy))) {
   "not reported"}
 else {conf},
 score = if (!("score" %in% colnames(hx_data_tidy))) {
   -1}
 else {score}
)

hx_data_tidy <- hx_data_tidy %>%
  rename(
    hx_sample = state,
    pep_start = start,
    pep_end = end,
    pep_sequence = sequence,
    pep_charge = charge,
    score = score,
    confidence = conf
  )

hx_data_tidy <- hx_data_tidy %>%
  mutate(
    pep_start = as.integer(pep_start),
    pep_end = as.integer(pep_end),
    hx_time = as.double(hx_time),
    replicate_cnt = as.integer(replicate_cnt),
    pep_charge = as.integer(pep_charge),
    confidence = tolower(confidence),
    confidence = factor(
      confidence,
      levels =
        c("low", "medium", "high", "not reported", "theoretical")
    ),
    d = as.double(d),
    score = as.double(score)
  )

# attach peptide ID numbers
hx_data_tidy <- hx_data_tidy %>%
  left_join(peptide_list, by = "pep_sequence") %>%
  mutate(pep_ID = as.integer(pep_ID))

# output files
hx_data_tidy %>% write_csv(out_file_csv)
hx_data_tidy %>% saveRDS(file = out_file_R)