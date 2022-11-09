library(tidyverse)
library(stringr)
library(stats)

# load user parameters
source("scripts/user_parameters.R")

if(rep_min < 2){
  message(c("ATTENTION: You have set minimum replicates to 1. ",
            "Significance testing and uncertainty propagation will be disabled " ,
            "unless an estimated value is provided for standard deviation estimate",
            " using the sd_est parameter."))
  invisible(readline(prompt="Press [enter] to continue."))
}

count_amide <- function(x, ignore.second = NULL) {
  # count backbone amides, ignore proline and first and/or second residue
  # return NA for sequences containing unrecognized amino acids
  amino_acids = "ACDEFGHIKLMNPQRSTVWY"
  white_list = str_c("^[", amino_acids, "]+$")
  library(stringr)
  x <- as.character(x)
  if (is.null(ignore.second))
    ignore.second <- FALSE
  stopifnot(length(ignore.second) == 1)
  ifelse(str_detect(x, white_list), {
    start <- if (ignore.second) 3 else 2
    amide_seq <- str_remove_all(str_sub(x, start), "P")
    str_length(amide_seq)
  },
  returnValue(NA))
}

# DEV: this function should be made more robust
stdev_pool <- function(s, n) {
  sqrt(
    sum(s ^ 2 * (n - 1)) / sum(n - 1)
  )
}

if (all(unique(hx_data_tidy$hx_time) != Inf)){
theo_NH <- hx_data_tidy %>% 
  distinct(hx_sample, pep_sequence, pep_charge, .keep_all = TRUE) %>% 
  mutate(
    hx_time=Inf, 
    d = count_amide(pep_sequence, ignore.second = !count_first_nh),
    confidence = "theoretical",
    score = NA
    )
# theoretical data must also have the minimum number of replicates
theo_NH <- theo_NH %>% slice(rep(1:n(), each = 2))
hx_data_tidy <- add_row(hx_data_tidy, theo_NH)
hx_data_tidy <- hx_data_tidy %>% 
  arrange(hx_sample,pep_ID,pep_charge,hx_time,replicate_cnt)
}

# end insertion

# ------------------------ BEGIN DATA FILTERING --------------------------------
# preserves the full list of peptides before any subsequent filtering
all_peptides <- hx_data_tidy %>% 
  select(pep_ID, pep_start, pep_end, pep_sequence) %>% 
  group_by(pep_ID, pep_start, pep_end, pep_sequence) %>% 
  summarise()

# number of HX times 
hx_time_cnt <- length(unique(hx_data_tidy$hx_time))
# longest HX time != Inf
max_time <- hx_data_tidy$hx_time %>% .[. != Inf] %>% unique() %>% max()

# filter out HX data that do not have the minimum number of replicates
# and that do not have complete HX time course in both protein states

hx_data_filtered <- hx_data_tidy %>% 
  filter(!is.na(d)) %>% 
  group_by(hx_sample, pep_ID, pep_sequence, pep_charge, hx_time) %>%
  mutate(n = n()) %>% 
  filter(n >= rep_min) %>% 
  select(-n) %>% 
  group_by(hx_sample, pep_ID, pep_sequence, pep_charge) %>% 
  filter(length(unique(hx_time)) == hx_time_cnt) %>% 
  group_by(pep_ID, pep_charge) %>% 
  mutate(
    n = length(unique(hx_sample))
    ) %>% 
  filter(n == 2) %>% 
  select(-n)

# IN DEVELOPMENT: added the if gate here 
# remove peptides that have an infinite HX that is significantly less than the
# longest hx_time by a one-tailed t-test if there are experimental data for 
# maximum HX

if (!any(select(filter(hx_data_tidy, hx_time == Inf), confidence) == 
         "theoretical")){
  longest_hx_data <- hx_data_filtered %>%
    filter(hx_time == max_time)
  inf_hx_data <- hx_data_filtered %>%
    filter(hx_time == Inf)
  
  inf_hx_test <- full_join(longest_hx_data, # this is x
                           inf_hx_data, # this is y
                           by = c("hx_sample", 
                                  "pep_ID", 
                                  "pep_charge", 
                                  "replicate_cnt")
                           ) %>%
    group_by(hx_sample, pep_ID, pep_charge) %>%
# Set use t-test if there is more than one replicate, otherwise set the t-test
# p value to infinity
        summarise(
          p_inf = ifelse(rep_min > 1, 
                         t.test(d.x, d.y, alternative = "greater")$p.value,
                         Inf)
              )
  
  hx_data_filtered <- hx_data_filtered %>%
    group_by(hx_sample, pep_ID, pep_charge) %>%
    left_join(inf_hx_test) %>%
    filter(p_inf > p_inf_limit) %>% 
    select(-p_inf) %>% 
    ungroup()
}
# ------------------------ BEGIN DATA SUMMARY/RESHAPE---------------------------
# averages and reshapes the data
# build a tibble of averaged data with data from both protein states 
# at each HX time treated as variables for each peptide
# add the number of amide H
# compute the summary stats: mean, standard deviation, and number of
# observations

hx_by_time <- hx_data_filtered %>%
  mutate(n_amide = count_amide(pep_sequence, count_first_nh)) %>%
  group_by(hx_sample,
           pep_ID,
           pep_charge,
           pep_start,
           pep_end,
           pep_sequence,
           n_amide,
           hx_time) %>%
  summarize(
    d_avg = mean(d, na.rm = TRUE),
    stdev = sd(d, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# get the vector of infinite HX data (assumes that same control used in both
# states)
one_sample <- unique(hx_by_time$hx_sample)[1]
inf_hx_data <- hx_by_time %>% 
  filter(
    hx_time == Inf, 
    hx_sample == one_sample
  ) %>%
  ungroup() %>% 
  select(-hx_sample, -hx_time, -pep_start, -pep_end, -pep_sequence, -n_amide) 

# modify the tibble to remove the infinite HX data as observations/rows
# temporarily unite the summary statistics to make it easier to rename the
# variables
sample_cols <- unique(hx_by_time$hx_sample)
hx_by_time <- hx_by_time %>%
  filter(hx_time != Inf) %>% 
  unite(united_var, d_avg, stdev, n) %>% 
  group_by(hx_sample,
           pep_ID,
           pep_charge,
           pep_start,
           pep_end,
           pep_sequence,
           n_amide,
           hx_time) %>% 
  spread(key = hx_sample, value = united_var) %>%
  ungroup() 

# rename the protein states to ref and samp
new_name <- c("ref","samp")
if (ref_state == 2) {
  new_name <- rev(new_name)
}

# standardize the variables names with ref and samp prefixes
# separate the united variables
for (i in 1:2) { 
  hx_by_time <- separate(
    hx_by_time, 
    sample_cols[i], 
    into = 
      str_c(new_name[i], c("_d_avg", "_stdev", "_n")),
    sep = "_",
    convert = TRUE
  )
}

# add the infinite HX as a variable
hx_by_time <- hx_by_time %>% 
  left_join(
    inf_hx_data,
    by = c("pep_ID", "pep_charge")
    ) %>% 
  rename(
    inf_d_avg = d_avg,
    inf_stdev = stdev,
    inf_n = n
  )

# filtering data for simultaneous completeness in both states
hx_by_time <- hx_by_time %>% 
  filter(!is.na(samp_d_avg), !is.na(ref_d_avg))

#----------------BEGIN COMPUTE DELTA HX AND PROPAGATED UNCERTAINTY--------------
delta_hx_summary <- hx_by_time %>%
  group_by(
    pep_ID, 
    pep_charge,
    pep_start,
    pep_end,
    pep_sequence,
    n_amide
    ) %>%
  summarise(
    delta_hx_avg =
      # all the d_avg_inf values in the group are the same, summarize does not
      # want a vector
      (mean(samp_d_avg) - mean(ref_d_avg)) / inf_d_avg[1], 
      samp_stdev = stdev_pool(samp_stdev, samp_n),
      samp_n_tot = sum(samp_n),
      samp_n_recip_sum = sum(1 / samp_n),
      ref_stdev = stdev_pool(ref_stdev, ref_n),
      ref_n_tot = sum(ref_n),
      ref_n_recip_sum = sum(1 / ref_n),
      inf_d_avg = inf_d_avg[1],
      inf_stdev = inf_stdev[1],
      inf_n = inf_n[1],
      n_hx_time = n(),
      # pooled standard deviations will be NaN for n = 1 replicates
      # replace with a manual value if defined, otherwise stays as NA
      samp_stdev = if(is.na(samp_stdev)){sd_est},
      ref_stdev = if(is.na(ref_stdev)){sd_est},
      inf_stdev = if(is.na(inf_stdev)){sd_est}
      ) %>% 
  filter(!is.na(delta_hx_avg)) %>% 
# propagate uncertainty based on experimental standard deviation or based on
# sd_est value, if supplied
  mutate(
    err_prop = 
      ifelse(
        is.na(sd_est),
    # based on the experimental data 
    # DEV: this case should be revised to use fully
    # pooled standard deviation across the entire data set
             1 / inf_d_avg * 
      sqrt(
        (ref_stdev ^ 2 * ref_n_recip_sum + samp_stdev ^ 2 * samp_n_recip_sum) 
        / n_hx_time ^ 2 + delta_hx_avg ^ 2 * inf_stdev ^ 2 / inf_n),
    # based on an estimated value for pooled standard deviation 
      sd_est/(n_hx_time * inf_d_avg) *
      sqrt(ref_n_recip_sum + samp_n_recip_sum + delta_hx_avg ^ 2 / inf_n)
      )
  ) %>% 
  right_join(all_peptides) %>% 
  select(
    -contains("recip") )


# handles single measurement use case without sd_est value
# reset NaN propagated uncertainty to zero if rep_min is 1 
# this means that all results will be treated as significant
delta_hx_summary <- delta_hx_summary %>% mutate(
  err_prop = ifelse(rep_min == 1 & is.na(err_prop), 0, err_prop)
)

# ----------- BEGIN SIGNIFICANCE TESTING ---------------------------------------
# Test for significance of DeltaHX at the peptide level
delta_hx_summary <- delta_hx_summary %>% mutate(
  significant = abs(delta_hx_avg) > err_limit * err_prop
  )
# ----------- HANDLE CASES WHERE PEPTIDE IS REPORTED IN MORE THAN ONE CHARGE ---
# the best result will be the one with smallest propagated uncertainty
# in the no-replicates case, choose the one with most extreme delta_hx_avg value
# this step retains the "best" observation of each peptide, removing others
# except for the unlikely event that two charge states have the same propagated
# uncertainty and also the same delta_hx value
delta_hx_summary <- delta_hx_summary %>% 
  group_by(pep_ID) %>% 
  filter(err_prop == min(err_prop)) %>% 
  filter(abs(delta_hx_avg) == max(abs(delta_hx_avg))) %>% 
  ungroup()

# DEV: This is sloppy but does the job of preventing multiple charge
# states of the same peptide
max_pep_ID <- delta_hx_summary %>% group_by(pep_ID) %>% count() %>% 
  pull(n) %>% max()
if(max_pep_ID != 1){
  message("ERROR: Multiple instances of same peptide detected. Halting.")
  stop
}