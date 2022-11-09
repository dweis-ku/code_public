# ------------------------------------------------------------------------------
# user input

# DATA IMPORT
# name of exported HDExaminer peptide pool results file
hx_data_file <- 
  "input/E1+RTA_one_rep_all.csv"

# name of peptide ID file
peptide_ID_file <- "input/RTA VHH peptide list Toth.csv"

# the amino acid sequence as a simple text file
seq_file_name <- "input/RTA_sequence.txt" 

# file of residue numbers and their SASA values
SASA_file_name <- "input/3srp_SASA.csv"

# offset to add to all residue numbers of the amino acid sequence, zero means
# the first residue in the sequence will be numbered 1
seq_start <- 0 

# DELTA HX CALCULATION
# minimum replicates per HX time for data acceptance
rep_min <- 1

# manual estimate of standard deviation
# set this to NA for standard deviation calculations based on the experimental
# data. Set to a numerical value to override the calculations, for example
# when n = 1 (no replicates). When sd_est is NA and rep_min is 1, all results
# will be treated as statistically significant
sd_est <- 0.25

# reference state, integer 1 or 2 that tells script which set of data is the
# reference state, in other words, the first (1) or second (2) set
ref_state <- 1

# the number of multiples of propagated error for the significance limit
err_limit <- 5 

# one-tail t test p value for Inf HX time lower than max HX time for data
# acceptance
p_inf_limit <- 0.001

# MAPPING OPTIONS
# count first amide (second residue) in pool of exchangeable amides?
count_first_nh <- TRUE

# map altered HX effects at the first amide (second residue)
map_first_nh <- FALSE

# minimum SASA value for surface residue
SASA_threshold <- 2.5

# CLUSTERING DELTA HX
# number of clusters of the absolute value of delta_hx
n_clus <- 3

# SOLVED EPITOPE OPTIONS
# the following options are only relevant if a pre-defined epitope will be
# loaded for comparison such as from X-ray crystallography

# optional: epitope interface file (PISA)
# enter NA if no epitope data will be included
epitope_file_name <- "input/E1_5BOZ_PISA.xml"

# offset to apply to the first residue of the amino acid sequence of the epitope
# file, zero means the first residue in the sequence will be numbered 1
epitope_offset <- 0

# threshold of buried surface area for a residue to be classified as epitope
BSA_threshold <- 2.5
# ------------------------------------------------------------------------------