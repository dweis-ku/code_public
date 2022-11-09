# cluster delta_hx_avg into n_clus clusters based on absolute values using
# kmeans algorithm then bin the delta_hx_avg values into 2n-1 clusters by
# preserving the sign the middle cluster will be centered on zero

library(cluster)
library(tidyverse)
library(stringr)
# load user parameters
source("scripts/user_parameters.R")

# build a residue-centric tibble

AA_seq <- read_lines(seq_file_name)
AA_seq <- unlist(str_split(AA_seq,""))
protein_data <- tibble(
  res_num = (0 + seq_start):(length(AA_seq) - 1 + seq_start),
  residue = AA_seq
)

# append X-ray crystallography data for surface residues

SASA <- read_csv(SASA_file_name)
SASA <- SASA %>% mutate(surf_resi = SASA > SASA_threshold)
protein_data <- protein_data %>% left_join(SASA)

# append a pre-defined epitope, if known
# THIS LINE IS A BUG FIX
if (is.na(epitope_file_name)) {epitope_file_name <- ""}
# END BUG FIX
if (str_sub(epitope_file_name, -3, -1) == "rds") {
  epitope_Xray <- readRDS(epitope_file_name)
} else {
  # if you don't show up with a pre-processed epitope file, you don't get an
  # epitope
  epitope_Xray <- tibble(
    res_num = protein_data$res_num,
    PISA_resi = NA,
    BSA = NA,
    epitope = FALSE
  )
}

protein_data <- protein_data %>% left_join(epitope_Xray)

# --------------BEGIN CLUSTERING ----------------------------------------------

x <- delta_hx_summary %>%
  filter(!is.na(delta_hx_avg)) %>%
  pull(delta_hx_avg)

ctrs <- sort(kmeans(abs(x), centers = n_clus)$center)

cuts <- sort(c(0.5 * (ctrs[-length(ctrs)] + ctrs[-1]),
               -0.5 * (ctrs[-length(ctrs)] + ctrs[-1]),
               -Inf,
               Inf))

delta_hx_summary <- delta_hx_summary %>%
  mutate(hx_cluster =
           as.double(cut(
             delta_hx_avg,
             breaks = cuts,
             labels = -(length(ctrs) - 1):(length(ctrs) - 1)
           ))
         - n_clus)

# --------------BEGIN MAPPING HX EFFECTS ONTO SEQUENCE-------------------------

# makes one column for each residue, each cell contains the residue number
add_resi <- function(x, start, end) {
  for (i in start:end) {
    x <- x %>%
      add_column(!!paste0("R_", as.character(i)) := i)
  }
  x
}

# generate a peptide-by-residue mapping structure
# by attaching a residue-by-residue map to the peptide-centric data

hx_cluster_map <- delta_hx_summary
hx_cluster_map <-
  add_resi(hx_cluster_map,
           min(protein_data$res_num),
           max(protein_data$res_num))

pep_map <- function(i, j, k, map.second = NULL)
#  maps 1 when i is between j and k limits or NA if outside limits with the
#  logic of HX
  {
    if (is.null(map.second))
      map.second <- TRUE
    offset <- if (map.second) 1 else 2
    ifelse((i >= j + offset) & (i <= k), 1, NA)
  }

resi_in_peptide <- function(i, j, k, map.second = NULL)
#  returns true when i is between j and k limits with the logic of HX
{
  if (is.null(map.second))
    map.second <- TRUE
  offset <- if (map.second) 1 else 2
  ifelse((i >= j + offset) & (i <= k), TRUE, FALSE)
}

proline_list <- protein_data %>% 
  filter(residue == "P") %>% 
  select(res_num) %>% 
  pull()

surf_resi <- protein_data %>% 
  filter(surf_resi) %>% 
  select(res_num) %>% 
  pull()

# map peptide HX cluster effect onto the residues 
# require that the residue is not proline (proline_list) 
# require that the residue is on the surface (surf_resi)
# the second residue in a peptide will be mapped or not depending on
# map_first_nh as set in user parameters
# mark residue as NA in the peptide if it does not meet these criteria 
# require that the absolute value of the HX effect exceeds the significance
# limit assign to zero if not significant

# BUG HERE: the significance testing is not being applied correctly
# in V1C7 there are -1 cluster peptides that do not exceed the significance limit

# BEGIN BUG FIX HERE
# the fix involved adding significance testing in delta_hx

hx_cluster_map <- hx_cluster_map %>%
  ungroup %>%
  mutate_at(vars(starts_with("R_")),
            ~ (
              # this steps applies the logic as NA or 1
              ifelse(
                !resi_in_peptide(
                  .,
                  hx_cluster_map$pep_start,
                  hx_cluster_map$pep_end,
                  map_first_nh
                ) |
                  . %in% proline_list |
                  !. %in% surf_resi,
                NA,
                1
              )
              * hx_cluster_map$hx_cluster * hx_cluster_map$significant
            ))

# END BUG FIX HERE, original code follows, it has been commented out

# hx_cluster_map <- hx_cluster_map %>%
#   ungroup %>%
#   mutate_at(vars(starts_with("R_")),
#             ~ (
#               # this steps applies the logic as NA or 1
#               ifelse(
#                 !resi_in_peptide(
#                   .,
#                   hx_cluster_map$pep_start,
#                   hx_cluster_map$pep_end,
#                   map_first_nh
#                 ) |
#                   . %in% proline_list |
#                   !. %in% surf_resi,
#                 NA,
#                 1
#               )
#               *
#                 # this part multiplies the logic by either zero or the HX
#                 # cluster depending on significance limit
#                 ifelse(
#                   abs(hx_cluster_map$delta_hx_avg) <
#                     err_limit * hx_cluster_map$err_prop,
#                   0,
#                   hx_cluster_map$hx_cluster
#                 )
#             ))

extreme_value <- function(x){
  # returns the value farthest from zero, NaN for ties unless tie is zero
  if (all(is.na(x))) return(NA)
  if (length(unique(x[!is.na(x)])) == 1) return(unique(x[!is.na(x)]))
  if (abs(min(x, na.rm = TRUE)) == abs(max(x, na.rm = TRUE))) return(NaN)
  if (abs(min(x, na.rm = TRUE)) < abs(max(x, na.rm = TRUE))) {
    max(x, na.rm = TRUE) }
  else{
    min(x, na.rm = TRUE)
  }
}

# HX cluster mapping
# generate a residue-centric map based on the most extreme cluster measured
# for each residue
hx_cluster_by_resi <- hx_cluster_map %>% 
  summarise_at(vars(starts_with("R_")), ~ extreme_value(.)
  ) %>% 
  gather(key = "res_num", value = "hx_cluster") %>% 
  mutate(
    res_num = as.integer(str_remove(res_num, "R_"))
  ) %>% 
  right_join(protein_data) %>% 
  mutate(
    hx_class = 
      ifelse(SASA & residue != "P",
             hx_cluster,
             NA)
  )

# HX effect mapping
# map peptide HX effect onto the residues 
# generate a peptide-by-residue mapping structure
# by attaching a residue-by-residue map to the peptide-centric data
delta_hx_map <- delta_hx_summary
delta_hx_map <- add_resi(delta_hx_map, min(protein_data$res_num), 
                         max(protein_data$res_num))

# transform the map into the HX effect values
# require that the residue is not proline (proline_list) 
# the second residue in a peptide will be mapped or not depending on
# map_first_nh as set in user parameters
# mark residue as NA in the peptide if it does not meet these criteria 
# require that the absolute value of the HX effect exceeds the significance
# limit assign to zero if not significant
delta_hx_map <- delta_hx_map %>%
  ungroup %>%
  mutate_at(
    vars(starts_with("R_")),
    ~ (
      delta_hx_map$delta_hx_avg
      * pep_map(
        .,
        delta_hx_map$pep_start,
        delta_hx_map$pep_end,
        map_first_nh
      )
      *
        ! (. %in% proline_list)
    )
  )

# HX cluster mapping
# generate a residue-centric map based on the most extreme HX effect measured
# for each residue
delta_hx_by_resi <- delta_hx_map %>% 
  summarise_at(vars(starts_with("R_")), ~ extreme_value(.)
  ) %>% 
  gather(key = "res_num", value = "delta_hx_avg") %>% 
  mutate(
    res_num = as.integer(str_remove(res_num, "R_"))
  ) %>% 
  right_join(protein_data) %>% 
  mutate(
    delta_hx_avg = 
      ifelse(surf_resi & residue != "P",
             delta_hx_avg,
             NA)
  )

# --------------BEGIN OUTPUT FILE GENERATION------------------------------------

out_file_prefix <- str_remove(hx_data_file, ".csv")
out_file_prefix <- str_remove(out_file_prefix, "input/")
out_file_resi <- str_c("output/", out_file_prefix, "_resi_list.csv")
hx_cluster_by_resi %>% write_csv(out_file_resi)

output <- hx_cluster_by_resi %>% 
  group_by(hx_cluster) %>% 
  summarise(
    resi_list = paste(res_num, collapse = "+")
  )

out_file_clus <- str_c("output/", out_file_prefix, "_clustered_resi.csv")
output %>% write_csv(out_file_clus)

delta_hx_summary %>% 
  write_csv(str_c("output/",out_file_prefix,"_peptide_summary.csv"))
delta_hx_summary %>% 
  saveRDS(file = str_c("output/",out_file_prefix,"_peptide_summary.RData"))
save.image(str_c("output/",out_file_prefix,"_workspace.RData"))

# if you need n_clus other than 3, you will need to define a different scheme or
# accept this one don't blame me if your protein looks like crap . . .
if (n_clus == 3) {
  pymol_colors <-  c("blue","cyan","white", "paleyellow","yellow","gray20",
                     "magenta")
} else {
  pymol_colors <- c("red", "orange", "yellow", "green", "blue", "cyan", 
                    "magenta", "teal", "salmon", "olive", "hotpink", 
                    "gray20")[1:(2 * n_clus)]
}

script_source <- tibble(
  hx_cluster = c(-(n_clus - 1):(n_clus - 1), NA, NaN),
  pymol_colors = pymol_colors
) %>%
  left_join(output) %>% 
  mutate(
    hx_cluster = ifelse(!is.na(hx_cluster) | is.nan(hx_cluster), 
                        as.character(hx_cluster), "NA")
  ) %>% 
  filter(!is.na(resi_list))

pymol_script <- if (file.exists("input/input_script.pml")) {
  str_c(
    scan("input/input_script.pml", character(), quote = "", sep = "\n",
         quiet = TRUE),
    collapse =   "\r\n"
    )
}
# THERE IS A PROBLEM HERE IF THE PDB FILE CONTAINS MORE THAN ONE CHAIN
# THE SYNTAX HERE WILL APPLY THE CLUSTER DATA TO ALL RESIDUES WITH THIS NUMBER
# NEED TO MODIFY AS SOMETHING LIKE
# "select ", script_source$hx_cluster, "_HX_cluster, ", antingen_chain_name, script_source$resi_list
# where antigen_chain_name is something like "/1GGQ//A/"
pymol_script <- str_c(pymol_script, 
  str_c(
  "# automated delta_hx_cluster script\r\n",
    str_c(
          str_c("select ", script_source$hx_cluster, "_HX_cluster, resi ", 
                script_source$resi_list, "; deselect"),
          str_c("color ", script_source$pymol_colors, ", ", 
                script_source$hx_cluster, "_HX_cluster"),
          sep = "\r\n", 
          collapse = "\r\n"
          ),
      "\r\ndeselect\r\n",
      "hide everything\r\n",
      "show surface\r\n",
      "ray 2000\r\n",
      "png hx_output_structure_surface.png\r\n",
      "hide everything\r\n",
      "show cartoon\r\n",
      "ray 2000\r\n",
      "png hx_output_structure_cartoon.png\r\n",
      "save hx_output.pse\r\n",
      "quit\r\n")
      , sep = "\r\n"
      )
script_file_name <- str_c("output/", out_file_prefix,
                          "_hx_epitope_pymol_script.pml")
write(pymol_script, script_file_name, append = FALSE)