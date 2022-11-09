library(cluster)
library(tidyverse)
library(stringr)

# WARNING: THIS SCRIPT IS NOW OUT-OF-DATE WITH SOME CHANGES MADE TO HOW USER
# PARAMETERS WERE FORMATTED AND TO HOW THE X-RAY EPITOPE DATA ARE GENERATED.
# THE SCRIPTING SHOULD BE CONSOLIDATED TO ANOTHER
# COMPONENT OF THE CODE. FOR NOW JUST RUN THE PORTION MARKED RUN ME!

# load user parameters
source("scripts/user_parameters.R")

# MOST OF THIS IS SUPERFLOUS IF THIS SCRIPT HAS BEEN RUN WITH ALL OF THE ENVIRONMENT
# VALUES IN MEMORY

# BEGIN SUPERFLOUS-----------------------------------------------------------

seq_file <- str_c("input/", seq_file_name)
AA_seq <- read_csv(seq_file, col_names = FALSE)
AA_seq <- unlist(str_split(AA_seq,""))
protein_data <- tibble(
  res_num = (0 + seq_start):(length(AA_seq) - 1 + seq_start),
  residue = AA_seq
)

epitope_file <- str_c("input/", epitope_file_name)

if (str_sub(epitope_file, -3, -1) == "xml") {
  library(xml2)
  pisa_xml <- read_xml(epitope_file)
  epitope_Xray <- tibble(
    PISA_res = str_sub(str_squish(xml_text(
      xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/STRUCTURE"))), 3, 5),
    res_num = as.integer(str_extract(xml_text(
      xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/STRUCTURE")), "[0-9]+")),
    ASA = as.numeric(xml_text(
      xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/SOLVENTACCESSIBLEAREA"))),
    BSA = as.numeric(xml_text(
      xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/BURIEDSURFACEAREA"))),
    DeltaG =  as.numeric(xml_text(
      xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/SOLVATIONENERGY")))
  )
} else {
  epitope_Xray <- read_fwf(
    epitope_file,
    fwf_cols(
      N = 6,
      I = 3,
      chain = 3,
      PISA_resi = 3,
      res_num = 4,
      junk = 6,
      ASA = 7,
      BSA = 7,
      DeltaG = 7
    ),
    skip = 3
  ) %>%
    select(PISA_resi, res_num, BSA)
}


epitope_Xray <- epitope_Xray  %>% 
  mutate(
    res_num = res_num + epitope_offset,
    epitope = BSA > BSA_threshold
  )


# 
# protein_data <- protein_data %>% 
#   left_join(SASA) %>% 
#   left_join(epitope_Xray)

protein_data <- protein_data %>% 
  left_join(epitope_Xray)

# END SUPERFLOUS------------------------------------------------------

# THIS PART SEEMS TO BE BROKEN, THE SCRIPT DOES NOT HAVE + SEPARATORS
# epitope_list <- epitope_Xray %>% 
#   filter(epitope) %>% 
#   select(res_num) %>% 
#   summarise(resi_list = paste(res_num, collapse = "+"))

# RUN ME! START RUNNING HERE IF ALL VALUES ARE IN THE ENVIRONMENT

epitope_list <- epitope_Xray %>% 
  filter(epitope) %>% 
  pull(res_num) %>% 
  paste(collapse = "+")

pymol_script <- if (file.exists("input/input_script.pml")) {
  str_c(
    scan("input/input_script.pml", character(), quote = "", sep = "\n",quiet = TRUE),
    collapse =   "\r\n"
    )
}

pymol_script <- str_c(pymol_script, 
    str_c(
          "\r\n# automated X-ray epitope script\r\n",
          "color white\r\n",
          "select X_ray_epitope, resi ", epitope_list, "; deselect\r\n",
          "color orange, X_ray_epitope\r\n",
          "deselect\r\n",
          "ray 2000\r\n",
          "png xray_epitope_structure.png\r\n",
          "save xray_epitope_output.pse\r\n",
          "quit\r\n"
      ),
    collapse = "\r\n")
write(pymol_script, "output/xray_epitope_pymol_script.pml", append = FALSE)