library(tidyverse)
library(stringr)
library(xml2)

source("scripts/user_parameters.R")

# compiles xml format output files from PISA protein-protein interface analysis
# for now this is meant for case-by-case pre-processing of the xml files

# which chain in binary interface is the one of interest?
resi_tag <- "RESIDUE2"
resi_tag <- str_c(resi_tag, "//RESIDUE/")

epitope_file_list <- dir(str_c(getwd(), "/input"), pattern = ".xml$", 
                         full.names = TRUE)

# it's really clunky to build the tibble this way, but I just didn't have more
# time to figure out a cleaner way to do this

pisa_xml <- read_xml(epitope_file_list[1])
epitope_Xray <- tibble(
  epitope_file = epitope_file_list[1],
  PISA_resi = str_sub(str_squish(xml_text(
    xml_find_all(pisa_xml, str_c(resi_tag,"STRUCTURE"))
  )), 3, 5),
  res_num = as.integer(str_extract(xml_text(
    xml_find_all(pisa_xml, str_c(resi_tag,"STRUCTURE"))
  ), "[0-9]+")),
  ASA = as.numeric(xml_text(
    xml_find_all(pisa_xml, str_c(resi_tag,"SOLVENTACCESSIBLEAREA"))
  )),
  BSA = as.numeric(xml_text(
    xml_find_all(pisa_xml, str_c(resi_tag,"BURIEDSURFACEAREA"))
  ))
)

first_interface <- epitope_Xray

for (f in epitope_file_list[-1]) {
  pisa_xml <- read_xml(f)
  epitope_Xray <- first_interface %>% add_row(
    epitope_file = f,
    PISA_resi = str_sub(str_squish(xml_text(
      xml_find_all(pisa_xml, str_c(resi_tag,"STRUCTURE"))
    )), 3, 5),
    res_num = as.integer(str_extract(xml_text(
      xml_find_all(pisa_xml, str_c(resi_tag,"STRUCTURE"))
    ), "[0-9]+")),
    ASA = as.numeric(xml_text(
      xml_find_all(pisa_xml, str_c(resi_tag,"SOLVENTACCESSIBLEAREA"))
    )),
    BSA = as.numeric(xml_text(
      xml_find_all(pisa_xml, str_c(resi_tag,"BURIEDSURFACEAREA"))
    ))
  )
}

epitope_Xray <- epitope_Xray %>% group_by(res_num, PISA_resi, ASA) %>%
  summarise(BSA = max(BSA), .groups = "keep")

epitope_Xray <- epitope_Xray  %>%
  mutate(res_num = res_num + epitope_offset,
         epitope = BSA > BSA_threshold)

saveRDS(epitope_Xray, file = "input/A9_epitope_Xray.rds")
rm(epitope_Xray, first_interface)
