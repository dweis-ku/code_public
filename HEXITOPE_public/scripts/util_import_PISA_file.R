epitope_file <- str_c("input/", epitope_file_name)

# WIP: this utility will read a PISA file
# it should ultimately generate an RDS file
# DEV: this is for a very specific PISA file structure; needs updating to be
# smarter

# 
 if (!is.na(epitope_file_name)) {
   if (str_sub(epitope_file,-3,-1) == "xml") {
     library(xml2)
     pisa_xml <- read_xml(epitope_file)
     epitope_Xray <- tibble(
       PISA_resi = str_sub(str_squish(xml_text(
         xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/STRUCTURE")
       )), 3, 5),
       res_num = as.integer(str_extract(xml_text(
         xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/STRUCTURE")
       ), "[0-9]+")),
       ASA = as.numeric(xml_text(
         xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/SOLVENTACCESSIBLEAREA")
       )),
       BSA = as.numeric(xml_text(
         xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/BURIEDSURFACEAREA")
       )),
       DeltaG =  as.numeric(xml_text(
         xml_find_all(pisa_xml, "RESIDUE2//RESIDUE/SOLVATIONENERGY")
       ))
     )
   } else {
# legacy for a white-space format
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
     mutate(res_num = res_num + epitope_offset,
            epitope = BSA > BSA_threshold)
 }