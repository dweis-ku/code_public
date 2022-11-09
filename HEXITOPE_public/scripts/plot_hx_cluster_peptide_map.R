# DEV: Error in seq.default(0, max(tiles$pep_ID) + 1, by = 20) : 
# 'to' must be a finite number
# when there are peptides in the data that don't have a pep_ID because 
# they are not on the list

HX_colors <- c("-2" = "blue", "-1" = "cadetblue1", "0" = "gray", "1" = "#FFFF99", "2" = "yellow", "NA" = "white")
vlines <- protein_data %>% filter(epitope) %>% pull(res_num)
tiles <- hx_cluster_map %>% 
  group_by(pep_ID) %>% 
  # DEV: avoiding the problem of multiple peptides with same index . . .
  # this step is no longer necessary because this issue is now trapped in delta_hx
  # filter(samp_stdev == min(samp_stdev)) %>% 
  select(pep_ID, pep_charge, starts_with("R_")) %>% 
  gather(key = "res_num", value = "hx_cluster", -pep_ID, -pep_charge) %>% 
  mutate(res_num = as.integer(str_remove_all(res_num, "R_")),
         hx_cluster = ifelse(!is.na(hx_cluster) | is.nan(hx_cluster), as.character(hx_cluster), "NA"),
         hx_cluster = factor(hx_cluster, 
                             levels = c(min(delta_hx_summary$hx_cluster, na.rm = TRUE):max(delta_hx_summary$hx_cluster, na.rm = TRUE),"NA","NaN"))
         ) 
  ggplot() +
  geom_tile(data = tiles, mapping = aes(y = pep_ID, x = res_num, fill = hx_cluster)) +
  scale_fill_manual(values = HX_colors, name = "HX Cluster") +
  # DEV: this next bit is incomplete but works OK in the present application
  scale_x_continuous(breaks = seq(min(tiles$res_num), max(tiles$res_num),by = 20)) +
  scale_y_continuous(breaks = seq(0, max(tiles$pep_ID) + 1, by = 20)) + 
  # 
  geom_rect(data = delta_hx_map, aes(xmin = pep_start - 0.5, xmax = pep_end + 0.5,
                                     ymin = pep_ID - 0.5, ymax = pep_ID + 0.5),
            # fill = NA, color = "black", size = 0.5) +
            fill = NA, color = "black", size = 0.15) +
  labs(x = "Residue Number", y = "Peptide Index") +
  coord_cartesian(expand = FALSE) +
  geom_vline(
      xintercept = vlines,
      linetype = 1,
      color = "orange",
      show.legend = TRUE,
      # size = 1
      size = 0.25
      ) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
  )
ggsave(str_c("output/", out_file_prefix, "_hx_cluster_peptide_map.pdf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_hx_cluster_peptide_map.wmf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_hx_cluster_peptide_map.png"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_hx_cluster_peptide_map.tiff"), dpi = 600)