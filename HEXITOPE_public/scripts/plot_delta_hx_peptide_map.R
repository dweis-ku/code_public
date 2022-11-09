library(scales)
vlines <- protein_data %>% filter(epitope) %>% pull(res_num)
tiles <- delta_hx_map %>% 
  group_by(pep_ID) %>% 
  # DEV: avoiding the problem of multiple peptides with same index . . .
  # this step is no longer necessary because this issue is now trapped in delta_hx
  # filter(samp_stdev == min(samp_stdev)) %>% 
  select(pep_ID, pep_charge, starts_with("R_")) %>% 
  gather(key = "res_num", value = "delta_hx_avg", -pep_ID, -pep_charge) %>% 
  mutate(res_num = as.integer(str_remove_all(res_num, "R_")))
  delta_hx_limits <-  c(
                  0.1*(floor(10*min(delta_hx_map$delta_hx_avg, na.rm = TRUE))),
                  0.1*(ceiling(10*max(delta_hx_map$delta_hx_avg, na.rm = TRUE)))
                      )
  ggplot() +
  geom_tile(data = tiles, mapping = aes(y = pep_ID, x = res_num, fill = delta_hx_avg)) +
  scale_fill_gradientn(colors = c("blue", "deepskyblue", "white", "lightyellow","yellow"), 
                       name = expression(paste(Delta,bar("HX"))), na.value = "white",
                       limits = delta_hx_limits,
                       values = rescale(c(delta_hx_limits[1],0,delta_hx_limits[2]))
                       ) +
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
      # original size = 1
      size = 0.25
      ) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
  )
ggsave(str_c("output/", out_file_prefix, "_delta_HX_residue_map.pdf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_residue_map.wmf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_residue_map.png"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_residue_map.tiff"), dpi = 600)