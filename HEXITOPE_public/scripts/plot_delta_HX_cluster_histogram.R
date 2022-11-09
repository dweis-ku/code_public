delta_hx_summary %>% 
  group_by(pep_ID) %>% 
  filter(!is.na(delta_hx_avg)) %>% 
  # avoiding the problem of multiple peptides with same index . . .
  filter(samp_stdev == min(samp_stdev, na.rm = TRUE)) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = delta_hx_avg), color = "black", fill = "gray", 
                 binwidth = abs((min(delta_hx_summary$delta_hx_avg, na.rm = TRUE) - max(delta_hx_summary$delta_hx_avg, na.rm = TRUE))/25)) + 
  geom_vline(xintercept = tail(head(cuts, -1),-1), linetype = 2) +
  labs(x = expression(paste(Delta,bar("HX"))), y = "Count") +
  coord_cartesian(xlim = 
                    round(
                      c(min(1.1*delta_hx_summary$delta_hx_avg, na.rm = TRUE),
                        max(1.1*delta_hx_summary$delta_hx_avg, na.rm = TRUE))
                        ,1),
                  expand = FALSE) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
  )
ggsave(str_c("output/", out_file_prefix, "_delta_HX_cluster_plot.pdf"), dpi = 300)
# there is a bug in ggplot2 3.3.4 passing a background color argument in emf
# and wmf with ggsave. this is a temporary workaround.
# ggsave(str_c("output/", out_file_prefix, "_delta_HX_cluster_plot.wmf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_cluster_plot.png"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_cluster_plot.tiff"), dpi = 600)
