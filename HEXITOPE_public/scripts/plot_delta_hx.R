HX_colors <- c("-2" = "blue", "-1" = "cadetblue1", "0" = "gray", "1" = "#FFFF99", "2" = "yellow")
delta_hx_summary %>% 
  group_by(pep_ID) %>% 
  # avoiding the problem of multiple peptides with same index . . .
  # this step is no longer necessary because this issue is now trapped in delta_hx
  # filter(samp_stdev == min(samp_stdev)) %>%
ggplot() +
 geom_col(mapping = aes(x = pep_ID, y = delta_hx_avg, fill = as.factor(hx_cluster)), color = "black", size = 0.25) +
  scale_fill_manual(values = HX_colors, name = "HX Cluster") +
  geom_line(mapping = aes(x = pep_ID, y = 5*err_prop), linetype = 3) +
  geom_line(mapping = aes(x = pep_ID, y = -5*err_prop), linetype = 3) +
  geom_hline(yintercept = 0) +
  labs(x = "Peptide Index", y = expression(paste(Delta,bar("HX")))) +
  # coord_cartesian(expand = FALSE) +
  coord_cartesian(expand = FALSE, ylim = c(-0.35,0.2)) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
    )
ggsave(str_c("output/", out_file_prefix, "_delta_HX_plot.pdf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix,"_delta_HX_plot.wmf"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_plot.png"), dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_delta_HX_plot.tiff"), dpi = 600)