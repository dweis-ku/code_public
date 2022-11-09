library(tidyverse)
library(gridExtra)
source("scripts/user_parameters.R")
SASA_plot <- ggplot(data = protein_data) +
 geom_col(mapping = aes(x = res_num, y = replace_na(SASA,0)), na.rm = TRUE, color = "black", size = 0.15) +
  labs(x = "Residue Number", y = expression("SASA (\u212b"^2*")")) +
  scale_x_continuous(breaks = seq(min(tiles$res_num), max(tiles$res_num),by = 20)) +
  coord_cartesian(expand = FALSE) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
  )
BSA_plot <- ggplot(data = protein_data) +
  geom_col(mapping = aes(x = res_num, y = replace_na(BSA,0)), na.rm = TRUE, color = "black", size = 0.15) +
  labs(x = "Residue Number", y = expression("Buried Surface Area (\u212b"^2*")")) +
  scale_x_continuous(breaks = seq(min(tiles$res_num), max(tiles$res_num),by = 20)) +
  coord_cartesian(expand = FALSE) +
  theme(
    aspect.ratio = 0.5,
    panel.grid = element_line(0),
    panel.background = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line()
  )
combined_plot <- grid.arrange(SASA_plot, BSA_plot, nrow = 2)

ggsave(str_c("output/", out_file_prefix, "_SASA_BSA_plot.pdf"), combined_plot, dpi = 300)
# there is a bug in ggplot2 3.3.4 passing a background color argument in emf
# and wmf with ggsave. this is a temporary workaround.
# ggsave(str_c("output/", out_file_prefix, "_SASA_BSA_plot.wmf"), combined_plot, dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_SASA_BSA_plot.png"), combined_plot, dpi = 300)
ggsave(str_c("output/", out_file_prefix, "_SASA_BSA_plot.tiff"), combined_plot, dpi = 600)
# this is another way to put the two plots on the same display
# library(grid)
# pushViewport(viewport(layout = grid.layout(2, 1)))
# print(SASA_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(BSA_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))