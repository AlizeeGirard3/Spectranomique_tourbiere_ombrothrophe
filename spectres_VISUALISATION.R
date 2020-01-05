# 13 août

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                        Visualisation des spectres 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)
library(tidyr)
library("RColorBrewer")
library(ggpubr) # ggarrange
library(spectrolab)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

setwd("~/Documents/Maîtrise/DonnéesAnalyses/spectres")

load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_large.RData") # rt_sg_large
# obtenu via manips sur script "nettoyage_spectres_lissage_correction.R" dans le document "scripts-NETTOYAGE"


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Manipulations préalables ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# POUR SPECTRES ABSOLUS 
refl <- rt_sg_large %>% 
  filter(propriete == "reflectance") %>% 
  gather(key = "wavelength", value = "valeur_sg", c("400":"2400"), -sample_id, -scientific_name, -site_class)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Visualisation ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

## moyenne par espèce, spectres absolus
refl_mean_sp <- refl %>% 
  group_by(scientific_name, wavelength) %>% 
  summarize(mean_value = mean(valeur_sg), 
            STD_value = sd(valeur_sg)) %>% 
  ungroup()
refl_mean_sp$wavelength <- as.numeric(paste(refl_mean_sp$wavelength))

## moyenne par espèce, spectres NORM
refl_mean_sp <- refl %>%
  group_by(scientific_name, wavelength) %>%
  summarize(mean_value = mean(valeur_sg),
            STD_value = sd(valeur_sg)) %>%
  ungroup()
refl_mean_sp$wavelength <- as.numeric(paste(refl_mean_sp$wavelength))

## mean+sd spectra per species
refl_mean_sp$sc = factor(refl_mean_sp$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                                  "Kalmia angustifolia", 'Rhododendron groenlandicum'), 
                         labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum")) 

regions <- data_frame(Region = factor(c('VIS', 'NIR', 'SWIR1', 'SWIR2'),
                                      levels = c('VIS', 'NIR', 'SWIR1', 'SWIR2')),
                      xmin = c(400, 720, 1530, 2000),
                      xmax = c(700, 1400, 1900, 2400),
                      ymin = 0,
                      ymax = 1)
# voir ?spectrolab::plot_regions() pour les régions selon Anna

regions_text <- data_frame(Region = factor(c('VIS', 'NIR', 'SWIR1', 'SWIR2'),
                                           levels = c('VIS', 'NIR', 'SWIR1', 'SWIR2')),
                           xmin = c(530, 1060, 1715, 2200))

p <- ggplot(refl_mean_sp) +
  geom_rect(data = regions, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.3) +
  geom_text(data = regions_text, aes(x = xmin, y = 0.026, label = Region), size = 3) +
  geom_line(aes(x = wavelength, y = mean_value, color = sc)) +
  geom_ribbon(data = refl_mean_sp, aes(fill = sc,
                                       x = wavelength,
                                       ymin = mean_value - STD_value, 
                                       ymax = mean_value + STD_value), 
              alpha = 0.3) +
  labs(y = 'Reflectance', x = 'Wavelength (nm)') +
  coord_cartesian(ylim = c(0, 0.7), expand = F) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  scale_fill_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  theme_light() +
  theme(legend.position = c(0.83, 0.78), legend.text = element_text(face = "italic"),
        panel.border = element_rect(colour = brewer.pal(9, "Greys")[5], fill = NA, size = 1),
        panel.grid = element_blank())
  # theme(legend.position = "none", legend.text = element_text(face = "italic"),
  #       panel.border = element_rect(colour = brewer.pal(9, "Greys")[5], fill = NA, size = 1),
  #       panel.grid = element_blank())
p
# ggsave('~/Desktop/spectra_sans_zone.png', width = 4, height = 3)
ggsave('figures/spectra_sp.png', width = 6, height = 4)


## mean+sd spectra per site separated species
# {
#   # noms pour graphiques
#   refl$sc = factor(refl$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
#                                                               "Kalmia angustifolia", 'Rhododendron groenlandicum'), 
#                         labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum")) 
#   
#   refl$site_ordered = factor(refl$site_class, levels = c("GPB_canal", "GTV_SE", "MBP_open", "MBP_Intermediate", "MBP_High"), 
#                                   labels = c("GPB_canal", "GTV_SE", "Background", "Intermediate", "High")) 
#   
#   refl$wavelength <- as.numeric(as.character(refl$wavelength))
#   
#   # définir les boîtes zoomées
#   library(grid) # for unit
#   library(gridExtra) # for grid.arrange
#   
#   fulltheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                      panel.background = element_blank(), plot.title = element_text(hjust = 0), 
#                      legend.position = "none", plot.margin = margin(0, 0, 0, 0, "mm")) 
#   
#   zoomtheme <- theme(legend.position="none", axis.line = element_blank(),axis.text.x=element_blank(),
#                      axis.text.y=element_blank(),axis.ticks=element_blank(),
#                      axis.title.x=element_blank(),axis.title.y=element_blank(),
#                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                      panel.background = element_rect(color = 'black', fill = "white"), 
#                      plot.title = element_text(hjust = 0.5), plot.margin = margin(0, 0, 0, 0, "mm"))
#   
#   # créer objets
#   refl_mean_site_EV <- refl %>% 
#     dplyr::filter(site_ordered %in% c("Background", "Intermediate", "High")) %>% 
#     dplyr::filter(sc == "E. vaginatum") %>% 
#     group_by(site_ordered, wavelength) %>% 
#     summarize(mean_value = mean(valeur_sg), 
#               STD_value = sd(valeur_sg)) 
#   
#   refl_mean_site_CC <- refl %>% 
#     dplyr::filter(site_ordered %in% c("Background", "Intermediate", "High")) %>% 
#     dplyr::filter(sc == "C. calyculata") %>%
#     group_by(site_ordered, wavelength) %>% 
#     summarize(mean_value = mean(valeur_sg), 
#               STD_value = sd(valeur_sg))
#   
#   refl_mean_site_KA <- refl %>% 
#     dplyr::filter(site_ordered %in% c("Background", "Intermediate", "High")) %>% 
#     dplyr::filter(sc == "K. angustifolia") %>% 
#     group_by(site_ordered, wavelength) %>% 
#     summarize(mean_value = mean(valeur_sg), 
#               STD_value = sd(valeur_sg)) 
#   refl_mean_site_RG <- refl %>% 
#     dplyr::filter(site_ordered %in% c("Background", "Intermediate", "High")) %>% 
#     dplyr::filter(sc == "R. groenlandicum") %>%
#     group_by(site_ordered, wavelength) %>%
#     summarize(mean_value = mean(valeur_sg), 
#               STD_value = sd(valeur_sg)) 
#   
#   CC.full <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(ylim = c(0.0015, 0.0390)) +
#     labs(y = 'Normalized reflectance', x = '', title = bquote(italic("C. calyculata"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     theme_bw() +
#     fulltheme
#   CC.zoom <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(xlim = c(475, 675), ylim = c(0.004, 0.0118)) +
#     labs(title = bquote(bold("VIS"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     zoomtheme
#   CC.zoom.grob <- ggplotGrob(CC.zoom)
#   CC.zoom.grob <- CC.full + 
#     annotation_custom(grob = CC.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.0245, ymax = 0.0399)
#   # ggsave('figures/spectra_CC_site.png', CC.zoom.grob, width = 6, height = 5)
#   
#   EV.full <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(ylim = c(0.0015, 0.0390)) +
#     labs(y = '', x = '', title = bquote(italic("E. vaginatum"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     theme_bw() +
#     fulltheme
#   EV.zoom <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(xlim = c(475, 675), ylim = c(0.004, 0.0108)) +
#     labs(title = bquote(bold("VIS"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     zoomtheme
#   EV.zoom.grob <- ggplotGrob(EV.zoom)
#   EV.zoom.grob <- EV.full + 
#     annotation_custom(grob = EV.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.026, ymax = 0.0399)
#   # ggsave('figures/spectra_EV_site.png', EV.zoom.grob, width = 6, height = 5)
#   
#   KA.full <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(ylim = c(0.0015, 0.0390)) +
#     labs(y = 'Normalized reflectance', x = 'Wavelength (nm)', title = bquote(italic("K. angustifolia"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     theme_bw() +
#     fulltheme
#   KA.zoom <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(xlim = c(475, 675), ylim = c(0.0037, 0.0121)) +
#     labs(title = bquote(bold("VIS"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     zoomtheme
#   KA.zoom.grob <- ggplotGrob(KA.zoom)
#   KA.zoom.grob <- KA.full + 
#     annotation_custom(grob = KA.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.023, ymax = 0.0399)
#   # ggsave('figures/spectra_KA_site.png', KA.zoom.grob, width = 6, height = 5)
#   
#   RG.full <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(ylim = c(0.0015, 0.0390)) +
#     labs(y = '', x = 'Wavelength (nm)', title = bquote(italic("R. groenlandicum"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     theme_bw() +
#     fulltheme
#   RG.zoom <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value)) +
#     geom_line(aes(color = site_ordered, group = site_ordered)) +
#     coord_cartesian(xlim = c(475, 675), ylim = c(0.0025, 0.0098)) +
#     labs(title = bquote(bold("VIS"))) +
#     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
#     scale_colour_manual(name = "Nitrogen fertilization level", values = brewer.pal(8, "YlOrRd")[c(3, 6, 8)]) +
#     zoomtheme
#   RG.zoom.grob <- ggplotGrob(RG.zoom)
#   RG.zoom.grob <- RG.full + 
#     annotation_custom(grob = RG.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.0245, ymax = 0.0399)
#   # ggsave('figures/spectra_RG_site.png', RG.zoom.grob, width = 6, height = 5)
# }
# 
# tous <- ggarrange(CC.zoom.grob, EV.zoom.grob, KA.zoom.grob, RG.zoom.grob,
#                   ncol = 2, nrow = 2, common.legend = T)
# # ggsave('figures/spectra_4sp_site.png', tous, width = 8, height = 8)


