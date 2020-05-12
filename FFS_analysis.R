# 24 mars 2020

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#            Data set for band-wise/FFS analysis on NORMALIZED spectra with plot        #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : Col 1 is rank (#1 is the most different band and so on), Col 2 is the **FFS metric (higher number = greater separability)**

library(dplyr)
library(tidyr) # gather
library(tibble) # rownames_to_column()
library(readtext) # read txtEdit files
library(ggplot2)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Create sharable data ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
setwd("~/Documents/Maîtrise/DonnéesAnalyses/spectres")
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_large.RData") # rt_sg_large
spectra_large_provisoire <- rt_sg_large

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData") # rt_sg_normalized_large
normalized_spectra_large_provisoire <- rt_sg_normalized_large

# add plot#
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/tous_s_id.RData") # tous_s_id

raw_bulk <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./tele_bulk_leaf_samples.csv", sep = ";") %>% # changer eriophorum à la main !
  select(c("plant_id","sample_id")) %>% 
  filter(sample_id %in% tous_s_id$sample_id)
raw_bulk$sample_id <- as.character(raw_bulk$sample_id) 

raw_plants <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./raw_plants.csv", sep = ";") %>% 
  select("plant_id", "plot_id")

spectra_large <- full_join(spectra_large_provisoire, raw_bulk, by = "sample_id") %>% 
  full_join(raw_plants) %>%
  select("sample_id","scientific_name","plot_id","site","site_class","propriete",
         "400":"2400", -"plant_id") %>% # réordonner et enlever impertinent
  filter(sample_id %in% tous_s_id$sample_id)
   ## vérifier si duplication de sample_id
   length(unique(spectra_large$sample_id))
   nrow(spectra_large) # 2* plus c chill (transmittance et réflectance)

normalized_spectra_large <- full_join(normalized_spectra_large_provisoire, raw_bulk, by = "sample_id") %>% 
  full_join(raw_plants) %>%
  select("sample_id","scientific_name","plot_id","site","site_class","propriete","normalization_magnitude",
         "400":"2400", -"plant_id") %>% # réordonner et enlever impertinent
  filter(sample_id %in% tous_s_id$sample_id)
   ## vérifier si duplication de sample_id
   length(unique(normalized_spectra_large$sample_id))
   nrow(normalized_spectra_large) # 2* plus c chill (transmittance et réflectance)

# large to long
spectra_long <- gather(spectra_large, wavelength, value, -c("sample_id":"propriete"))
normalized_spectra_long <- gather(normalized_spectra_large, wavelength, value, -c("sample_id":"normalization_magnitude"))

# save in csv
spectra_large.2 <- as.data.frame(spectra_large)
# write.table(spectra_large.2, file = "~/Desktop/spectra_large.csv")

spectra_long.2 <- as.data.frame(spectra_long)
# write.table(spectra_long.2, file = "~/Desktop/spectra_long.csv")

normalized_spectra_large.2 <- as.data.frame(normalized_spectra_large)
# write.table(normalized_spectra_large.2, file = "~/Desktop/normalized_spectra_large.csv")

normalized_spectra_long.2 <- as.data.frame(normalized_spectra_long)
# write.table(normalized_spectra_long.2, file = "~/Desktop/normalized_spectra_long.csv")

# visualize (bonnes données?)
library("ggplot2")
RG <- spectra_long.2 %>% 
  dplyr::filter(scientific_name == "Rhododendron groenlandicum") %>%
  dplyr::filter(propriete == "reflectance") 
RG$wavelength <- as.numeric(RG$wavelength)
ggplot(RG) + geom_line(aes(x = wavelength, y = value)) 
# ok, pas d'artefact

RG <- normalized_spectra_long.2 %>% 
  dplyr::filter(scientific_name == "Rhododendron groenlandicum") %>%
  dplyr::filter(propriete == "reflectance") 
RG$wavelength <- as.numeric(RG$wavelength)
ggplot(RG) + geom_line(aes(x = wavelength, y = value))


#-=-=-=-=-=-=-=-=-=-=-=-=-
# FSS results ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
                 
## Load FSS data and create NN criterion graph
setwd("~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis") 

## Open FFS data + adjust columns
RG <- read.table("RGnorm.txt", sep = ",", row.names = NULL)
RG$sc <- "R. groenlandicum" # Add scientific name column
colnames(RG) <- c("rank", "FFS.metric", "band.number", "sc")
# first global maxima in FFS metric = which rank (equals line number) ?
RG_highlight <- filter(RG, RG$rank <= which.max(RG$FFS.metric)) # rank = 4 (rank is line number)
rg.plot <- ggplot(RG) + geom_line(aes(rank, FFS.metric)) + geom_vline(xintercept = which.max(RG$FFS.metric) , color = "blue") + geom_hline(yintercept = max(RG$FFS.metric), color = "blue") + geom_text(aes(label = paste0("rank = ", as.character(which.max(RG$FFS.metric)), ", FFS = ", max(RG$FFS.metric))), x = 700, y = 0.6) + labs(title = bquote(italic("Rhododendron groenlandicum"))) + theme_bw()
ggsave("~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/NNcriterion_RG.png", rg.plot, width = 4, height = 3)

EV <- read.table("evnorm.txt", sep = ",", row.names = NULL)
EV$sc <- "E. vaginatum"  # Add scientific name column
colnames(EV) <- c("rank", "FFS.metric", "band.number", "sc")
# first global maxima in FFS metric = which rank (equals line number) ?
EV_highlight <- filter(EV, EV$rank <= which.max(EV$FFS.metric)) # rank = 4 (rank is line number)
ev.plot <- ggplot(EV) + geom_line(aes(rank, FFS.metric)) + geom_vline(xintercept = which.max(EV$FFS.metric) , color = "blue") + geom_hline(yintercept = max(EV$FFS.metric), color = "blue") + geom_text(aes(label = paste0("rank = ", as.character(which.max(EV$FFS.metric)), ", FFS = ", max(EV$FFS.metric))), x = 500, y = 0.45) + labs(title = bquote(italic("Eriophorum vaginatum"))) + theme_bw()
ggsave("~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/NNcriterion_EV.png", ev.plot, width = 4, height = 3)

KA <- read.table("kanorm.txt", sep = ",", row.names = NULL)
KA$sc <- "K. angustifolia" # Add scientific name column
colnames(KA) <- c("rank", "FFS.metric", "band.number", "sc")
# first global maxima in FFS metric = which rank (equals line number) ?
KA_highlight <- filter(KA, KA$rank <= which.max(KA$FFS.metric)) # rank = 4 (rank is line number)
ka.plot <- ggplot(KA) + geom_line(aes(rank, FFS.metric)) + geom_vline(xintercept = which.max(KA$FFS.metric) , color = "blue") + geom_hline(yintercept = max(KA$FFS.metric), color = "blue") + geom_text(aes(label = paste0("rank = ", as.character(which.max(KA$FFS.metric)), ", FFS = ", max(KA$FFS.metric))), x = 500, y = 0.4) + labs(title = bquote(italic("Kalmia angustifolia"))) + theme_bw()
ggsave("~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/NNcriterion_KA.png", ka.plot, width = 4, height = 3)

CC <- read.table("ccnorm.txt", sep = ",", row.names = NULL)
CC$sc <- "C. calyculata" # Add scientific name column
colnames(CC) <- c("rank", "FFS.metric", "band.number", "sc")
# first global maxima in FFS metric = which rank (equals line number) ?
CC_highlight <- filter(CC, CC$rank <= which.max(CC$FFS.metric)) # rank = 4 (rank is line number)
cc.plot <- ggplot(CC) + geom_line(aes(rank, FFS.metric)) + geom_vline(xintercept = which.max(CC$FFS.metric) , color = "blue") + geom_hline(yintercept = max(CC$FFS.metric), color = "blue") + geom_text(aes(label = paste0("rank = ", as.character(which.max(CC$FFS.metric)), ", FFS = ", max(CC$FFS.metric))), x = 500, y = 0.3) + labs(title = bquote(italic("Chamaedaphne calyculata"))) + theme_bw()
ggsave("~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/NNcriterion_CC.png", cc.plot, width = 4, height = 3)


## Combine all FFS results 
FFS.res.long <- rbind(RG_highlight, EV_highlight, KA_highlight, CC_highlight)

# verify : add 399 or 400 to get back WVL ?
RG$wvl399 <- RG$band.number + 399 # add 399 so first wvl is 400 and not 401
RG$wvl400 <- RG$band.number + 400

FFS.res.long$wvl <- FFS.res.long$band.number + 399
# write.csv(FFS.res.long, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/FFS_important_bands.csv", row.names = FALSE)


## Load and adjust spectral data for graphing
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData") # rt_sg_normalized_large
refl_norm_large <- rt_sg_normalized_large %>% dplyr::filter(propriete == "reflectance") %>%  # from the "spectre_normalisation.R" script 
  select(-c("site", "propriete"))
refl_norm_long <- gather(refl_norm_large, wavelength, value, -c("sample_id":"normalization_magnitude")) # from large to long format
refl_norm_long$wavelength <- as.numeric(refl_norm_long$wavelength) # adjust column type

# visualize (make sure I have the right data)
RG_spectra <- refl_norm_long %>% dplyr::filter(scientific_name == "Rhododendron groenlandicum") 
RG_spectra$wavelength <- as.numeric(RG_spectra$wavelength)
ggplot(RG_spectra) + geom_line(aes(x = wavelength, y = value)) + theme_bw() # ok, normalized and no artefacts

# adjust species and site naming
refl_norm_long$sn_graph = factor(refl_norm_long$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                            "Kalmia angustifolia", 'Rhododendron groenlandicum'),
                      labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum"))

refl_norm_long$site_graph = factor(refl_norm_long$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open", "MBP_Intermediate", "MBP_High"), 
                                    labels = c("GTV", "GPB", "MBP.No.N", "MBP.Mid.N", "MBP.High.N"))


#### FFS results graphics ####
library("RColorBrewer")
library(grid) # unit()
library(gridExtra) # grid.arrange()

## Species mean spectra
{refl_mean_site_RG <- refl_norm_long %>% 
  dplyr::filter(sn_graph == "R. groenlandicum") %>%
  group_by(wavelength, site_graph) %>% 
  summarize(mean_value = mean(value), 
            STD_value = sd(value)) %>% 
  ungroup()

refl_mean_site_EV <- refl_norm_long %>% 
  dplyr::filter(sn_graph == "E. vaginatum") %>%
  group_by(wavelength, site_graph) %>% 
  summarize(mean_value = mean(value), 
            STD_value = sd(value)) %>% 
  ungroup()

refl_mean_site_KA <- refl_norm_long %>% 
  dplyr::filter(sn_graph == "K. angustifolia") %>%
  group_by(wavelength, site_graph) %>% 
  summarize(mean_value = mean(value), 
            STD_value = sd(value)) %>% 
  ungroup()

refl_mean_site_CC <- refl_norm_long %>% 
  dplyr::filter(sn_graph == "C. calyculata") %>%
  group_by(wavelength, site_graph) %>% 
  summarize(mean_value = mean(value), 
            STD_value = sd(value)) %>% 
  ungroup()}


## Themes : one for the entire graph (full), one for the VIS inlet (zoom)
fulltheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), plot.title = element_text(hjust = 0),
                   legend.position = "none", plot.margin = margin(0, 3.5, 0, 0, "mm"))

zoomtheme <- theme(axis.line = element_blank(),axis.text.x=element_blank(), legend.position="none", 
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),axis.title.y=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(color = 'black', fill = "white"),
                   plot.margin = margin(0, 0, 0, -2, "mm")) #  plot.title = element_text(hjust = 0.5),


## RG
RGnorm_FFSrankings.full <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) + # color = brewer.pal(9, "Greys")[7], 
  geom_ribbon(aes(fill = site_graph, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.04), expand = F, clip = "off") +
  labs(y = "Normalized reflectance", x = 'Wavelength (nm)', title = bquote(italic("Rhododendron groenlandicum"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(data = RG_highlight, xintercept = (RG_highlight$band.number + 399), alpha = 0.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
RG.zoom <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) +
  coord_cartesian(xlim = c(510, 620), ylim = c(0.001, 0.012), expand = F) +
  annotate("text", label = 'bold("VIS")',  x = 565, y = 0.005, parse = T) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
RG.zoom.grob <- ggplotGrob(RG.zoom)
RG.zoom.grob <- RGnorm_FFSrankings.full +
  annotation_custom(grob = RG.zoom.grob, xmin = 1800, xmax = 2400, ymin = 0.0225, ymax = 0.04)
RG.zoom.grob
# ggsave('~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/RGnorm_FFS_VIS.png', RG.zoom.grob, width = 4, height = 3)
# ggsave('~/Documents/Maîtrise/Rédaction/Figures(.zip)/Fig.4.FFS-bott.right.png', RG.zoom.grob, width = 4, height = 3)


## EV
EVnorm_FFSrankings.full <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) + # color = brewer.pal(9, "Greys")[7], 
  geom_ribbon(aes(fill = site_graph, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.04), expand = F, clip = "off") +
  labs(y = "Normalized reflectance", x = 'Wavelength (nm)', title = bquote(italic("Eriophorum vaginatum"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(data = EV_highlight, xintercept = (EV_highlight$band.number + 399), alpha = 0.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
EV.zoom <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) +
  coord_cartesian(xlim = c(510, 620), ylim = c(0.001, 0.012), expand = F) +
  annotate("text", label = 'bold("VIS")',  x = 565, y = 0.005, parse = T) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
EV.zoom.grob <- ggplotGrob(EV.zoom)
EV.zoom.grob <- EVnorm_FFSrankings.full +
  annotation_custom(grob = EV.zoom.grob, xmin = 1800, xmax = 2400, ymin = 0.0235, ymax = 0.04)
EV.zoom.grob
# ggsave('~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/EVnorm_FFS_VIS.png', EV.zoom.grob, width = 4, height = 3)
# ggsave('~/Documents/Maîtrise/Rédaction/Figures(.zip)/Fig.4.FFS-topright.png', EV.zoom.grob, width = 4, height = 3)


## KA
KAnorm_FFSrankings.full <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) + # color = brewer.pal(9, "Greys")[7], 
  geom_ribbon(aes(fill = site_graph, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.04), expand = F, clip = "off") +
  labs(y = "Normalized reflectance", x = 'Wavelength (nm)', title = bquote(italic("Kalmia angustifolia"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(data = KA_highlight, xintercept = (KA_highlight$band.number + 399), alpha  = 0.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
KA.zoom <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) +
  coord_cartesian(xlim = c(510, 620), ylim = c(0.001, 0.0155), expand = F) +
  annotate("text", label = 'bold("VIS")',  x = 565, y = 0.005, parse = T) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
KA.zoom.grob <- ggplotGrob(KA.zoom)
KA.zoom.grob <- KAnorm_FFSrankings.full +
  annotation_custom(grob = KA.zoom.grob, xmin = 1800, xmax = 2400, ymin = 0.0225, ymax = 0.04)
KA.zoom.grob
# ggsave('~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/KAnorm_FFS_VIS.png', KA.zoom.grob, width = 4, height = 3)
# ggsave('~/Documents/Maîtrise/Rédaction/Figures(.zip)/Fig.4.FFS-bott.left.png', KA.zoom.grob, width = 4, height = 3)


## CC
CCnorm_FFSrankings.full <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) + # color = brewer.pal(9, "Greys")[7], 
  geom_ribbon(aes(fill = site_graph, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.04), expand = F, clip = "off") +
  labs(y = "Normalized reflectance", x = 'Wavelength (nm)', title = bquote(italic("Chamaedaphne calyculata"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(data = CC_highlight, xintercept = (CC_highlight$band.number + 399), alpha = 0.5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
CC.zoom <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value, color = site_graph)) +
  geom_line(size = 0.7, aes(group = site_graph, linetype = site_graph)) +
  coord_cartesian(xlim = c(510, 620), ylim = c(0.001, 0.012), expand = F) +
  annotate("text", label = 'bold("VIS")',  x = 565, y = 0.005, parse = T) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
CC.zoom.grob <- ggplotGrob(CC.zoom)
CC.zoom.grob <- CCnorm_FFSrankings.full +
  annotation_custom(grob = CC.zoom.grob, xmin = 1800, xmax = 2400, ymin = 0.0235, ymax = 0.04)
CC.zoom.grob
# ggsave('~/Documents/Maîtrise/DonnéesAnalyses/spectres/FFS_analysis/Figures/CCnorm_FFS_VIS.png', CC.zoom.grob, width = 4, height = 3)
# ggsave('~/Documents/Maîtrise/Rédaction/Figures(.zip)/Fig.4.FFS-topleft.png', CC.zoom.grob, width = 4, height = 3)

