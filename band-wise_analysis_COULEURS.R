# 14 oct 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#            Band-wise analysis on NORMALIZED spectra per sites per species             #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)
library(tidyr) # gather
library("RColorBrewer")
library(ggpubr) # for ggarrange
library(tibble) # rownames_to_column()


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

setwd("~/Documents/Maîtrise/DonnéesAnalyses/spectres")

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData") # rt_sg_normalized_large
refl_norm_large <- rt_sg_normalized_large %>% 
  dplyr::filter(propriete == "reflectance") %>% 
  dplyr::select(-"sample_id", -"propriete", -"site", -"normalization_magnitude")
# issu du script "spectre_normalisation.R" dans document scripts-NETTOYAGE


# prevent scientific notation and characters to be converted to factors
options(scipen = F, stringsAsFactors = F)


# -=-=-=-=-=-=-=-=-=-=-=-=-
# Calcul sur les bandes ----
# -=-=-=-=-=-=-=-=-=-=-=-=-

# séparer les espèces
EV <- refl_norm_large %>%
  dplyr::filter(scientific_name == "Eriophorum vaginatum")
KA <- refl_norm_large %>%
  dplyr::filter(scientific_name == "Kalmia angustifolia")
CC <- refl_norm_large %>%
  dplyr::filter(scientific_name == "Chamaedaphne calyculata")
RG <- refl_norm_large %>%
  dplyr::filter(scientific_name == "Rhododendron groenlandicum")

# définir la fonction
aov_func_sites <- function(y) aov(y ~ site_class, data = dada)

# évaluer les différences/wvl/sp
{
  ## EV
  dada <- EV
  aov.res <- list()
  inVars <- colnames(dada[, 3:2002])
  aov.res.df <- matrix(data = NA, nrow = length(inVars), ncol = 2) # colnames : wvl, p-val lm

  for (var in seq(inVars)) {
    print(var)
    vars <- inVars[var]
    aov.res[[var]] <- aov_func_sites(dada[, vars])
    aov.res.df[var, 1] <- paste0("X", var + 399)
    aov.res.df[var, 2] <- round(as.numeric(paste(summary(aov.res[[var]])[[1]][["Pr(>F)"]][[1]])), 4) # extract p-value from aov.res
  }
  aov.res.df <- as.data.frame(aov.res.df, stringsAsFactors = F)
  colnames(aov.res.df)[1:2] <- c("wavelength", "p.value")
  EV_abs_spec_stats <- aov.res.df %>%
  mutate(`corr.p-valueBonFer` = p.adjust(p.value, method = "bonferroni")) %>%
  mutate(`corr.p-value.benjHochberg` = p.adjust(p.value, method = "BH"))
  # write.csv(EV_abs_spec_stats,"BAND-WISE_ANALYSIS/EV_aov_per_wvl_abs_reflectance.csv", row.names = F)

  ## KA
  dada <- KA
  aov.res <- list()
  inVars <- colnames(dada[, 3:2002])
  aov.res.df <- matrix(data = NA, nrow = length(inVars), ncol = 2) # colnames : wvl, p-val lm

  for (var in seq(inVars)) {
    print(var)
    vars <- inVars[var]
    aov.res[[var]] <- aov_func_sites(dada[, vars])
    aov.res.df[var, 1] <- paste0("X", var + 399)
    aov.res.df[var, 2] <- round(as.numeric(paste(summary(aov.res[[var]])[[1]][["Pr(>F)"]][[1]])), 4) # extract p-value from aov.res
  }
  aov.res.df <- as.data.frame(aov.res.df, stringsAsFactors = F)
  colnames(aov.res.df)[1:2] <- c("wavelength", "p.value")
  KA_abs_spec_stats <- aov.res.df %>%
    mutate(`corr.p-valueBonFer` = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(`corr.p-value.benjHochberg` = p.adjust(p.value, method = "BH"))
  # write.csv(KA_abs_spec_stats,"BAND-WISE_ANALYSIS/KA_aov_per_wvl_abs_reflectance.csv", row.names = F)

  ## CC
  dada <- CC
  aov.res <- list()
  inVars <- colnames(dada[, 3:2002])
  aov.res.df <- matrix(data = NA, nrow = length(inVars), ncol = 2) # colnames : wvl, p-val lm

  for (var in seq(inVars)) {
    print(var)
    vars <- inVars[var]
    aov.res[[var]] <- aov_func_sites(dada[, vars])
    aov.res.df[var, 1] <- paste0("X", var + 399)
    aov.res.df[var, 2] <- round(as.numeric(paste(summary(aov.res[[var]])[[1]][["Pr(>F)"]][[1]])), 4) # extract p-value from aov.res
  }
  aov.res.df <- as.data.frame(aov.res.df, stringsAsFactors = F)
  colnames(aov.res.df)[1:2] <- c("wavelength", "p.value")
  CC_abs_spec_stats <- aov.res.df %>%
    mutate(`corr.p-valueBonFer` = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(`corr.p-value.benjHochberg` = p.adjust(p.value, method = "BH"))
  # write.csv(CC_abs_spec_stats,"BAND-WISE_ANALYSIS/CC_aov_per_wvl_abs_reflectance.csv", row.names = F)

## RG
  dada <- RG
  aov.res <- list()
  inVars <- colnames(dada[, 3:2002])
  aov.res.df <- matrix(data = NA, nrow = length(inVars), ncol = 2) # colnames : wvl, p-val lm

  for (var in seq(inVars)) {
    print(var)
    vars <- inVars[var]
    aov.res[[var]] <- aov_func_sites(dada[, vars])
    aov.res.df[var, 1] <- paste0("X", var + 399)
    aov.res.df[var, 2] <- round(as.numeric(paste(summary(aov.res[[var]])[[1]][["Pr(>F)"]][[1]])), 4) # extract p-value from aov.res
  }
  aov.res.df <- as.data.frame(aov.res.df, stringsAsFactors = F)
  colnames(aov.res.df)[1:2] <- c("wavelength", "p.value")
  RG_abs_spec_stats <- aov.res.df %>%
    mutate(`corr.p-valueBonFer` = p.adjust(p.value, method = "bonferroni")) %>%
    mutate(`corr.p-value.benjHochberg` = p.adjust(p.value, method = "BH"))
  # write.csv(RG_abs_spec_stats,"BAND-WISE_ANALYSIS/RG_aov_per_wvl_abs_reflectance.csv", row.names = F)

}


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Affichage par espèce p.val BENJ HOCHBERG ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
{
  refl_long <- refl_norm_large %>%
    gather(key = "wavelength", value = "valeur_abs_sg", c("400":"2400"), -scientific_name, -site_class)
  refl_long$wavelength <- as.numeric(paste(refl_long$wavelength))
  
  # # noms pour graphiques
  refl_long$sc = factor(refl_long$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                              "Kalmia angustifolia", 'Rhododendron groenlandicum'),
                        labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum"))
  
  refl_long$site_class = factor(refl_long$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open",  "MBP_Intermediate", "MBP_High"), 
                                labels = c("GTV", "GPB", "MBP.No.N", "MBP.Mid.N", "MBP.High.N"))
  
  ## EV, moyenne
  refl_mean_site_EV <- refl_long %>% 
    dplyr::filter(sc == "E. vaginatum") %>%
    group_by(site_class, wavelength) %>%
    summarize(mean_value = mean(valeur_abs_sg), 
              STD_value = sd(valeur_abs_sg)) %>% 
    ungroup()
  
  # obtenir les wvl significatifs
  EV_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/EV_aov_per_wvl_abs_reflectance.csv")
  sign_wvl_prov <- as.vector(ifelse(EV_abs_spec_stats$corr.p.value.benjHochberg < 0.05, EV_abs_spec_stats$wavelength, NA))
  sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  
  # affichage
  EV_abs_sign.bands <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value, color = site_class)) +
    geom_line(size = 0.5, aes(linetype = site_class)) + # color = brewer.pal(9, "Greys")[7], 
    geom_ribbon(aes(fill = site_class, 
                    ymin = mean_value - STD_value,
                    ymax = mean_value + STD_value, 
                    color = NA),
                alpha = 0.2) +
    coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
    labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Eriophorum vaginatum"))) +
    scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
    scale_linetype_manual(name = "Site",values = c("dotdash", "solid", "solid", "solid", "solid")) +
    scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    geom_vline(xintercept = sign_wvl, alpha = .015) +
    theme_bw() +
    theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank(), legend.position = "none") 
  # EV_abs_sign.bands
  # ggsave('BAND-WISE_ANALYSIS/BenjHochEV_abs_sign.bands_COULEURS_LÉGENDE.png', EV_abs_sign.bands, width = 4, height = 3)
  ggsave('BAND-WISE_ANALYSIS/BenjHochEV_abs_sign.bands_COULEURS.png', EV_abs_sign.bands, width = 4, height = 3)
  
  ## KA, moyenne
  refl_mean_site_KA <- refl_long %>% 
    dplyr::filter(sc == "K. angustifolia") %>%
    group_by(site_class, wavelength) %>%
    summarize(mean_value = mean(valeur_abs_sg), 
              STD_value = sd(valeur_abs_sg)) %>% 
    ungroup()
  
  # obtenir les wvl significatifs
  KA_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/KA_aov_per_wvl_abs_reflectance.csv")
  sign_wvl_prov <- as.vector(ifelse(KA_abs_spec_stats$corr.p.value.benjHochberg < 0.05, KA_abs_spec_stats$wavelength, NA))
  sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  
  # affichage
  KA_abs_sign.bands <-  ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value, color = site_class)) +
    geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) + # color = brewer.pal(9, "Greys")[7], 
    geom_ribbon(aes(fill = site_class, 
                    ymin = mean_value - STD_value,
                    ymax = mean_value + STD_value, 
                    color = NA),
                alpha = 0.2) +
    coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
    labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Kalmia angustifolia"))) +
    scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
    scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
    # scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
    geom_vline(xintercept = sign_wvl, alpha = .015) +
    theme_bw() +
    theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank(), legend.position = "none")
  # KA_abs_sign.bands
  ggsave('BAND-WISE_ANALYSIS/BenjHochKA_abs_sign_COULEURS.bands.png', KA_abs_sign.bands, width = 4, height = 3)
  
  ## CC, moyenne
  refl_mean_site_CC <- refl_long %>% 
    dplyr::filter(sc == "C. calyculata") %>%
    group_by(site_class, wavelength) %>%
    summarize(mean_value = mean(valeur_abs_sg), 
              STD_value = sd(valeur_abs_sg)) %>% 
    ungroup()
  
  # obtenir les wvl significatifs
  CC_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/CC_aov_per_wvl_abs_reflectance.csv")
  sign_wvl_prov <- as.vector(ifelse(CC_abs_spec_stats$corr.p.value.benjHochberg < 0.05, CC_abs_spec_stats$wavelength, NA))
  sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  
  # affichage
  CC_abs_sign.bands <-  ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value, color = site_class)) +
    geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) + # color = brewer.pal(9, "Greys")[7], 
    geom_ribbon(aes(fill = site_class, 
                    ymin = mean_value - STD_value,
                    ymax = mean_value + STD_value, 
                    color = NA),
                alpha = 0.2) +
    coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
    labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Chamaedaphne calyculata"))) +
    scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
    scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
    # scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
    geom_vline(xintercept = sign_wvl, alpha = .015) +
    theme_bw() +
    theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank(), legend.position = "none")
  # CC_abs_sign.bands
  ggsave('BAND-WISE_ANALYSIS/BenjHochCC_abs_sign.bands_COULEURS.png', CC_abs_sign.bands, width = 4, height = 3)
  
  ## RG, moyenne
  refl_mean_site_RG <- refl_long %>% 
    dplyr::filter(sc == "R. groenlandicum") %>%
    group_by(site_class, wavelength) %>%
    summarize(mean_value = mean(valeur_abs_sg), 
              STD_value = sd(valeur_abs_sg)) 
  
  # obtenir les wvl significatifs
  RG_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/RG_aov_per_wvl_abs_reflectance.csv")
  sign_wvl_prov <- as.vector(ifelse(RG_abs_spec_stats$corr.p.value.benjHochberg < 0.05, RG_abs_spec_stats$wavelength, NA))
  sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  
  # affichage
  RG_abs_sign.bands <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value, color = site_class)) +
    geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) + # color = brewer.pal(9, "Greys")[7], 
    geom_ribbon(aes(fill = site_class, 
                    ymin = mean_value - STD_value,
                    ymax = mean_value + STD_value, 
                    color = NA),
                alpha = 0.2) +
    coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
    labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Rhododendron groenlandicum"))) +
    scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
    scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
    scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
    # scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
    geom_vline(xintercept = sign_wvl, alpha = .015) +
    theme_bw() +
    theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) #, legend.position = "none")
  # RG_abs_sign.bands
  ggsave('BAND-WISE_ANALYSIS/Légende_COULEURS.png', RG_abs_sign.bands, width = 4, height = 3)
}

# tous <- ggarrange(EV_abs_sign.bands, KA_abs_sign.bands, CC_abs_sign.bands, RG_abs_sign.bands,
#                   ncol = 2, nrow = 2, legend = "right", common.legend = T)
# tous <- annotate_figure(tous, fig.lab = "p-value corrigée avec Benjamini-Hochberg", fig.lab.size = 10)
# # tous
# ggsave('BAND-WISE_ANALYSIS/BenjHoch_ABS_sign.bands_spectra_4sp_site.png', tous, width = 8, height = 8)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Affichage par espèce p.val BENJ HOCHBERG + VIS zoomé ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

  # obtenir objets et noms pour graphiques
refl_long <- refl_large %>%
  gather(key = "wavelength", value = "valeur", c("400":"2400"), -scientific_name, -site_class)
refl_long$wavelength <- as.numeric(paste(refl_long$wavelength))

# # noms pour graphiques
refl_long$sc = factor(refl_long$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                            "Kalmia angustifolia", 'Rhododendron groenlandicum'),
                      labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum"))

refl_long$site_class = factor(refl_long$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open",  "MBP_Intermediate", "MBP_High"), 
                              labels = c("GTV", "GPB", "MBP.No.N", "MBP.Mid.N", "MBP.High.N"))

# définir les boîtes zoomées
library(grid) # for unit
library(gridExtra) # for grid.arrange

fulltheme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), plot.title = element_text(hjust = 0),
                   legend.position = "none", plot.margin = margin(0, 0, 0, 0, "mm"))

zoomtheme <- theme(legend.position="none", axis.line = element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),axis.title.y=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(color = 'black', fill = "white"),
                   plot.title = element_text(hjust = 0.5), plot.margin = margin(2, 2, 0, 0, "mm"))

# créer objets
refl_long_mean_site_EV <- refl_long %>%
  dplyr::filter(sc == "E. vaginatum") %>%
  group_by(site_class, wavelength) %>%
  summarize(mean_value = mean(valeur),
            STD_value = sd(valeur))

refl_long_mean_site_CC <- refl_long %>%
  dplyr::filter(sc == "C. calyculata") %>%
  group_by(site_class, wavelength) %>%
  summarize(mean_value = mean(valeur),
            STD_value = sd(valeur))

refl_long_mean_site_KA <- refl_long %>%
  dplyr::filter(sc == "K. angustifolia") %>%
  group_by(site_class, wavelength) %>%
  summarize(mean_value = mean(valeur),
            STD_value = sd(valeur))

refl_long_mean_site_RG <- refl_long %>%
  dplyr::filter(sc == "R. groenlandicum") %>%
  group_by(site_class, wavelength) %>%
  summarize(mean_value = mean(valeur),
            STD_value = sd(valeur))


## CC
# obtenir les wvl significatifs
CC_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/CC_aov_per_wvl_abs_reflectance.csv")
sign_wvl_prov <- as.vector(ifelse(CC_abs_spec_stats$corr.p.value.benjHochberg < 0.05, CC_abs_spec_stats$wavelength, NA))
sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]

CC.full <- ggplot(refl_long_mean_site_CC, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Chamaedaphne calyculata"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(xintercept = sign_wvl, alpha = .015) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
CC.zoom <- ggplot(refl_long_mean_site_CC, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class,
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value,
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(470, 650), ylim = c(0.06, 0.2)) +
  labs(title = bquote(bold("VIS"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
CC.zoom.grob <- ggplotGrob(CC.zoom)
CC.zoom.grob <- CC.full +
  annotation_custom(grob = CC.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.42, ymax = 0.7)
# ggsave('figures/spectra_CC_VIS.png', CC.zoom.grob, width = 4, height = 3)

#   ## EV
# obtenir les wvl significatifs
EV_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/EV_aov_per_wvl_abs_reflectance.csv")
sign_wvl_prov <- as.vector(ifelse(EV_abs_spec_stats$corr.p.value.benjHochberg < 0.05, EV_abs_spec_stats$wavelength, NA))
sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]

EV.full <- ggplot(refl_long_mean_site_EV, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Eriophorum vaginatum"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(xintercept = sign_wvl, alpha = .015) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
EV.zoom <- ggplot(refl_long_mean_site_EV, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class,
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value,
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(470, 650), ylim = c(0.04, 0.18)) +
  labs(title = bquote(bold("VIS"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
EV.zoom.grob <- ggplotGrob(EV.zoom)
EV.zoom.grob <- EV.full +
  annotation_custom(grob = EV.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.42, ymax = 0.7)
# ggsave('figures/spectra_EV_VIS.png', EV.zoom.grob, width = 4, height = 3)


#   ## KA
# obtenir les wvl significatifs
KA_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/KA_aov_per_wvl_abs_reflectance.csv")
sign_wvl_prov <- as.vector(ifelse(KA_abs_spec_stats$corr.p.value.benjHochberg < 0.05, KA_abs_spec_stats$wavelength, NA))
sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]

KA.full <- ggplot(refl_long_mean_site_KA, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Kalmia angustifolia"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(xintercept = sign_wvl, alpha = .015) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
KA.zoom <- ggplot(refl_long_mean_site_KA, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class,
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value,
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(470, 650), ylim = c(0.04, 0.18)) +
  labs(title = bquote(bold("VIS"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
KA.zoom.grob <- ggplotGrob(KA.zoom)
KA.zoom.grob <- KA.full +
  annotation_custom(grob = KA.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.42, ymax = 0.7)
# ggsave('figures/spectra_KA_VIS.png', KA.zoom.grob, width = 4, height = 3)


#   ## RG
# obtenir les wvl significatifs
RG_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/RG_aov_per_wvl_abs_reflectance.csv")
sign_wvl_prov <- as.vector(ifelse(RG_abs_spec_stats$corr.p.value.benjHochberg < 0.05, RG_abs_spec_stats$wavelength, NA))
sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]

RG.full <- ggplot(refl_long_mean_site_RG, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class, 
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value, 
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.8), expand = F, clip = "off") +
  labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Rhododendron groenlandicum"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  geom_vline(xintercept = sign_wvl, alpha = .015) +
  theme_bw() +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), panel.grid = element_blank()) + #, legend.position = "none")
  fulltheme
RG.zoom <- ggplot(refl_long_mean_site_RG, aes(x = wavelength, y = mean_value, color = site_class)) +
  geom_line(size = 0.5, aes(group = site_class, linetype = site_class)) +
  geom_ribbon(aes(fill = site_class,
                  ymin = mean_value - STD_value,
                  ymax = mean_value + STD_value,
                  color = NA),
              alpha = 0.2) +
  coord_cartesian(xlim = c(470, 650), ylim = c(0.04, 0.18)) +
  labs(title = bquote(bold("VIS"))) +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_fill_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
  scale_linetype_manual(name = "Site", values = c("dotdash", "solid", "solid", "solid", "solid")) +
  zoomtheme
RG.zoom.grob <- ggplotGrob(RG.zoom)
RG.zoom.grob <- RG.full +
  annotation_custom(grob = RG.zoom.grob, xmin = 1500, xmax = 2400, ymin = 0.46, ymax = 0.8)
# ggsave('figures/spectra_RG_VIS.png', RG.zoom.grob, width = 4, height = 3)


  
  ############################# ARCHIVES ############################# 
  
  # #-=-=-=-=-=-=-=-=-=-=-=-=-
  # # Affichage par espèce p.val non corrigée ----
  # #-=-=-=-=-=-=-=-=-=-=-=-=-
  # {
  #   refl_long <- refl_large %>%
  #     gather(key = "wavelength", value = "valeur_abs_sg", c("400":"2400"), -scientific_name, -site_class)
  #   refl_long$wavelength <- as.numeric(paste(refl_long$wavelength))
  #   
  #   # noms pour graphiques
  #   refl_long$sc = factor(refl_long$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
  #                                                                         "Kalmia angustifolia", 'Rhododendron groenlandicum'), 
  #                              labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum")) 
  #   
  #   refl_long$site_class = factor(refl_long$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open",  "MBP_Intermediate", "MBP_High"), 
  #                                      labels = c("GTV.Background", "GPB.Background", "MBP.Background", "MBP.Intermediate", "MBP.High")) 
  #   
  #   ## EV, moyenne
  #   refl_mean_site_EV <- refl_long %>% 
  #     dplyr::filter(sc == "E. vaginatum") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   EV_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/EV_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(EV_abs_spec_stats$p.value < 0.05, EV_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   EV_abs_sign.bands <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class, 
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Eriophorum vaginatum"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw() +
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # EV_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/EV_abs_sign.bands.png', EV_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## KA, moyenne
  #   refl_mean_site_KA <- refl_long %>% 
  #     dplyr::filter(sc == "K. angustifolia") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   KA_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/KA_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(KA_abs_spec_stats$p.value < 0.05, KA_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   KA_abs_sign.bands <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class, 
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Kalmia angustifolia"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # KA_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/KA_abs_sign.bands.png', KA_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## CC, moyenne
  #   refl_mean_site_CC <- refl_long %>% 
  #     dplyr::filter(sc == "C. calyculata") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   CC_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/CC_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(CC_abs_spec_stats$p.value < 0.05, CC_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   CC_abs_sign.bands <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class, 
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Chamaedaphne calyculata"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # CC_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/CC_abs_sign.bands.png', CC_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## RG, moyenne
  #   refl_mean_site_RG <- refl_long %>% 
  #     dplyr::filter(sc == "R. groenlandicum") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) 
  #   
  #   # obtenir les wvl significatifs
  #   RG_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/RG_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(RG_abs_spec_stats$p.value < 0.05, RG_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   RG_abs_sign.bands <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class, 
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Rhododendron groenlandicum"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # RG_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/RG_abs_sign.bands.png', RG_abs_sign.bands, width = 4, height = 3)
  #   
  # }
  # 
  # tous <- ggarrange(EV_abs_sign.bands, KA_abs_sign.bands, CC_abs_sign.bands, RG_abs_sign.bands,
  #                   ncol = 2, nrow = 2)
  # tous <- annotate_figure(tous, fig.lab = "p-value de l'anova", fig.lab.size = 10)
  # # tous
  # ggsave('BAND-WISE_ANALYSIS/ABS_sign.bands_spectra_4sp_site.png', tous, width = 8, height = 8)
  # 
  # 
  # #-=-=-=-=-=-=-=-=-=-=-=-=-
  # # Affichage par espèce p.val BON FERRONI ----
  # #-=-=-=-=-=-=-=-=-=-=-=-=-
  # {
  #   refl_long <- refl_large %>%
  #     gather(key = "wavelength", value = "valeur_abs_sg", c("400":"2400"), -scientific_name, -site_class)
  #   refl_long$wavelength <- as.numeric(paste(refl_long$wavelength))
  #   
  #   # noms pour graphiques
  #   refl_long$sc = factor(refl_long$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
  #                                                                         "Kalmia angustifolia", 'Rhododendron groenlandicum'), 
  #                              labels = c("C. calyculata","E. vaginatum","Kalmia angustifolia", "R. groenlandicum")) 
  #   
  #   refl_long$site_class = factor(refl_long$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open",  "MBP_Intermediate", "MBP_High"), 
  #                                      labels = c("GTV.Background", "GPB.Background", "MBP.Background", "MBP.Intermediate", "MBP.High")) 
  #   
  #   ## EV, moyenne
  #   refl_mean_site_EV <- refl_long %>% 
  #     dplyr::filter(sc == "E. vaginatum") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   EV_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/EV_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(EV_abs_spec_stats$corr.p.valueBonFer < 0.05, EV_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   EV_abs_sign.bands <- ggplot(refl_mean_site_EV, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class, 
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +  
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Eriophorum vaginatum"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw() +
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # EV_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/BonFerrEV_abs_sign.bands.png', EV_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## KA, moyenne
  #   refl_mean_site_KA <- refl_long %>% 
  #     dplyr::filter(sc == "K. angustifolia") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   KA_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/KA_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(KA_abs_spec_stats$corr.p.value.benjHochberg < 0.05, KA_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   KA_abs_sign.bands <- ggplot(refl_mean_site_KA, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class,
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Kalmia angustifolia"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # KA_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/BonFerrKA_abs_sign.bands.png', KA_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## CC, moyenne
  #   refl_mean_site_CC <- refl_long %>% 
  #     dplyr::filter(sc == "C. calyculata") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) %>% 
  #     ungroup()
  #   
  #   # obtenir les wvl significatifs
  #   CC_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/CC_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(CC_abs_spec_stats$corr.p.value.benjHochberg < 0.05, CC_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   CC_abs_sign.bands <- ggplot(refl_mean_site_CC, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class,
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Chamaedaphne calyculata"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # CC_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/BonFerrCC_abs_sign.bands.png', CC_abs_sign.bands, width = 4, height = 3)
  #   
  #   ## RG, moyenne
  #   refl_mean_site_RG <- refl_long %>% 
  #     dplyr::filter(sc == "R. groenlandicum") %>%
  #     group_by(site_class, wavelength) %>%
  #     summarize(mean_value = mean(valeur_abs_sg), 
  #               STD_value = sd(valeur_abs_sg)) 
  #   
  #   # obtenir les wvl significatifs
  #   RG_abs_spec_stats <- read.csv("BAND-WISE_ANALYSIS/RG_aov_per_wvl_abs_reflectance.csv")
  #   sign_wvl_prov <- as.vector(ifelse(RG_abs_spec_stats$corr.p.value.benjHochberg < 0.05, RG_abs_spec_stats$wavelength, NA))
  #   sign_wvl_prov.2 <- as.numeric(gsub("X", "", sign_wvl_prov))
  #   sign_wvl <- sign_wvl_prov.2[!is.na(sign_wvl_prov.2)]
  #   
  #   # affichage
  #   RG_abs_sign.bands <- ggplot(refl_mean_site_RG, aes(x = wavelength, y = mean_value)) +
  #     geom_line(color = brewer.pal(9, "Greys")[7], size = 0.26, aes(group = site_class)) +
  #     geom_ribbon(aes(fill = site_class,
  #                     ymin = mean_value - STD_value,
  #                     ymax = mean_value + STD_value),
  #                 alpha = 0.3) +
  #     coord_cartesian(xlim = c(400, 2400), ylim = c(0.001, 0.7), expand = F, clip = "off") +
  #     labs(y = "Reflectance", x = 'Wavelength (nm)', title = bquote(italic("Rhododendron groenlandicum"))) +
  #     scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2400)) +
  #     # scale_colour_manual(name = "", values = c(brewer.pal(9, "Greys")[c(4,9)],"#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     scale_fill_manual(name = "", values = c(rep(brewer.pal(9, "Greys")[7], times = 5))) + # , "#FEB24C", "#FD8D3C", "#FC4E2A")) +
  #     geom_vline(xintercept = sign_wvl, alpha = .04) +
  #     theme_bw()+
  #     theme(plot.margin = unit(c(0.75,0.75,0.75,0.75), "lines"), legend.position = "none")
  #   # RG_abs_sign.bands
  #   # ggsave('BAND-WISE_ANALYSIS/BonFerrRG_abs_sign.bands.png', RG_abs_sign.bands, width = 4, height = 3)
  #   
  # }
  # 
  # tous <- ggarrange(EV_abs_sign.bands, KA_abs_sign.bands, CC_abs_sign.bands, RG_abs_sign.bands,
  #                   ncol = 2, nrow = 2) # , common.legend = T)
  # tous <- annotate_figure(tous, fig.lab = "p-value corrigées avec Bon Ferroni", fig.lab.size = 10)
  # # # tous
  # ggsave('BAND-WISE_ANALYSIS/BonFerr_ABS_sign.bands_spectra_4sp_site.png', tous, width = 8, height = 8)
  