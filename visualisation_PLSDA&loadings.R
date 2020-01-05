# 29 août

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                        Visualisation PLSDA et loadings                                #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : 

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)
library("RColorBrewer") # coloblnd friendly colors
library(tidyverse) # column_to_rownames
library(ggcorrplot)

setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/DA")


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Matrice de confusion  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

tabs_perc <- read.csv("R_output/PLSDA_confuperc_especes_4comps.csv")
# obtenue via manips sur "PLS-DA.R"

tabs_perc <- column_to_rownames(tabs_perc, var = "X") 

labs <- c(expression(italic('C. calyculata')), expression(italic('E. vaginatum')),
  expression(italic("K. angustifolia")), expression(italic('R. groenlandicum')))

colors <- c("#66A61E", "#969696")

corr <- ggcorrplot(tabs_perc[,], method = "circle", lab = T, digits = 1, tl.cex = 10, 
                   lab_size = 2.5, outline.color = "#969696", tl.col = "#969696") +
  scale_fill_gradient2(low = colors[2], high = colors[1], midpoint = 50, 
                            limit = c(0, 100), name = "Classification\naccuracy (%)") + 
  scale_y_discrete(labels = labs) +
  scale_x_discrete(labels = labs) +
  theme(legend.title = element_text(size = 9))
corr
# ggplot2::ggsave("figures_finales/confusion_70_40_4ncomps.pdf", corr, width = 4, height = 4)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Loadings ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

ldngs <- read.csv(file = "4ncomps_ldngs.csv", as.is = T)[-1,]

# n_comps <- 3 # vérifier dans le nom !
ldngs$wvl <- as.numeric(as.character(ldngs$wvl)) 

# trouver les wvl où se trouvent les pics
hgst_ldngs <- ldngs[order(-ldngs$mean),][1:100,]

range1 <- seq(713, 735)
hgst_ldngs1 <- filter(hgst_ldngs, wvl %in% range1)
hgst_ldngs1_hgst <- head(sort(hgst_ldngs1$mean, decreasing = T), 1)
which(hgst_ldngs$mean == hgst_ldngs1) # ligne 323
ifelse(hgst_ldngs$mean == hgst_ldngs1_hgst, hgst_ldngs$wvl, NA)  # wvl 722

range2 <- seq(1410, 1437)
hgst_ldngs2 <- filter(hgst_ldngs, wvl %in% range2)
hgst_ldngs2_hgst <- head(sort(hgst_ldngs2$mean, decreasing = T), 1)
ifelse(hgst_ldngs$mean == hgst_ldngs2_hgst, hgst_ldngs$wvl, NA)  # wvl 1419

range3 <- seq(1881, 1888)
hgst_ldngs3 <- filter(hgst_ldngs, wvl %in% range3)
hgst_ldngs3_hgst <- head(sort(hgst_ldngs3$mean, decreasing = T), 1)
ifelse(hgst_ldngs$mean == hgst_ldngs3_hgst, hgst_ldngs$wvl, NA)  # wvl 1881

range4 <- seq(2000, 2400)
hgst_ldngs4 <- filter(hgst_ldngs, wvl %in% range4)
hgst_ldngs4_hgst <- head(sort(hgst_ldngs4$mean, decreasing = T), 1)
ifelse(hgst_ldngs$mean == hgst_ldngs4_hgst, hgst_ldngs$wvl, NA)  # wvl 1881

peaks <- c(722, 1419, 1881, 2162)

bla.2 <- ggplot(ldngs, aes(x = wvl, y = mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin = mean - (sd * 1.96),
                  ymax = mean + (sd * 1.96)), 
              alpha = 0.45, fill = "grey") +
  labs(x = "Wavelength (nm)", y = "Absolute loading") +
  coord_cartesian(ylim = c(0, 0.02), expand = F) +
  theme_light() +
  geom_vline(xintercept = peaks, alpha = 0.6, color = "#FB6A4A") +
  annotate("text", label = peaks, x = peaks, y = 0.0025, angle = 90, size = 3) +
  theme(panel.border = element_rect(colour = brewer.pal(9, "Greys")[5], fill = NA, size = 1)) #,  
        # panel.grid = element_blank())
bla.2
# ggsave('/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/DA/figures_finales/loadings_4comps_7 nov.png', width = 6, height = 4)
# dev.off()

