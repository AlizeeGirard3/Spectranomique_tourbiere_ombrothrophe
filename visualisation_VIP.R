# 8 oct 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                             # Tests pour VIP et Loadings                             #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)
library("RColorBrewer") # coloblnd friendly colors
library(tidyverse) # column_to_rownames
library(grid)
library(plotly)
library("varhandle") # unfactor


#-=-=-=-=-=-=-=-=-=-=-=-=-
# VIP (PLS-R)  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

setwd("~/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/vip")

# importer toutes les données de VIP dans l'environnement
temp = list.files(pattern = "*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))),
         read.csv), envir = .GlobalEnv)

# faire une liste des tableaux de données
ft_VIP <- mget(ls(envir = .GlobalEnv, pattern = "vipAggr"))

# faire des noms cutes pour les graphs ATTENTION CE DOIT ÊTRE LE MM ORDRE QUE names(FT_VIP)
noms <- as.list(c('C:N ratio', 'Total C (%)', bquote("Carotenoids (mg g"^-1~")"), 'Cellulose (%)', 
                   bquote("Chl"~italic(a)~"(mg g"^-1~")"), bquote("Chl"~italic(a)~"(mg m"^-2~")"), 
                   bquote("Chl"~italic(b)~"(mg g"^-1~")"), bquote("Chl"~italic(b)~"(mg m"^-2~")"),
                   'EWT (cm)', 'Hemicellulose (%)', bquote("LDMC (mg g"^-1~")"), bquote("LMA (g m"^-2~")"),
                   bquote("LWC (mg g"^-1~")"), 'Lignin (%)', 'N (%)', 'recalcitrant', 'Soluble carbon (%)',
                   bquote("SLA (m"^2~"kg"^-1~")")))

  # vérification de l'ordre des noms 
bla <- as.data.frame(names(ft_VIP))
bla$deux <- noms
  # oui, c'est dans le même ordre

# pour trouver les wvl associés aux pics
# i<-
# plot_ly(data = ft_VIP[[i]], x = ~ band, y = ~ mean_VIP, type = "pointcloud") %>%
#   layout(title = names(ft_VIP)[i]) 

# df avec pics
pics <- matrix(data = NA, nrow = nrow(bla), ncol = 9)
pics[, 1] <- names(ft_VIP)
{
  pics[1,2:8]  <- c(400, 454, 555, 714, 1922, 2312, 2400)
  pics[2,2:5] <- c(1233, 1926, 2309, 2400)
  pics[3,2:6] <- c(400, 477, 533, 681, 760)
  pics[4,2:7] <- c(1254, 1667, 1925, 2158, 2320, 2400)
  pics[5,2:8] <- c(400, 425, 451, 510, 556, 689, 760)
  pics[6,2:7] <- c(425, 510, 553, 677, 688, 760)
  pics[7,2:8] <- c(400, 427, 451, 509, 555, 677, 760)
  pics[8,2:5] <- c(400, 553, 677, 760)
  pics[9,2:5] <- c(1405, 1714, 1929, 2400)
  pics[10,2:9] <- c(1200, 1395, 1666, 1881, 1927, 2158, 2313, 2400)
  pics[11,2:5] <- c(1407, 1713, 1927, 2400)
  pics[12,2:6] <- c(700, 711, 1115, 1211, 1350)
  pics[13,2:5] <- c(1407, 1713, 1927, 2400)
  pics[14,2:7] <- c(1236, 1395, 1716, 1924, 2293, 2398)
  pics[15,2:9] <- c(400, 454, 552, 715, 1718, 1922, 2313, 2398)
  pics[17,2:7] <- c(1273, 1503, 1667, 1925, 2157, 2400)
  pics[18,2:6] <- c(700, 712, 962, 1212, 1350)
}
pics <- as.data.frame(pics, stringsAsFactors = F, col.names = F) %>% 
  mutate_at(funs = is.character, vars(c(2:9)), as.numeric)


for (i in 1:length(ft_VIP)) {
  print(names(ft_VIP)[i])
  
  bla <- ggplot(ft_VIP[[i]], aes(x = band, y = mean_VIP)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean_VIP - (stdv * 1.96),
                    ymax = mean_VIP + (stdv * 1.96)),
                alpha = 0.45, fill = "grey") +
    labs(x = "Wavelength (nm)", y = "VIP", title = noms[[i]]) +
    geom_vline(xintercept = as.vector(pics[i, 2:9], mode = "numeric"), color = "#FB6A4A", na.rm = T, alpha = 0.65) +
    annotate("text", label = as.vector(pics[i, 2:9], mode = "character"), 
             x = as.vector(pics[i, 2:9], mode = "numeric"), y = c(1.6,1.4,1.6,1.4,1.6,1.4,1.6,1.4), angle = 90, size = 3) +
    coord_cartesian(clip = "off") + # xlim = c(min(ft_VIP[[i]]$band)-5, max(ft_VIP[[i]]$band)+5), expand = F,
    theme_bw() +
    theme(plot.title = element_text(hjust = 0))
  # ggsave(bla , file = paste0("~/Documents/Maîtrise/DonnéesAnalyses/PLS/R/graphs_finaux/",
  #                            names(ft_VIP)[i], "_VIP.jpg"), width = 3.3, height = 3.5)
  ggsave(bla , file = paste0("~/Desktop/", names(ft_VIP)[i], "_VIP.jpg"), width = 3.3, height = 3.5)
}

