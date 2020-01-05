# 27 août

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                              Visualisation et NRMSE, PLSR cute
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : sainte !


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
library(dplyr)
library(ggplot2)
library("RColorBrewer") # adjust species color
library(stringr)
library(grid)
# library(data.table)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# obtenu via nettoyage toutes données
load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData") # ft_all_vX2
ft_all <- ft_all_vX2 # version vX2 -> enlevé les nég dans recalcitrant (version précédante : "ARCHIVEft_all")
ft_all$sample_id <- as.numeric(ft_all$sample_id)

# obtenus via section Manip. préalables dans ce script
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/pred_meas/df_93.RData")
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/pred_meas/df_94.RData")
NRMSE_dt <- read.csv("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/NRMSE.csv")


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Manipulations préalables ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# importer toutes les données de prédictions dans l'environnement
setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/pred_meas")
temp = list.files(pattern = "*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))),
         read.csv), envir = .GlobalEnv)
ft_pred.meas_prov <- mget(ls(envir = .GlobalEnv, pattern = "preds"))
ft_pred.meas <- list()

inVars_prov <- c(gsub("comps_abs_preds", "", names(ft_pred.meas_prov)))
inVars_prov2 <- c(gsub(".{2}$", "", inVars_prov))
inVars <- c(gsub("_$", "", inVars_prov2))

# for each, joindre l'espèce et le NRSME
for (i in 1:length(ft_pred.meas_prov)) {
  inVar <- inVars[i]
  print(c(i, inVar))
  ft_pred.meas[[i]] <- ft_pred.meas_prov[[i]][, grep(inVar, colnames(ft_pred.meas_prov[[i]]))] # ne garder que la colonne qui correspond exactement à "inVar" pour chaque DF dans la liste
  ft_pred.meas[[i]][,"sample_id"] <- ft_pred.meas_prov[[i]][, "sample_id"]
  ft_pred.meas[[i]] <- left_join(ft_pred.meas[[i]], ft_all[, c("scientific_name", "sample_id")])
  ft_pred.meas[[i]][, "scientific_name"] = factor(ft_pred.meas[[i]][, "scientific_name"], levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                                    "Kalmia angustifolia", 'Rhododendron groenlandicum'), 
                           labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum")) 
  
  colnames(ft_pred.meas[[i]])[c(2,3)] <- c("predMean", "predStdv")
}
names(ft_pred.meas) <- inVars

# importer toutes les statistiques des modèles dans l'environnement
# setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/calcul_NRMSE")
# 
# temp2 = list.files(pattern = "*.csv")
# list2env(
#   lapply(setNames(temp2, make.names(gsub("*.csv$", "", temp2))),
#          read.csv), envir = .GlobalEnv)
# ft_stats <- mget(ls(envir = .GlobalEnv, pattern = "stats"))
# 
# inVars_prov <- c(gsub("[0-9]*", "", names(ft_pred.meas)))
# inVars <- c(gsub("_comps_abs_preds", "", inVars_prov))
NRMSE_dt <- read.csv("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/NRMSE.csv")
colnames(NRMSE_dt) <- c("Functionnal.trait", "Spectrum.range", "Ncomps","R2.cal", "R2.val",
                        "RMSE.cal", "RMSE.val", "model.NRMSE", "Range", "Mean.sd") 

# vérifier l'ordre
bla <- as.data.frame(NRMSE_dt[,1])
bla$deux <- names(ft_pred.meas)
bla$trois <- as.vector(inVars)


#-=-=-=-=-=-=-=-=-=-=-=-=-
#  Affichage du graphique de prédiction ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R")

# créer un df avec les position de chaque label colnames(annotate 1x, annotate1y, annotate 2x, annotate 2y) 
df_pos <- as.data.frame(matrix(data = NA, nrow = length(ft_pred.meas), ncol = 6), stringsAsFactors = FALSE)
df_pos[,1] <- inVars 
  #              1   2    3       4   5    6    7    8     9    10    11   12  13    14   15   16   17   18
df_pos[,2] <- c(43, 55, 2.2,    25, 5.7, 800, 2.1, 300, 0.024, 30,   518, 154,645, 20.2, 1.88, 1,  70, 12.2)
df_pos[,3] <- c(24, 47, 0.3,    5.8,  1.5, 160, 0.9, 70,  0.009, 6,    325, 94, 492, 2.5,  1.2, 0.01,40, 6.4)
df_pos[,4] <- c(43, 55, 2.2,    25, 5.7, 800, 2.1, 300, 0.024, 30,   518, 154, 645,20.2, 1.88, 1,  70, 12.2)
df_pos[,5] <- c(26.7,46, 0.475, 8.5, 2,  235, 1.07, 95, 0.011, 9.7, 349., 102, 460, 5,   1.3,0.13,45.3, 5.5)
df_pos[,6] <- c("8", "5", "6", "7", "6", "6", "6", "6", "6", "10", "5", "6", "5", "7", "8", "5", "11", "6")
colnames(df_pos) <- c("TF", "annotate1x", "annotate1y", "annotate1y", "annotate2y", "lbl.comps")


## NOMS des FT (doivent fiter)
quatre <- as.list(c('C:N ratio', 'Total C (%)', bquote("Carotenoids (mg g"^-1~")"), 'Cellulose (%)', 
                    bquote("Chl"~italic(a)~"(mg g"^-1~")"), bquote("Chl"~italic(a)~"(mg m"^-2~")"), 
                    bquote("Chl"~italic(b)~"(mg g"^-1~")"), bquote("Chl"~italic(b)~"(mg m"^-2~")"),
                    'EWT (cm)', 'Hemicellulose (%)', bquote("LDMC (mg g"^-1~")"), bquote("LMA (g m"^-2~")"),
                    bquote("LWC (mg g"^-1~")"), 'Lignin (%)', 'N (%)', 'recalcitrant', 'Soluble carbon (%)',
                    bquote("SLA (m"^2~"kg"^-1~")")))

# vérificationde l'ordre des TF
bla.quatre <- as.vector(c('C:N ratio', 'Total C (%)', "Carotenoids (mg g)", 'Cellulose (%)', 
                "Chl(mg g)", "Chl(mg m)", 
                "Chl(mg g)", "Chl(mg m)",
                'EWT (cm)', 'Hemicellulose (%)', "LDMC (mg g)", "LMA (g m)",
                "LWC (mg g)", 'Lignin (%)', 'N (%)', 'recalcitrant', 'Soluble carbon (%)',
                "SLA (mkg)"))
bla$quatre <- bla.quatre

# generate the graphs
for (i in 1:length(ft_pred.meas)) {
  inVar <- inVars[i]
  print(c(i, inVar))
  lbl.coefs <- paste0('atop(italic(R)^2 ==', round(mean(NRMSE_dt[i, "R2.val"]), 2),
                      ',\nRMSE == ', round(mean(NRMSE_dt[i, "RMSE.val"]), 2),')')
  # ',\nNRMSE ==', NRMSE_dt[i, "model.NRMSE"], 
  
  bla <- ggplot(ft_pred.meas[[i]], aes(x = ft_pred.meas[[i]][, inVar], y = ft_pred.meas[[i]][,"predMean"], color = scientific_name)) +
    geom_point(aes(shape = scientific_name), size = 2) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(method = "lm", se = F) +
    geom_smooth(method = "lm", se = F, inherit.aes = F, aes(x = ft_pred.meas[[i]][, inVar], y = ft_pred.meas[[i]][,"predMean"]), 
                color = "black", linetype = "dashed", size = 0.6) +
    # coord_cartesian(expand = F) +
    annotate("label", label = lbl.coefs, parse = TRUE, label.size = 0, x = df_pos[i,2], y = df_pos[i,3]) +
    annotate("label", label = paste("N comps :",  df_pos[i, 6]), label.size = 0, x = df_pos[i,4], y = df_pos[i,5]) +
    geom_errorbar(data = ft_pred.meas[[i]], 
                  aes(x = ft_pred.meas[[i]][, inVar],
                      ymin = ft_pred.meas[[i]][, "predMean"] - ft_pred.meas[[i]][, "predStdv"],
                      ymax = ft_pred.meas[[i]][, "predMean"] + ft_pred.meas[[i]][, "predStdv"],
                      color = scientific_name)) +
    scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
    scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
    labs(y = 'Predicted', x = 'Measured', title = quatre[[i]]) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0), legend.text = element_text(face = "italic")) #, legend.position = "none",
  bla
    # ggsave(bla , file = paste0("graphs_finaux/", inVar, ".jpg"), width = 4, height = 4.3)
    ggsave(bla , file = "graphs_finaux/LÉGENDE.jpg", width = 4, height = 4.3)
  
}

