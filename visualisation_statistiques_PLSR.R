# 14 oct

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                              Visualisation et NRMSE
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# obtenu via nettoyage toutes données
load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData") # ft_all_vX2
ft_all <- ft_all_vX2 %>% # version vX2 -> enlevé les négatifs de recalcitrants, anyway, on s'en fou on s'en sert pas !
  mutate_at(funs = is.character, vars(c(9:29)), as.numeric) %>% 
  dplyr::select(-c("Phenols_mg_g", "Tannins_mg_g"))
  
# obtenus via section Manip. préalables dans ce script
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/tous_sites_21aout_16sept/donnees_pour_graphs/df_93.RData")
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/tous_sites_21aout_16sept/donnees_pour_graphs/df_94.RData")


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Manipulations préalables ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# absent de quelques tables : les outliers !
# setdiff(chlA_chlB_6comps_preds$sample_id,equivalent_water_thickness_cm_5comps_preds$sample_id)

# importer toutes les données de prédictions dans l'environnement
setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/pred_meas")
temp = list.files(pattern = "*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))),
         read.csv), envir = .GlobalEnv)

# importer toutes les statistiques des modèles dans l'environnement
setwd("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/resultats/calcul_NRMSE")

temp2 = list.files(pattern = "*.csv")
list2env(
  lapply(setNames(temp2, make.names(gsub("*.csv$", "", temp2))),
         read.csv), envir = .GlobalEnv)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Calcul NRMSE ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

## obtenir données 
inVars <- sort(names(ft_all)[c(8:23, 26:27)]) # ordre alphabetique
rmse_measured <- list()
ft_names_prov <- gsub("\\_[0-9]*", "", temp2)
ft_names <- gsub("compsabsstats.csv", "", ft_names_prov)
NRMSE_dt <- data.frame(ft = ft_names, matrix(NA, nrow = length(ft_names), ncol = 7), stringsAsFactors = F)
List = lapply(temp2, read.csv)

# vérification ** IL FAUT QUE CE SOIENT LES MÊMES NOMS, DANS LE MÊME ORDRE !
bla <- as.data.frame(inVars)
bla$bladeux <- ft_names

## créer df avec valeurs d'intérêt
for (i in 1:length(ft_names)) {
  inVar <- inVars[i]
  print(c(i, inVar)) 
  NRMSE_dt[i,4] <- round(as.numeric(mean(List[[i]][, 1])), 2) # "fitR2"
  NRMSE_dt[i,5] <- round(as.numeric(mean(List[[i]][, 4])), 2) #  "valR2"
  NRMSE_dt[i,6] <- round(as.numeric(mean(List[[i]][, 2])), 2) # "fitRMSE"
  
  # NRMSE_dt[i,6] <- ifelse(i == 9, round(as.numeric(mean(List[[i]][, 2])), 5), round(as.numeric(mean(List[[i]][, 2])), 2)) # pour EWT, arrondir à 5 décimales
  
  NRMSE_dt[i,7] <- round(as.numeric(mean(List[[i]][, 5])), 2) # "valRMSE"
  NRMSE_dt[i,8] <- round(as.numeric(mean(List[[i]][, 5], na.rm = T) / mean(as.numeric(ft_all[, inVar]), na.rm = T)) * 100, 2)
  NRMSE_dt[i,9] <- paste(round(as.numeric(range(ft_all[, inVar], na.rm = T)[1]), 2), 
                         round(as.numeric(range(ft_all[, inVar], na.rm = T)[2]), 2), sep = "-")
  NRMSE_dt[i, 10] <- paste(round(mean(ft_all[, inVar], na.rm = T), 2), " (", round(sd(ft_all[, inVar], na.rm = T), 2), ")", sep = "")
}

# EWT 5 décimales
(c(round(as.numeric(mean(List[[9]][, 2])), 5),
  round(as.numeric(mean(List[[9]][, 5])), 5),
  round(as.numeric(mean(List[[9]][, 5], na.rm = T) / mean(as.numeric(ft_all[, "equivalent_water_thickness_cm"]), na.rm = T)) * 100, 5)))
# 0.00178  0.00245 16.66499

# ajouter colonnes "spectrum range" & "ncomps"
NRMSE_dt[,2] <- as.vector(c("400-2400", "1200-2400", "400-760", "1200-2400", "400-760", "400-760","400-760", "400-760",
                            "800-2400", "1200-2400", "800-2400", "800-2400", "800-2400", "1200-2400", "400-2400", "1200-2400", "1200-2400", 
                            "800-2400"))
NRMSE_dt[,3] <- as.vector(c("8", "5", "6", "7", "6", "6", "6", "6", "6", "10", "5", "6", "5", "7", "8", "5", "11", "6"))# solC 11, c% 7comps, carot 6, 

# ajuster les noms de colones
colnames(NRMSE_dt) <- c("Functionnal trait", "Spectrum range (nm)", "Number of components",  "R2 cal", "R2 val",
                        "RMSE* cal", "RMSE val", "model NRMSE (%)", "Range", "Mean (sd)") # écrire en nbp que RMSE sont en unités originales

# write.csv(NRMSE_dt, file = "/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/R/NRMSE.csv", row.names = FALSE)

