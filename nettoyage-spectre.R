# 6 MARS

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                                  NETTOYAGE SPECTRES
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : 
# données issues de l'algorithme : http://data.caboscience.org/field-data/projects/2018-Girard-MSc-UdeM/spectra/ 
# ou par feuilles : http://data.caboscience.org/field-data/projects/2018-Girard-MSc-UdeM/spectra/processed/project_leaves_combined.csv      

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Bibliothèques
#-=-=-=-=-=-=-=-=-=-=-=-=-

#install.packages("")

library(tidyverse) 
library(plyr)
library(dplyr)
library(ggplot2)
library(magrittr)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# APPEL DATAFRAME "rt"
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# PLUS RÉCENTES VERSIONS (overlap retiré manuellement)
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_all_no_ovrlp.RData") # rt_all_no_ovrlp (tous wvl) / PRÉALABLEMENT rt_v6_no_overlap, changé 18 juin (mtn ARCHIVErt_v6_no_overlp)
# **********et/ou**********
load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long_no_ovrlp.RData") # rt_long_no_ovrlp (long, 400:2400)
# tous deux issu de  manips dans ce script sur all_combined (initialement; issu directement du web), dernière version du calcul = 15 mai 2019


# VERSIONS RÉCENTES
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_v6_large.RData") # rt_v6_large (tous wvl)
# rt_v6_large issu de  manips dans ce script sur all_combined (issu directement du web), dernière version du calcul (15 mai 2019)
load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long.RData") # rt_long (long, 400:2400)


# VIELLES VERSION
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/ref_trans_v5.RData") # rt_v5
# rt_v5 est issu de manips dans ce script sur all_combined (issu directement du web), dernière version du calcul (8 avril 2019)
# calculs du 8 avril on laissé 2 problèmes corrigés le 12 avril (13837104 et 12086180)

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/ref_trans_v3.RData") # rt_v3
# rt_v3 est issu de manips dans ce script sur rt_v2

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/ref_trans.RData") # rt_v2
# rt_v2 est issu du code "analyses_spectres2.R" et des scripts qui ont mené à lui ("analyses_spectres.R" et "archives_nouvelles_donnees.R", etc)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# CODE DU DATAFRAME "rt"
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# NOTES : 
# travailler avec ordi puissant pour obtenir "rt" !

# project_all_combined -> les sample_id par moyenne des 6 feuilles
all_combined <- read.csv('http://data.caboscience.org/field-data/projects/2018-Girard-MSc-UdeM/spectra/processed/project_all_combined.csv')

tous_spectres_v4 <- all_combined %>%
  dplyr::select(-c(record_id, measured_by, spectroradiometer_start_time, spectroradiometer_id, instrumentation_id, leaf_side_measured))

tous_spectres_v4$sample_id <- as.character(tous_spectres_v4$sample_id)
# cols = c("sample_id")
# tous_spectres_v4[,cols] %<>% lapply(function(x) as.character(x))

# simplifier noms d'espèces
unique(tous_spectres_v4$scientific_name)
tous_spectres_v4$scientific_name <- revalue(tous_spectres_v4$scientific_name, # revalue() => pkg "plyr"
                                            c("Eriophorum vaginatum subsp. spissum (Fernald) Hultén" = "Eriophorum vaginatum",
                                              "Eriophorum vaginatum Linnaeus subsp. vaginatum" = "Eriophorum vaginatum",
                                              'Kalmia angustifolia Linnaeus' = 'Kalmia angustifolia',
                                              'Kalmia angustifolia Linnaeus var. angustifolia' = 'Kalmia angustifolia',
                                              'Rhododendron groenlandicum (Oeder) Kron & Judd' = 'Rhododendron groenlandicum',
                                              'Chamaedaphne calyculata (Linnaeus) Moench' = 'Chamaedaphne calyculata'))

# ajouter site
raw_bulk <- read.csv("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./tele_bulk_leaf_samples.csv", sep = ";")
site <- raw_bulk %>%
  dplyr::filter(scientific_name != "Sphagnum capillifolium (Ehrh.) Hedw.") %>%
  dplyr::filter(! site_id %in% c("MBP_temporal", "Labo-B234")) %>%
  dplyr::filter(! status == "deleted") %>%
  dplyr::filter(! sample_id %in% c("9236658", "9249843", "9320621", "9324461", "9335280", "9338745", "9342418", "9668337", "9671728", "9674948")) %>%
  dplyr::select(sample_id, site_id) # ne pas mettre scientfic names ! si j'en ai besoin, faire un autre objet
site$site_id = droplevels(site$site_id)
site$sample_id <- as.character(site$sample_id)

# joindre le dataframe site avec le dataframe tous_spectres_vX / ne pas changer v4 pour v3, sert à rien
tous_spectres_v7 <- left_join(tous_spectres_v4, site, by = 'sample_id') # left_join enlèvera les sample_id qui sont dans site et pas dans tous_spectres
#save(tous_spectres_v7, file = '/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/ref_trans_v7.RData') # v5 faite le 10 mai 2019 / # v6 faite le 15 mai, tout corrigé / v7 fait le 20 juin



# Reflectance et transmittance en 2 colonnes
refl_v6 <- tous_spectres_v7 %>% # refl_v4 fait le 10 mai pour vérifier nouveaux calculs # refl_v5 fait 15 mai, tout fini / v6 20 juin
  dplyr::filter(reflectance_transmittance =='reflectance') %>%
  dplyr::select(R_T_Average, R_T_STD, sample_id, scientific_name, wavelength, site_id)
colnames(refl_v6) <- c("R.Ave", "R.STD", "sample_id", "scientific_name", "wavelength", "site")

trans_v6 <- tous_spectres_v7 %>%
  dplyr::filter(reflectance_transmittance =='transmittance') %>%
  dplyr::select(R_T_Average, R_T_STD, sample_id, scientific_name, wavelength, site_id)
colnames(trans_v6) <- c("T.Ave", "T.STD", "sample_id", "scientific_name", "wavelength",  "site")
trans_v6[,'T.Ave'] <- 1-trans_v6[,'T.Ave'] # fait 15 mai

# joindre le dataframe reflectance avec le dataframe transmittance
# .
# ...
# .....
# .......   ATTENTION !
# ......... ne pas changer pour rt_v3 ! ou changer ci-bas rt_v2 -> rt_v3 et rt_v3 -> rt_v4 faire avec crtl+f et replace !! ****
# .......   changer pour rt_v5 et + / 10 mai
# .....
# ...
# .
rt_v7_large <- left_join(refl_v6, trans_v6) # rt_v5 pour vérifications, avec refl_v4 et trans_v4, 10 mai / rt_v6 fait avec refl_v5 et trans_v5 15 mai

#save(rt_v7_large, file = '/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_v7_large.RData') # fait 8 mars / v5 fait 10 mai / v6 15 mai, fini corriger / 21 juin OFFICIEL
#load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_v7_large.RData")


#=======================================# RETIRER LES CHEVAUCHEMENTS #=======================================#                    ----

# modifié d'Etienne Github : warren-spectra/scripts/data_manipulations/Step3_spectra_refl_calculations.R)

# *****sur RT_VX_LARGE ****

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Détermination des longueurs d'ondes dupliquées ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

s_id_13837104 <- rt_v7_large %>%
  dplyr::filter(sample_id == "13837104") %>%
  #dplyr::select(-c("T.Ave", "T.STD")) %>%
  dplyr::select(-c("R.Ave", "R.STD"))

which(duplicated(s_id_13837104$wavelength) == T) # √ LIGNES 633-804 -> wvl 970- 1013 // √ LIGNES 1687-1754 -> wvl 1895-1911

# faire objet avec wvls
all_wvls <- dplyr::filter(rt_v7_large, sample_id == "13837104")$wavelength

# créer objet avec wvl à enlever
all_wvls[634:805] # enlever
rem1 <- all_wvls[513:523]
all_wvls[1688:1754] # enlever
rem2 <- all_wvls[769:773]
rem <- c(rem1, rem2)

# enlever wvls
rt_all_no_ovrlp_provisoire <- rt_v7_large %>%
  dplyr::filter(!wavelength %in% rem)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Interpolation linéaire sur l'objet ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# interpolation linéaire aux 1 nm pour tous les spectres
inter_wvls <- ceiling(min(rt_all_no_ovrlp_provisoire$wavelength)):floor(max(rt_all_no_ovrlp_provisoire$wavelength))
interpolate <- function(x) {
  wvls <- x$wavelength
  refl <- x$R.Ave
  trans <- x$T.Ave
  new_refl <- approx(wvls, refl, xout = inter_wvls)$y  # Return a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation.
  new_trans <- approx(wvls, trans, xout = inter_wvls)$y
  tmp <- data_frame(wvl = inter_wvls, refl = new_refl, trans = new_trans)
  return(tmp)
}

# appliquer la fonction
rt_all_no_ovrlp_provisoire2 <- rt_all_no_ovrlp_provisoire %>%
  group_by(sample_id) %>%
  do(interpolate(.)) %>%
  ungroup() # nécessaire parce que le group_by et ensuite le join by sur sample_id (ci-bas) faisait bogger mon R et mon ordi !!!

# vérification ----
cinq_echantillons_INTERPOL <- rt_all_no_ovrlp_provisoire2 %>%
  dplyr::filter(sample_id %in% c("13383299", "11995935", "13224081", "14009324", "11929171"))

cinq_echantillons_ORIGINAL <- rt_v7_large %>%
  dplyr::filter(sample_id %in% c("13383299", "11995935", "13224081", "14009324", "11929171"))

rem_ovrlp_verif <- ggplot() +
  geom_line(data = cinq_echantillons_INTERPOL, aes(x = wvl, y = refl,
                                                   group = sample_id)) +
  geom_line(data = cinq_echantillons_INTERPOL,aes(x = wvl, y = trans,
                                                  group = sample_id)) +
  geom_line(data = cinq_echantillons_ORIGINAL, aes(x = wavelength, y = R.Ave,
                                                     group = sample_id),
            linetype = "dotted", col = "coral", alpha = 0.5) +
  geom_line(data = cinq_echantillons_ORIGINAL, aes(x = wavelength, y = T.Ave,
                                                     group = sample_id),
            linetype = "dotted", col = "coral", alpha = 0.5) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  theme_bw()
rem_ovrlp_verif # satisfaisant, plus aucun chevauchement
#ggsave("rem_ovrlp_verif.pdf", plot = rem_ovrlp_verif, height = 10, width = 10)

# rassembler avec métadonnées ----
rt_join <- rt_v7_large %>%
  select(sample_id, scientific_name, site) %>% # , wavelength
  distinct()

rt_all_no_ovrlp <- rt_all_no_ovrlp_provisoire2 %>%  # ** sur rt_all...prov2 UNGROUPED ** sinon ça bogue
  full_join(rt_join, by = 'sample_id') %>%   # vérifier quelques valeurs inchangées
  dplyr::rename(wavelength = "wvl", reflectance = "refl", transmittance = "trans")
unique(rt_all_no_ovrlp$sample_id) # verif, 94 sample_id ok √
#save(rt_all_no_ovrlp, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_all_no_ovrlp.RData") # rt_all_no_ovrlp, fait le 29 mai 2019 / 21 juin !!!


#-=-=-=-=-=-=-=-=-=-=-=-=-
# création de rt_long_no_ovrlp à partir de rt_vx ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# format long tidy (reflectance et transmittance)
rt_long_no_ovrlp_provisoire <- rt_all_no_ovrlp %>%
  #plyr::rename(c("R.Ave" = "reflectance", "T.Ave" = "transmittance")) %>%
  #dplyr::select(-c(R.STD, T.STD)) %>%
  gather(key = propriete, value = valeur, reflectance, transmittance) %>%
  arrange(sample_id, scientific_name, site, propriete, valeur)

# retirer régions
rt_long_no_ovrlp <- rt_long_no_ovrlp_provisoire %>%
  dplyr::filter(wavelength %in% c(400:2400)) # filter les wvl
unique(sort(rt_long_no_ovrlp$wavelength)) # vérification
# [1]  400  401  402  403
max(sort(rt_long_no_ovrlp$wavelength)) # vérification
# [1] 2400
save(rt_long_no_ovrlp, file = "/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long_no_ovrlp.RData") # fait le 17 mai / refait le 29 mai avec rt_all_no_ovrlp pour enlever l'overlap *
