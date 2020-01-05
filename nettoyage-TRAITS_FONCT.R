# 8 MARS

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                                  NETTOYAGE DONNÉES
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : AVANT 11 JUIN, travaillé avec ft_all devenu ARCHIVE ***, nouveau s'apelle ft_all pour fin de simplicité

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Library
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(tidyverse) 
library(dplyr)
library(plyr)
library(googlesheets) # pour importer Google Sheets directement
library(magrittr)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données à charger et ménage ####
#-=-=-=-=-=-=-=-=-=-=-=-=-

# # charger raw_LAWS
# raw_LAWS <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./tele_leaf_area_and_water_samples.csv") # changer eriophorum à la main !
# # joindre sample et sample2
# paste_noNA <- function(x,sep=", ")
#   gsub(", " ,sep, toString(x[!is.na(x) & x!="" & x!="NA"] ) )
# raw_LAWS$sample_bon <- apply( raw_LAWS[ , c(19:20) ] , 1 , paste_noNA , sep="")
# raw_LAWS$sample_id <- as.character(raw_LAWS$sample_id)
# raw_LAWS <- raw_LAWS %>%
#   dplyr::select(c("sample_id", "specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2", "leaf_dry_matter_content_mg_g",
#          "leaf_water_content_mg_g", "equivalent_water_thickness_cm")) %>%
#   filter(! sample_id == 11998442) %>%  # retirer, pour l'instant, le S_ID 11998442 car récolté dans le buffer zone, et pas de données de LAWS dessus
#   dplyr::filter(!is.na(specific_leaf_area_m2_kg))
# 
# # vérifier si duplication de sample_id
# length(unique(raw_LAWS$sample_id))
# nrow(raw_LAWS)
# # pas same et pas pareil comme les autres
# # explcation : pcq on a des sphaignes la dedans et des temporal
# 
# 
# # charger raw_bulk (BUT : savoir le nom de l'espèce)
# raw_bulk <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./tele_bulk_leaf_samples.csv", sep = ";") # changer eriophorum à la main !
# raw_bulk$sample_id <- as.character(raw_bulk$sample_id)
# raw_bulk <- raw_bulk %>%
#   filter(! sample_id == 11998442) %>%  # retirer, pour l'instant, le S_ID 11998442 car récolté dans le buffer zone, et pas de données de LAWS dessus
#   dplyr::filter(! site_id == "MBP_temporal") %>%
#   dplyr::filter(! sample_id %in% c("9236658","9249843", "9320621", "9324461", "9335280",
#                                    "9338745", "9342418", "9668337", "9671728", "9674948")) %>%
#   dplyr::filter(! scientific_name == "Spinacia oleracea Linnaeus") %>% 
#   dplyr::select(c("plant_id", "scientific_name", "site_id", "sample_id", "date_sampled"))
#   
# 
# # charger CN
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/C:N/C_N_transfo.RData")
# C_N$sample_id <- as.character(C_N$sample_id) # convertir en numérique
# 
# # vérifier si duplication de sample_id
# length(unique(C_N$sample_id))
# nrow(C_N)
# # pas same et pas autant que les autres
# 
# 
# # charger carbon_fractions
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/carbon_fractions/c_fractions.RData")
# c_fractions$sample_id <- as.character(c_fractions$sample_id) # convertir en numerique
# 
# # vérifier si duplication de sample_id
# length(unique(c_fractions$sample_id))
# nrow(c_fractions)
# # same
# 
# 
# # charger raw_plants (BUT : savoir le plot et site (?) )
# raw_plants <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./raw_plants.csv", sep = ";")
# raw_plants <- raw_plants %>%
#   dplyr::select(c("plant_id", "plot_id", "plot_field_id"))
# # ajouter le nom du plot ? sinon juste inutile de mettre plot id vue qu'on a tt le reste de l'info anywy ?
# 
# 
# # charger  phénols tannins
# # ph_tann_ws <- gs_url("https://docs.google.com/spreadsheets/d/1NMtZMYLq6qkLOLtAR_gT9zpuPga_6Giu6tZ_zgkbOo8/edit#gid=0")
# # ph_tann_valeurs <- gs_read_csv(ph_tann_ws, ws = "valeurs")
# # ph_tann <- ph_tann_valeurs[1:122,] %>%
# #   select(-c("position_in_box", "box", "id_analyse/phénol", "run_no", "MARS sur quel ordi ?" , "statut_faire test_batch",
# #             "volume_extrait_ml","X33", "site_id", "actual_leaf_dry_matter_content_perc", "masse_g", "Tannins_A735",
# #             "Tannins_intercept_GAE", "Tannins_slope_GAE", "Tannins_r2_standard_curve", "Tannins_range_aborbance_std_curve",
# #             "Ph_conc_AFTER_precPVP_mgL", "Ph_conc_AFTER_precPVP_xfd_mgL", "Phenols_A735", "Phenols_intercept_GAE", "Phenols_slope_GAE",
# #             "Phenols_r2_standard_curve", "Ph_range_aborbance_std_curve", "Phenols_concentration_mgL", "Ph+Tann_concentration_xfd_mgL",
# #             "Tann_concentration_mgL", "volume_extrait_ml", "Tannins_pourc_verif", "Phenols_pourc_verif")) %>%
# #   dplyr::filter(!date_analyse %in% c("2019-03-05", NA))
# #save(ph_tann, file = "~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses/phenol_tannins_valeurs_pour_analyses.RData")
# load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/phenol_tannins_valeurs_pour_analyses.RData")
# 
# # vérifier si duplication de sample_id
# length(unique(ph_tann$sample_id))
# nrow(ph_tann) # ya des dupl
# 
# 
# # charger pigments
# # gs_auth()
# # pigments_ws <- gs_url("https://docs.google.com/spreadsheets/d/1CwTwJCtZ6HZmkJlPD1czNW3NOj8Ji308uQgHZ7FpidM/edit#gid=722152891")
# # pigments_valeurs <- gs_read_csv(pigments_ws, ws = "valeurs")
# # pigments <- pigments_valeurs %>%
# #   dplyr::filter(!chlB_mg_g == "-1.845") %>%
# #   dplyr::filter(!carotenoides_mg_g == "-1.530") %>%
# #   select(-c("position_in_box", "box", "id_analyse", "run_no", "MARS sur quel ordi ?" , "statut", "volume_extrait_ml","Notes",
# #             "site_id", "actual_dry_matter_cont_pourc", "masse_g", "id_analyse", "run_no",
# #             "MARS sur quel ordi ?", "actual_dry_matter_cont_pourc", "masse_g", "A470_mean_micropl", "pathlength_470_cm",
# #             "A652_mean_micropl", "pathlength_652_cm", "A665_mean_micropl", "pathlength_665_cm", "A710_mean_micropl",
# #             "pathlength_710_cm", "A470_cm", "A652_cm", "A665_cm", "A710_cm", "chlA_brut_μg_ml", "chlB_brut_μg_ml",
# #             "carotenoides_brut_μg_ml", "fact_dilution", "volume_extrait_ml")) %>%
# #            dplyr::filter(!sample_id %in% c("ChlA", "Maison (épinard)", "S1-ChlA", "S2-Maison (épinard0", "S2-Maison(épinard)", NA))
# # save(pigments, file = "~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses/pigments_valeurs_pour_analyses.RData")
# load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/pigments_valeurs_pour_analyses.RData")
# 
# # vérifier si duplication de sample_id
# length(unique(pigments$sample_id))
# nrow(pigments)
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# # Joindre tous les traits fonctionnels et colonnes d'intérêt ####
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# 
# # joindre les tables
# ft_all_provisoire <- full_join(raw_plants, raw_bulk, by = "plant_id") %>% # obtenir plot_id et plot field id
#   full_join(raw_LAWS, by = "sample_id") %>%
#   full_join(C_N, by = "sample_id") %>%
#   full_join(c_fractions, by = "sample_id") %>%
#   full_join(pigments, by = "sample_id") %>%
#   full_join(ph_tann, by = "sample_id") %>%
# #  full_join(phosph_oligo, by = "sample_id") %>%
#   filter(! sample_id == 11998442) %>%  # retirer, pour l'instant, le S_ID 11998442 car récolté dans le buffer zone, et pas de données de LAWS dessus
#   dplyr::filter(! scientific_name.y == "Sphagnum capillifolium (Ehrh.) Hedw.") %>%
#   dplyr::filter(! site_id == "MBP_temporal") %>%
#   dplyr::select(- c("scientific_name.x", "date_analyse.x","date_analyse.y", "scientific_name.y"))
# 
# 
# # mettre tout en numérique
# cols = c("specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2","leaf_dry_matter_content_mg_g",
#                     "leaf_water_content_mg_g","equivalent_water_thickness_cm", "N_pourc","C_pourc",
#                     "soluble_contents","hemicellulose","cellulose","lignine","recalcitrant","chlA_mg_g",
#                     "chlB_mg_g","carotenoides_mg_g", "Tannins_mg_g","Phenols_mg_g")
# ft_all_provisoire[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# # Assurer que les sample id sont uniques et correction problèmes ####
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# 
# # vérifier si duplication de sample_id
# length(unique(ft_all_provisoire$sample_id))
# nrow(ft_all_provisoire) # ya des duplica !!!
# duplicated(ft_all_provisoire$sample_id)
# 
# 
# # problème 1. (ligne 10)
# # 11667486 mesuré 2 fois le spectre, pas repris de bulk
# # solution temporaire : retirer ce réplicat
# # MBP_open, plot 1
# ft_all_provisoire_2 <- ft_all_provisoire %>%
#   dplyr::filter(! sample_id == 11667486)
# 
# 
# # problème 2. (lignes 68-72)
# # 13320387 j'ai 2 cryotubes de leaf disks donc pigments et ph&tann ont été faits 2 fois.
# # solution aucune raison de croire que l'un ou l'autre est meilleur, explication pour avoir 2 tubes, premier tube remplis, j'en ai remplis 2.
# # je fais une moyenne des traits (pigments, ph/tann) pour cet échantillon
# 
# # vérifier si duplication de sample_id
# length(unique(ft_all_provisoire_2$sample_id))
# nrow(ft_all_provisoire_2) # ya des duplica !!!, lesquels ?
# duplicated(ft_all_provisoire_2$sample_id) # il n'y a bien que 13320387 qui est répliqué
# 
# ft_all_provisoire_3 <- ft_all_provisoire_2 %>%
#   group_by(sample_id) %>%
#   mutate_at(c("chlA_mg_g","chlB_mg_g","carotenoides_mg_g","Tannins_mg_g","Phenols_mg_g"), mean) %>%
#   filter(!duplicated(sample_id)) # enlever les duplica ddu sample id restant
# 
# duplicated(ft_all_provisoire_3$sample_id)
# # aucun duplicat
# 
# # vérifier si la moyenne a fonctionné pour échantillon 13320387
# ft_all_13320387 <- ft_all_provisoire_3 %>%
#   filter(sample_id == 13320387)
# mean(ft_all_13320387$chlA_mg_g)
# # [1] 2.3375
# ft_all_13320387$chlA_mg_g
# # [1] 2.3375
# # c'est beau
# 
# 
# # problème 3. 
# # échantillon # 11677595 et 12193094 mm plante, 2 mesures valides de spectres
# # faire moyenne des traits
# ft_11677595_12193094 <- ft_all_provisoire_3 %>% 
#   dplyr::filter(sample_id %in% c("11677595","12193094"))
# # traits soit identiques (doublés) soit NA soit similaires (disques)
# 
# ft_11677595_12193094_mean <- ft_11677595_12193094 %>% 
#   group_by(plant_id) %>% 
#   mutate_at(c("specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2" ,"leaf_dry_matter_content_mg_g", 
#               "leaf_water_content_mg_g", "equivalent_water_thickness_cm", "N_pourc", "C_pourc", 
#               "soluble_contents" ,"hemicellulose", "cellulose", "lignine", "recalcitrant", "chlA_mg_g",
#               "chlB_mg_g", "carotenoides_mg_g", "Tannins_mg_g", "Phenols_mg_g"), mean, na.rm = T) %>% 
#   dplyr::filter(sample_id == 11677595)
# 
# # rejoindre ft_all_provisoire 3 à ft_11677595_12193094_mean
# ft_all_provisoire_4_provisoire <- ft_all_provisoire_3 %>% 
#   dplyr::filter(!plant_id == 11677538)
# 
# ft_all_provisoire_4 <- full_join(ft_all_provisoire_4_provisoire, ft_11677595_12193094_mean)
# unique(ft_all_provisoire_4$plant_id)
# duplicated(ft_all_provisoire_4$plant_id)
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# # Version finale ####
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# #
# # renommer, version finale
# ft_all <- ft_all_provisoire_4 %>%
#   dplyr::select("site_id", "plant_id","sample_id","date_sampled","scientific_name", "plot_id", "plot_field_id", "specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2","leaf_dry_matter_content_mg_g",
#            "leaf_water_content_mg_g","leaf_water_content_mg_g", "equivalent_water_thickness_cm", "N_pourc","C_pourc",
#            "soluble_contents","hemicellulose","cellulose","lignine","recalcitrant","chlA_mg_g",
#            "chlB_mg_g","carotenoides_mg_g","Phenols_mg_g", "Tannins_mg_g") # mettre en ordre cute pour aucune autre raison
# ft_all$site_id = droplevels(ft_all)$site_id # retirer les "levels" de site_id non-utilisés
# 
#
# #-=-=-=-=-=-=-=-=-=-=-=-=-
# # Obtenir nouvelles variables ####
# #-=-=-=-=-=-=-=-=-=-=-=-=-
#
# ouvrir dernière version de ft_all
#load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses.RData")
#
# # obtenir C:N ####
# ft_all <- ft_all %>% 
#   dplyr::mutate(C_N_ratio = C_pourc/N_pourc)
# 
#
# # obtenir ChlA:ChlB ####
# ft_all <- ft_all %>% 
#   dplyr::mutate(chlA_chlB = chlA_mg_g/chlB_mg_g)
#
#
# obtenir [pigment chla, chlb] content per unit area (µmol cm-2) ####

# ouvrir dernière version de ft_all
load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/ARCHIVEtous_analyses.RData") # ARCHIVE_ft_all

# constantes et extraction des valeurs à convertir
A_content_mg_g <- ARCHIVE_ft_all %>% 
  dplyr::select(chlA_mg_g, sample_id)
B_content_mg_g <- ARCHIVE_ft_all %>% 
  dplyr::select(chlB_mg_g, sample_id)
SLA <- ARCHIVE_ft_all %>% 
  dplyr::select(specific_leaf_area_m2_kg, sample_id)
SLA$SLA_m2_g <- SLA$specific_leaf_area_m2_kg/1000

whole <- full_join(A_content_mg_g, B_content_mg_g) %>% 
  full_join(SLA)

# calcul pigments (mg) / aire (m2)
whole$chlA_mg_m2 <- whole$chlA_mg_g / whole$SLA_m2_g 
whole$chlB_mg_m2<- whole$chlB_mg_g / whole$SLA_m2_g 
whole <- whole %>% dplyr::select(chlA_mg_m2, chlB_mg_m2, sample_id)

# joindre nouvelles colonnes A et B à ft_all (créer ft_all_vX1)
ft_all_vX1 <- ARCHIVE_ft_all %>% 
  full_join(whole) 

ft_all_vX1[, "recalcitrant"] <- if_else(as.numeric(ft_all_vX1[, "recalcitrant"]) < 0, 0, (as.numeric(ft_all_vX1[, "recalcitrant"])))
ft_all_vX2 <- ft_all_vX1

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Enregistrement final (save) ####
#-=-=-=-=-=-=-=-=-=-=-=-=-

# NE PAS TOUCHER ARCHIVE_ft_all, ARCHIVE // SI JE VEUX LE MODIFIER, rechanger son nom pour ft_all, rerouler TOUT le script
#save(ARCHIVE_ft_all, file = "~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/ARCHIVEtous_analyses.RData") # refait 7 mai, ajouté tf "WATER"
# refait 9 mai, ajouté droplevels pour enlever "labo-xxx" 
# refait 28 mai, ajouté ratio ChlA/ChlB
# refait au complet 10 juin, arrangé problèmes de plot_id à GPB / et ajouté "plot_id" pour analyses LMM avec environnement

# NOUVEAU ft_all
# save(ft_all_vX2, file = "~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData") # fait 6 sept (ajouté) pigment/aire / CALICE refait 22oct
# load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData")


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Autres vérifications ####
#-=-=-=-=-=-=-=-=-=-=-=-=-

# # retouver mes 6 Rhodo à MBP √ c beau 1 retiré car tt laid (plot1)
# verif <- raw_plants %>% 
#   full_join(raw_bulk) %>% 
#   dplyr::select(plant_id, plot_id, plot_field_id, plant, plant2, scientific_name, site_id, sample_id) %>% 
#   dplyr::filter(!scientific_name %in% c("Sphagnum capillifolium (Ehrh.) Hedw.", "Spinacia oleracea Linnaeus")) %>% 
#   dplyr::filter(!sample_id == "NA") %>% 
#   dplyr::filter(scientific_name == "Rhododendron groenlandicum (Oeder) Kron & Judd") %>% 
#   dplyr::filter(site_id == "MBP_open")
# 
# # réorganiser plots 1-2 GPB
# GPB <- ft_all %>% 
#   filter(site_id == "GPB_canal") %>% 
#   filter(plot_field_id == "Plot1") %>% 
#   filter(! scientific_name == "Rhododendron groenlandicum")
# # si j'ai 6 items, voir sur fulcrum les BONS no de plots, pas dans les fichiers de téléchargment !
# # les changer manuellement...

