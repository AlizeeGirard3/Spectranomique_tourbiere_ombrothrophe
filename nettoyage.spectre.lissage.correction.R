# 16 mai 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                   Lissage données spectrales
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : voir https://github.com/elaliberte/warren-spectra 

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
#if (!'signal' %in% installed.packages() ) install.packages("signal")

library(data.table)
library(signal) # sgolay()
                # masked from ‘package:stats’: filter, poly
library(plyr)
library(tidyverse) 
library(dplyr)
library(magrittr)
library(spectrolab) # as.spectra(), plot(),  / Meireles J, Schweiger A, Cavender-Bares J (2018). spectrolab: Class
# and Methods for Hyperspectral Data. R package version 0.0.8, <URL:
#  https://CRAN.R-project.org/package=spectrolab>.
# masked from ‘package:signal’: resample **
library(mgsub) # multiple pattern, multiple replacements mgsub()
library(reshape2) # dcast


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_lisse.RData") # rt_sg_lisse
# obtenu via manips sur ce script sur # rt_all_no_ovrlp et/ou rt_long_no_overlp

load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_large.RData") # rt_sg_large  
# obtenu via manips sur ce script sur # rt_sg_lisse (son équivalent, mais long)


load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_all_no_ovrlp.RData") # rt_all_no_overlp (données originales sans overlap)
load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long_no_ovrlp.RData") # rt_long_no_overlp (données originales sans overlap)
# fichiers issus du script : " nettoyage-spectre.R" dans document "DonnéesAnalyse/scripts-NETTOYAGE"

# PLUS VIEUX
load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long_lisse.RData") # rt_long_lisse (données lissée avec smooth())
load("~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_lisse.RData") # rt_sg_lisse (données lissée avec smooth() + s-golay)
# issus des manips sur rt_long et/ou rt_v6 dans ce script *


#=======================================# SMOOTH  #=======================================#                    ----

# les régions de transition entre détecteurs sont laids, il faut les smoother de façon  aggresive
# EST-CE QUE JE L'APPLIQUE AVEC SGOLAY OU SEULEMENT SGOLAY ? SGOLAY seulement n'est pas satisfaisant !

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Format large ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# large, réflectance
spectra_large_refl_provisoire <- dcast(rt_all_no_ovrlp, site + scientific_name + sample_id ~ wavelength, 
                                       fun.aggregate = mean,
                                       value.var = "reflectance")
spectra_large_refl_provisoire2 <- spectra_large_refl_provisoire %>% 
  dplyr::select(c("site", "scientific_name", "sample_id", "400":"2400")) # filter les wvl


# large, transmittance
spectra_large_trans_provisoire <- dcast(rt_all_no_ovrlp, site + scientific_name + sample_id ~ wavelength, 
                                        fun.aggregate = mean,
                                        value.var = "transmittance")
spectra_large_trans_provisoire2 <- spectra_large_trans_provisoire %>% 
  dplyr::select(c("site", "scientific_name", "sample_id", "400":"2400")) # filter les wvl

# transformer les spectra_large en objet "spectra"
spectra_large_refl <- as.spectra(spectra_large_refl_provisoire2, name_idx = 3, meta_idxs = c(1,2)) # utiliser no de colone *
spectra_large_trans <- as.spectra(spectra_large_trans_provisoire2, name_idx = 3, meta_idxs = c(1,2)) 

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Lissage ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

spectra_large_refl_smooth = smooth(spectra_large_refl[ , 400:2400 ]) # smooth {spectrolab}, using = 'spline' (défaut) : 
#plot_interactive(spectra_large_refl_smooth) # j'ai gossé pi ça change rien 

spectra_large_trans_smooth = smooth(spectra_large_trans[ , 400:2400 ]) 
#plot_interactive(spectra_large_trans_smooth) # j'ai gossé pi ça change rien 


#=======================================# LISSAGE + SAVITZKY-GOLAY #=======================================#          ----

# ce filtre pour le reste du spectre, conserve mieux la forme
# * These filters are particularly good at preserving lineshape while removing high frequency squiggles.

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Format "matrice" (à partir de l'objet "spectra" lissé) ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# reflectance
refl_liss <- as.matrix(spectra_large_refl_smooth) 
sample_id_provisoire <- rownames(refl_liss) # spec_xxxxxxx en "rowname", retransformer en vecteur
sample_id <- gsub( "spec_", "", as.character(sample_id_provisoire)) # enlever "spec_"
refl_liss.1 <- data.table(refl_liss)
refl_liss2 <- cbind(sample_id, refl_liss.1) # ajouter au df

# transmittance
trans_liss <- as.matrix(spectra_large_trans_smooth) 
sample_id_provisoire <- rownames(trans_liss) # spec_xxxxxxx en "rowname", retransformer en vecteur
sample_id <- gsub( "spec_", "", as.character(sample_id_provisoire)) # enlever "spec_"
trans_liss.1 <- data.table(trans_liss)
trans_liss2 <- cbind(sample_id, trans_liss.1) # ajouter au df


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Un fichier RT (à partir de l'objet lissé) ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# colonne réfl dans le df réfl
propriete <- rep("reflectance", times = nrow(refl_liss2))  
refl_liss3 <-  refl_liss2 %>% 
  cbind(propriete) %>% 
  dplyr::select("sample_id", "propriete", c("400":"2400"))

# colonne trans dans le df trans
propriete <- rep("transmittance", times = nrow(trans_liss2))  
trans_liss3 <-  trans_liss2 %>% 
  cbind(propriete) %>% 
  dplyr::select("sample_id", "propriete", c("400":"2400"))

# joindre reflectance et transmittance
rt_lisse_provisoire1 <- refl_liss3 %>% 
  full_join(trans_liss3)

# obtenir site 
raw_bulk <- read.csv("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/DONNÉES_BRUTES_fulcrum_ggsheets/interface_telech./tele_bulk_leaf_samples.csv", sep = ";")
site <- raw_bulk %>%
  dplyr::filter(scientific_name != "Sphagnum capillifolium (Ehrh.) Hedw.") %>%
  dplyr::filter(! site_id %in% c("MBP_temporal", "Labo-B234")) %>%
  dplyr::filter(! status == "deleted") %>%
  dplyr::filter(! sample_id %in% c("9236658", "9249843", "9320621", "9324461", "9335280", "9338745", "9342418", "9668337", "9671728", "9674948")) %>%
  dplyr::select(sample_id, site_id) %>%  # ne pas mettre scientfic names ! si j'en ai besoin, faire un autre objet
  dplyr::rename(site = "site_id")

site$site = droplevels(site$site) # éliminer levels inutiles
levels(site$site) # vérif
site$sample_id <- as.character(site$sample_id)

# obtenir scientific_name
scientific_name <- raw_bulk %>%
  dplyr::filter(! scientific_name %in% c("Sphagnum capillifolium (Ehrh.) Hedw.", "Spinacia oleracea Linnaeus")) %>%
  #dplyr::filter(! site_id %in% c("MBP_temporal", "Labo-B234")) %>%
  dplyr::filter(! status == "deleted") %>%
  dplyr::filter(! sample_id %in% c("9236658", "9249843", "9320621", "9324461", "9335280", "9338745", "9342418", "9668337", "9671728", "9674948")) %>%
  dplyr::select(sample_id, scientific_name) 

# simplifier noms d'espèces
unique(scientific_name$scientific_name)
scientific_name$scientific_name <- plyr::revalue(scientific_name$scientific_name, # revalue() => pkg "plyr"
                                                 c("Eriophorum vaginatum Linnaeus subsp. Vaginatum" = "Eriophorum vaginatum",
                                                   "Eriophorum vaginatum Linnaeus subsp. vaginatum" = "Eriophorum vaginatum",
                                                   #'Kalmia angustifolia Linnaeus' = 'Kalmia angustifolia',
                                                   'Kalmia angustifolia Linnaeus var. angustifolia' = 'Kalmia angustifolia',
                                                   'Rhododendron groenlandicum (Oeder) Kron & Judd' = 'Rhododendron groenlandicum',
                                                   'Chamaedaphne calyculata (Linnaeus) Moench' = 'Chamaedaphne calyculata'))


scientific_name$scientific_name = droplevels(scientific_name$scientific_name)
levels(scientific_name$scientific_name) # vérif
scientific_name$sample_id <- as.character(scientific_name$sample_id)
cols = c("sample_id", "scientific_name")
scientific_name[,cols] %<>% lapply(function(x) as.character(x))

# joindre site et scientific_name à rt_lisse
rt_lisse_provisoire2 <- left_join(rt_lisse_provisoire1, site, by = 'sample_id') # left_join enlèvera les sample_id qui sont dans gauche et pas le droite
rt_lisse <- left_join(rt_lisse_provisoire2, scientific_name, by = 'sample_id') %>% 
  dplyr::select(sample_id, scientific_name, site, propriete, "400":"2400")


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Format long, lissé ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

rt_long_lisse <- gather(rt_lisse, key = "wavelength", value = "valeur", c("400":"2400"), -sample_id, -scientific_name, -site)
rt_long_lisse$wavelength <- as.numeric(rt_long_lisse$wavelength)
# save(rt_long_lisse, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_long_lisse.RData") # rt_long_lisse / refait le 29 mai avec matched_sensors() / refait 21 juin données officielles

#-=-=-=-=-=-=-=-=-=-=-=-=-
# VÉRIFICATION DES PARAMÈTRES  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# NOTE : * These filters are particularly good at preserving lineshape while removing high frequency squiggles.
# DONNÉES LISSÉES ("spectrolab::smooth()") SEULEMENT

# déterminer région et "n" (largeur de fenêtre)
cinq_echantillons_r <- rt_long_lisse %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "13224081", "14009324", "11929171"))
cinq_echantillons_t <- rt_long_lisse %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))

cinq_echantillons_r_ORIGINAL <- rt_long_no_ovrlp %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "13224081", "14009324", "11929171"))
cinq_echantillons_t_ORIGINAL <- rt_long_no_ovrlp %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))


cinq_echantillons_rt <- ggplot() +
  geom_line(data = cinq_echantillons_r, aes(x = wavelength, y = valeur, 
                                            group = sample_id)) +
  geom_line(data = cinq_echantillons_t, aes(x = wavelength, y = valeur, 
                                            group = sample_id)) +
  geom_line(data = cinq_echantillons_r_ORIGINAL, aes(x = wavelength, y = valeur, 
                                                     group = sample_id),
            linetype = "dotted", col = "coral", alpha = 0.5) +
  geom_line(data = cinq_echantillons_t_ORIGINAL, aes(x = wavelength, y = valeur, 
                                                     group = sample_id), 
            linetype = "dotted", col = "coral", alpha = 0.5) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  theme_bw() 
cinq_echantillons_rt
#ggsave("rt_lisse.pdf", plot = cinq_echantillons_rt, height = 10, width = 10)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Savitzky-golay (sur données lissées) ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# fonction modifiée d'Étienne                  # permet de travailler sur sous-ensemble du df et d'obtenir directement nouvelles colonnes
filtre_sg <- function(x, p, n) {              
  x$valeur_sg <- sgolayfilt(x$valeur, p, n)    # sur matrice "x", fait une nouvelle colone "sg_valeur" avec la fonction sgolayfilt et les paramètres à préciser "p" et "n"
  return(x)                                    # renvoie la matrice "x" contenant la nouvelle colonnes contenant les données lissées
}                                              # odd number : chiffre impair

# créer sous-ensemble de longueur d'ondes et appliquer fonction
sg_data_VIS_NIR1_lis <- rt_long_lisse %>%
  dplyr::filter(wavelength < 900) %>% 
  group_by(sample_id, site, scientific_name, propriete) %>% 
  do(filtre_sg(., p = 3, n = 21))                                          # n = largeur de la fenêtre (Sadeghi, M., & Behnia, F. (2018). Optimum window length of Savitzky-Golay filters with arbitrary order. arXiv preprint arXiv:1808.10489.)

# SOUS ENSEMBLE -> détecteur NIR
sg_data_NIR_detect_lis <- rt_long_lisse %>%
  dplyr::filter(wavelength >= 900,
                wavelength <= 1090 ) %>%
  group_by(sample_id, site, scientific_name, propriete) %>% 
  do(filtre_sg(., p = 3, n = 55))   
# vérif
#unique(sg_data_NIR_detect_lis$wavelength)

# SOUS-ENSEMBLE -> après le détecteur NIR                     
sg_data_NIR3_lis <- rt_long_lisse %>%
  dplyr::filter(wavelength > 1090,                         
                wavelength <= 1390) %>%
  group_by(sample_id, site, scientific_name, propriete) %>% 
  do(filtre_sg(., p = 3, n = 21))                                                              
# vérif
#unique(sg_data_NIR3_lis$wavelength)

sg_data_SWIR1_lis <- rt_long_lisse %>%
  dplyr::filter(wavelength > 1390,
                wavelength < 1880) %>%   # SI PAS ACCEPTABLE, RETOURNER À 1880 et refaire figure !
  group_by(sample_id, site, scientific_name, propriete) %>% 
  do(filtre_sg(., p = 3, n = 21))        # AVEC p = 4, n = 175 (plus petite fenêtre), lisse inutilement ! 
# vérif
#unique(sg_data_SWIR1_lis$wavelength)

sg_data_SWIR2_lis <- rt_long_lisse %>%
  dplyr::filter(wavelength > 1880) %>%       
  group_by(sample_id, site, scientific_name, propriete) %>% 
  do(filtre_sg(., p = 4, n = 175))             # n = 175 min avec p = 5... marche pas mm avec p = 3...
#unique(sg_data_SWIR2_lis$wavelength)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# obtenir rt_sg_lisse ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

rt_sg_lisse <- bind_rows(sg_data_VIS_NIR1_lis, sg_data_NIR_detect_lis,
                         sg_data_NIR3_lis,sg_data_SWIR1_lis,sg_data_SWIR2_lis) %>%
  dplyr::select(sample_id, scientific_name, site, propriete, wavelength, valeur, valeur_sg) 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# VISU effet LISSAGE + SGOLAY ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# VIS-NIR1 et NIR3 ----

VIS_NIR1_5_ech_r_lis <- sg_data_VIS_NIR1_lis %>%                 #√ CADUQUE SI PARAMÈTRES VIS=NIR1
  dplyr::filter(propriete == 'reflectance',
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))
NIR3_5_ech_r_lis <- sg_data_NIR3_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))

VIS_NIR1_5_ech_t_lis <- sg_data_VIS_NIR1_lis %>%                   #√ CADUQUE SI PARAMÈTRES VIS=NIR1
  dplyr::filter(propriete == 'transmittance',
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))
NIR3_5_ech_t_lis <- sg_data_NIR3_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))

VN13_21_rt <- ggplot() +    
  geom_line(data = NIR3_5_ech_r_lis, aes(x = wavelength,                     # REFL
                                         y = valeur_sg, 
                                         group = sample_id), col = "black") +     
  geom_line(data = NIR3_5_ech_r_lis, aes(x = wavelength, 
                                         y = valeur,                                            # REFL
                                         group = sample_id), linetype = 'dotted', col = "blue") +    
  geom_line(data = NIR3_5_ech_t_lis, aes(x = wavelength,
                                         y = valeur_sg,                                           # TRANS
                                         group = sample_id), col = "black") + 
  geom_line(data = NIR3_5_ech_t_lis, aes(x = wavelength,                      
                                         y = valeur,                                             # TRANS
                                         group = sample_id), linetype = 'dotted', col = "blue") +   
    geom_line(data = VIS_NIR1_5_ech_r_lis, aes(x = wavelength,                     # REFL
                                           y = valeur_sg,
                                           group = sample_id), col = "black") +
    geom_line(data = VIS_NIR1_5_ech_r_lis, aes(x = wavelength,
                                           y = valeur,                                            # REFL
                                           group = sample_id), linetype = 'dotted', col = "blue") +
    geom_line(data = VIS_NIR1_5_ech_t_lis, aes(x = wavelength,
                                           y = valeur_sg,                                           # TRANS
                                           group = sample_id), col = "black") +
    geom_line(data = VIS_NIR1_5_ech_t_lis, aes(x = wavelength,
                                           y = valeur,                                             # TRANS
                                           group = sample_id), linetype = 'dotted', col = "blue") +
   xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  #coord_cartesian(ylim = c(0.3, 0.75), xlim = c(900, 1000)) +  
  theme_bw() 
VN13_21_rt  # très fidèle à la courbe


# NIR DÉTECTEUR ----

NIRd_5_ech_r_lis <- sg_data_NIR_detect_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))
NIRd_5_ech_t_lis <- sg_data_NIR_detect_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))

Nd_55_rt <- ggplot() +    
  geom_line(data = NIRd_5_ech_r_lis, aes(x = wavelength,                     # REFL
                                         y = valeur_sg, 
                                         group = sample_id), col = "black") +     
  geom_line(data = NIRd_5_ech_r_lis, aes(x = wavelength, 
                                         y = valeur,                                            # REFL
                                         group = sample_id), linetype = 'dotted', col = "blue") +    
  geom_line(data = NIRd_5_ech_t_lis, aes(x = wavelength,
                                         y = valeur_sg,                                           # TRANS
                                         group = sample_id), col = "black") + 
  geom_line(data = NIRd_5_ech_t_lis, aes(x = wavelength,                      
                                         y = valeur,                                             # TRANS
                                         group = sample_id), linetype = 'dotted', col = "blue") +   
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  #coord_cartesian(ylim = c(0.3, 0.75), xlim = c(900, 1000)) +  
  theme_bw() 
Nd_55_rt         # 55 semble stisfaisant....
#ggsave("NIR 1 section.pdf", plot = N55_rt, height = 10, width = 10)

# VIS-NIR1, NIR détecteur et NIR3 (ALL) ----

VIS_NIR1_5_ech_r_lis <- sg_data_VIS_NIR1_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))
NIR3_5_ech_r_lis <- sg_data_NIR3_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))
NIRd_5_ech_r_lis <- sg_data_NIR_detect_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12173174", "14009324", "11929171"))
NIR1d3_5_ech_r_lis_provis <- full_join(VIS_NIR1_5_ech_r_lis, NIR3_5_ech_r_lis)
NIR1d3_5_ech_r_lis <- full_join(NIR1d3_5_ech_r_lis_provis, NIRd_5_ech_r_lis)

VIS_NIR1_5_ech_t_lis <- sg_data_VIS_NIR1_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))
NIR3_5_ech_t_lis <- sg_data_NIR3_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))
NIRd_5_ech_t_lis <- sg_data_NIR_detect_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("12007476", "12086180","13224081", "13301329", "13318392"))
NIR1d3_5_ech_t_lis_provis <- full_join(VIS_NIR1_5_ech_t_lis, NIR3_5_ech_t_lis)
NIR1d3_5_ech_t_lis <- full_join(NIR1d3_5_ech_t_lis_provis, NIRd_5_ech_t_lis)


Nall <- ggplot() +    
  geom_line(data = NIR1d3_5_ech_r_lis, aes(x = wavelength,                     # REFL
                                           y = valeur_sg, 
                                           group = sample_id), col = "black") +     
  geom_line(data = NIR1d3_5_ech_r_lis, aes(x = wavelength, 
                                           y = valeur,                                            # REFL
                                           group = sample_id), linetype = 'dotted', col = "blue") +    
  geom_line(data = NIR1d3_5_ech_t_lis, aes(x = wavelength,
                                           y = valeur_sg,                                           # TRANS
                                           group = sample_id), col = "black") + 
  geom_line(data = NIR1d3_5_ech_t_lis, aes(x = wavelength,                      
                                           y = valeur,                                             # TRANS
                                           group = sample_id), linetype = 'dotted', col = "blue") +   
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  #coord_cartesian(ylim = c(0.3, 0.75), xlim = c(900, 1000)) +  
  theme_bw() 
Nall  
#ggsave("VIS-NIR1_NIRd_NIR3-3sections.pdf", plot = Nall, height = 10, width = 10)

# NIR COMPLET CONCLUSION : 3 SECTIONS, ben beau                                 <<= <<= <<= <<= <<= <<= <<= <<=    
# VOIR PARAMÈTRES SECTION "# Savitzky-golay (sur données lissées)"      <<= <<= <<= <<= <<= <<= <<= <<=


# SWIR 1 ----

SWIR1_10_ech_r_lis <- sg_data_SWIR1_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR1_10_ech_t_lis <- sg_data_SWIR1_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))

S1_21 <- ggplot() +    # rouler fonction et rouler "sg_data_NIR" seulement, puis valeur_sg ici
  geom_line(data = SWIR1_10_ech_r_lis, 
            aes(group = sample_id, x = wavelength, y = valeur_sg), 
            col = "black") + 
  geom_line(data = SWIR1_10_ech_t_lis, 
            aes(group = sample_id, x = wavelength, y = valeur_sg), 
            col = "black") +           # valeur originale réflectance lissées
  geom_line(data = SWIR1_10_ech_r_lis, 
            aes(group = sample_id, x = wavelength, y = valeur), 
            linetype = "dotted", col = "coral") + 
  geom_line(data = SWIR1_10_ech_t_lis, 
            aes(group = sample_id, x = wavelength, y = valeur), 
            linetype = "dotted", col = "coral") + 
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  #coord_cartesian(ylim = c(0, 1), xlim = c(1390, 1880)) +  
  theme_bw() 
S1_21           # sert pas mal à rien avec smooth, mm ne suit pas complètement la ligne...
#ggsave("SWIR1.pdf", plot = S1_21, height = 10, width = 10)


# SWIR1 CONCLUSION : n = 21 <<= <<= <<= <<= <<= <<= <<= <<= <<= <<= <<=    

# SWIR 2 ----

SWIR2_10_ech_r_lis <- sg_data_SWIR2_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR2_10_ech_t_lis <- sg_data_SWIR2_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))

S2_175 <- ggplot(SWIR2_10_ech_r_lis,
                 aes(group = sample_id)) +    # rouler fonction et rouler "sg_data_NIR" seulement, puis valeur_sg ici
  geom_line(aes(x = wavelength, 
                y = valeur_sg), col = "black") + 
  geom_line(aes(x = wavelength, 
                y = valeur), linetype = 'dotted', col = "blue") +                     # valeur originale réflectance
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  # coord_cartesian(ylim = c(0, 1), xlim = c(1920, 1930)) +  
  theme_bw() 
S2_175
# parfait pour refl

S2_175_t <- ggplot(SWIR2_10_ech_t_lis,
                   aes(group = sample_id)) +    # rouler fonction et rouler "sg_data_NIR" seulement, puis valeur_sg ici
  geom_line(aes(x = wavelength, 
                y = valeur_sg), col = "black") + 
  geom_line(aes(x = wavelength, 
                y = valeur), linetype = 'dotted', col = "blue") +                     # valeur originale réflectance
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  #coord_cartesian(ylim = c(0, 1), xlim = c(1390, 1880)) +  
  theme_bw() 
S2_175_t
# sweet aussi pour trans

# SWIR2 CONCLUSION : n = 175 <<= <<= <<= <<= <<= <<= <<= <<= <<= <<= <<=    

# SWIR1 et SWIR2 ----

SWIR1_10_ech_r_lis <- sg_data_SWIR1_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR2_10_ech_r_lis <- sg_data_SWIR2_lis %>% 
  dplyr::filter(propriete == 'reflectance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR21_10_ech_r_lis <- full_join(SWIR1_10_ech_r_lis, SWIR2_10_ech_r_lis)

SWIR1_10_ech_t_lis <- sg_data_SWIR1_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR2_10_ech_t_lis <- sg_data_SWIR2_lis %>% 
  dplyr::filter(propriete == 'transmittance', 
                sample_id %in% c("13383299", "11995935", "12193094", "14009324", "11929171", 
                                 "12007476", "12086180","13224081", "13301329", "13318392"))
SWIR21_10_ech_t_lis <- full_join(SWIR1_10_ech_t_lis, SWIR2_10_ech_t_lis)

Sall <- ggplot() +    
  geom_line(data = SWIR21_10_ech_r_lis, aes(x = wavelength,                     # REFL
                                            y = valeur_sg, 
                                            group = sample_id), linetype = 'dotdash', col = "black") +     
  geom_line(data = SWIR21_10_ech_r_lis, aes(x = wavelength, 
                                            y = valeur,                                            # REFL
                                            group = sample_id), linetype = 'dotted', col = "coral") +    
  geom_line(data = SWIR21_10_ech_t_lis, aes(x = wavelength,
                                            y = valeur_sg,                                           # TRANS
                                            group = sample_id), linetype = 'dotdash', col = "black") + 
  geom_line(data = SWIR21_10_ech_t_lis, aes(x = wavelength,                      
                                            y = valeur,                                             # TRANS
                                            group = sample_id), linetype = 'dotted', col = "coral") +   
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  geom_vline(aes(xintercept = 1880), linetype = 'dashed', size = 0.5) +
  theme_bw() 
Sall        # p = 4, n = 175 jonction pourrait être mieux, mais pas si pire !, j'ai exclu 1880.... ça fait + cute... et ça parraît pas     
#ggsave("SWIR12.PDF", plot = Sall, height = 10, width = 10)


# SWIR1 et 2 CONCLUSION : p = 3, n = 21, p = 4, n = 175                 <<= <<= <<= <<= <<= <<= <<= <<=    
# VOIR PARAMÈTRES SECTION "# Savitzky-golay (sur données lissées)"      <<= <<= <<= <<= <<= <<= <<= <<=


#=======================================# VALIDATION #=======================================#                                ----

# ROULER AVEC GROS ORDI
# ici réalisé avec données lissées ("spectrolab::smooth()")

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Visualisation tous spectres
#-=-=-=-=-=-=-=-=-=-=-=-=-
reflectance <- rt_sg_lisse %>%
  dplyr::filter(propriete == 'reflectance') 

transmittance <- rt_sg_lisse %>%
  dplyr::filter(propriete == 'transmittance')

x <- ggplot() +
  geom_line(aes(x = wavelength, y = valeur_sg), data = reflectance) +                           
  geom_line(aes(x = wavelength, y = valeur_sg), data = transmittance) +
  xlab('Wavelength (nm)') +
  ylab('Reflectance (or 1 - Transmittance)') +
  labs(title = "Tous spectres avec données lissées, 23 mai 2019") +
  facet_wrap(~ sample_id)
x


#=======================================# GRAPH PARAMÈTRES FINAUX #=======================================#                          ----

# ici réalisé avec données lissées ("spectrolab::smooth()") et sans chevauchement ("spectrolab::match_sensors()")

#-=-=-=-=-=-=-=-=-=-=-=-=-
# GRAPH AVEC PARAMÈTRES FINAUX  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
# séparer R et T
rt_sg_lisse$wavelength <- as.numeric(as.character(rt_sg_lisse$wavelength))

reflectance <- rt_sg_lisse %>%
  dplyr::filter(propriete == 'reflectance')
# transmittance <- rt_sg_lisse %>%
#   dplyr::filter(propriete == 'transmittance')

# afficher premier échantillon
x11670314_r <- rt_sg_lisse %>%
  dplyr::filter(propriete == 'reflectance',
                sample_id == 11670314)
# x11670314_t <- rt_sg_lisse %>%
#   dplyr::filter(propriete == 'transmittance',
#                 sample_id == 11670314)

rt_long_lisse$wavelength <-  as.numeric(as.character(rt_long_lisse$wavelength))

x11670314_r_ORIGINAL <- rt_long_lisse %>%
  dplyr::filter(propriete == 'reflectance',
                sample_id == 11670314)
# x11670314_t_ORIGINAL <- rt_long_lisse %>%             # DONNÉES NON LISSÉES pour comparaison ******
#   dplyr::filter(propriete == 'transmittance',
#                 sample_id == 11670314)

# définir paramètres (régions de couleur, texte, etc)
regions <- data_frame(Region = factor(c('VIS-NIR', "NIR", 'SWIR1', 'SWIR2'),
                             levels = c('VIS-NIR', "NIR", 'SWIR1', 'SWIR2')),
                      xmin = c(250, 900, 1080, 1880),
                      xmax = c(900, 1080, 1880, 2700),
                      ymin = -0.1,
                      ymax = 1.2)
regions_text <- regions %>% dplyr::filter(Region != 'VIS-NIR1')       # "!="vis..." " dit de ne pas mettre texte "250" pour définir le début de cette région 
sgolay_text <- data_frame(x = c(560, 990, 1575, 2200), 
                          y = c(0.34, 0.34, 0.4, 0.34),
                          text = c('p = 3\nn = 21',
                                   'p = 3\nn = 55',
                                   'p = 3\nn = 21',
                                   # 'p = 3\nn = 21',
                                   'p = 4\nn = 175'))
library("RColorBrewer")

x11670314_plot_sg <- ggplot(x11670314_r) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Region), alpha = 0.55, data = regions) +
  geom_line(aes(x = wavelength, y = valeur_sg), linetype = 'solid', colour = 'black') + # valeur lissée réflectance
  geom_line(data = x11670314_r_ORIGINAL, aes(x = wavelength, y = valeur), linetype = 'dotted') +                     # valeur originale réflectance
  labs(x = 'Wavelength (nm)', y = 'Reflectance') +
  # labs(title = "Paramètres S-GOLAY sur données lissées et sans chevauchement (spectrolab::smooth() et match_sensors()), 2 juin 2019") +
  coord_cartesian(ylim = c(-0.025, 0.55), xlim = c(400, 2400), expand = F) +
  theme_bw() +
  theme(legend.position = "none", legend.text = element_text(face = "italic"),
        panel.border = element_rect(colour = brewer.pal(9, "Greys")[5], fill = NA, size = 1),
        panel.grid = element_blank()) +
  scale_fill_manual(values = brewer.pal(8, "Dark2")[c(1,2,5,6)]) +
  geom_text(aes(x = c(250, 880, 1085, 1880), y = c(0, 0.005, 0.005, 0.005), label = xmin, angle = 45), data = regions_text) +
  geom_text(aes(x = x, y = y, label = text), angle = 45, data = sgolay_text) # +
  # geom_vline(aes(xintercept = 2300), linetype = 'dashed')
x11670314_plot_sg
ggsave('~/Desktop/paramètres_S-GOLAY_lissées.pdf', plot = x11670314_plot_sg, height = 4, width = 4)


#=======================================# Format large avec nouvelles valeurs et Enregistrements #=======================================#                          ----

# ajout colonne des nouveaux noms de site : Intermediate (5N + 10N) et High
table(rt_sg_large$site)/2 # refl et trans = 2x
# GPB_canal    GTV_SE   MBP_10N   MBP_20N    MBP_5N  MBP_open 
#    24        24         3          11         9        23   

site_class_prov <- as.character(rt_sg_lisse$site)
site_class <- mgsub(site_class_prov, pattern = c("MBP_5N|MBP_10N", "MBP_20N"), 
                    replacement = c("MBP_Intermediate", "MBP_High")) # long

unique(site_class)  
rt_sg_lisse$site_class <- site_class

cols <- c("sample_id", "scientific_name", "site", "site_class", "propriete", "wavelength")
rt_sg_lisse[cols] <- lapply(rt_sg_lisse[cols] , as.factor)
levels(rt_sg_lisse$site_class)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Enregistrements
#-=-=-=-=-=-=-=-=-=-=-=-=-

# Enregistrement de rt_sg_lisse
# save(rt_sg_lisse, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_lisse.RData") # rt_sg_lisse_NO_doublon, fait 2 juin / refait 15 août ($site_class)/ refait 21 aout problèmes perdu !

# obtenir rt_sg_large
rt_sg_large <- dcast(rt_sg_lisse, sample_id + scientific_name + site + site_class + propriete ~ wavelength, 
                       value.var = "valeur_sg")

# vérification de la nouvelle colonne site_class
table(rt_sg_large$site_class)/2 # refl et trans = 2x
# GPB_canal           GTV_SE        MBP_High   MBP_Intermediate         MBP_open 
# 24               24               11               12               23 
# MBP_High -> 11 encore √, MBP_Intermediate = 12 = 9 + 3 √

# Enregistrement de rt_sg_large
# save(rt_sg_large, file = "~/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_large.RData") # rt_sg_large, fait 2 juin / refait 21 juin ! / refait 15 août (ajout $site_class) / refait 21 aout
# save(rt_sg_large, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_large.RData") # rt_sg_large, fait 2 juin / refait 21 juin ! / refait 15 août (ajout $site_class) / refait 21 aout

