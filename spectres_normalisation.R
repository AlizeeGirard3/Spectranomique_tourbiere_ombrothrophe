# 8 OCT 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                              NORMALISATION TOUS SPECTRES                              #  
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES :  NORMALIZATION : Feilhauer, Asner, Martin, & Schmidtlein, 2010 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(spectrolab)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# définir le répertoire
setwd("~/Documents/Maîtrise/DonnéesAnalyses/scripts-NETTOYAGE")

# load("~/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData")
# issu des manipulations dans ce script

# load("~/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_large.RData") # rt_sg_large
# obtenu via manips sur script "nettoyage_spectres_lissage_correction.R" dans le document "scripts-NETTOYAGE"


#=======================================# Manipulations préalables #=======================================#

#-=-=-=-=-=-=-=-=-=-=-=-=-
# objet "spectra" ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

refl <- rt_sg_large %>% 
  dplyr::filter(propriete == "reflectance")
trans <- rt_sg_large %>% 
  dplyr::filter(!propriete == "reflectance")

# objet spectra TOUS
refl_spobj <- as.spectra(refl, name_idx = 1, meta_idxs = c(2:5))
trans_spobj <- as.spectra(trans, name_idx = 1, meta_idxs = c(2:5))

# normaliser 
refl_spobj_norm <- normalize(refl_spobj)
trans_spobj_norm <- normalize(trans_spobj)
#   ******    √   Note that y values will not be true reflectances anymore!   ******

# remettre ensemble et arranger
refl_spobj_norm <- as.data.frame(refl_spobj_norm)
trans_spobj_norm <- as.data.frame(trans_spobj_norm)
rt_sg_normalized_large <- refl_spobj_norm %>% full_join(trans_spobj_norm) %>% 
  dplyr::select(sample_id = "sample_name", everything())
rt_sg_normalized_large$sample_id <- gsub( "spec_", "", rt_sg_normalized_large$sample_id) # enlever "spec_"

# enregistrer aux 3 endroits pertinents (Ça remplace les anciennes versions ! Attention !)
# save(rt_sg_normalized_large, file = "~/Documents/Maîtrise/DonnéesAnalyses/spectres/rt_sg_normalized_large.RData")
# save(rt_sg_normalized_large, file = "~/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData")
