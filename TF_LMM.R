# 21 nov 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                   VISUALISATION des LMM sur données de traits fonctionnels
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# notes du workshop : https://wiki.qcbs.ca/lmm.workshop
#                    https://docs.google.com/document/d/1bdgcstVYnxNRk9fbC1fWne0NCkz7jLATv3QPsvzzaVM/edit


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(nlme)
library(multcomp)
library("emmeans") # for posthoc tests
library(dplyr)
library(tidyverse) # rown_to_colmun
library("RColorBrewer")
library(plotrix) # std.error


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Importation données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

setwd("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/graphs&LMM")

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData")
# obtenu du script "LMM~tt_sites.R"

## manips préalables sur ft_all
ft_all <- ft_all_vX2 %>% 
  dplyr::select("site_class", "sample_id", "scientific_name", "plot_field_id", "specific_leaf_area_m2_kg",
                "leaf_mass_per_area_g_m2","leaf_dry_matter_content_mg_g","leaf_water_content_mg_g",
                "equivalent_water_thickness_cm", "N_pourc", "C_pourc", "hemicellulose", "cellulose", "lignine",
                "soluble_contents", "chlA_mg_g", "chlA_mg_m2", "chlB_mg_g", "chlB_mg_m2",
                "carotenoides_mg_g","C_N_ratio", "chlA_chlB") %>% # enlevé récalcitrants
  mutate_at(funs = is.character, vars(c(5:22)), as.numeric)
# bon noms pour les graphiques
ft_all$scientific_name = factor(ft_all$scientific_name, levels = c('Chamaedaphne calyculata', 'Eriophorum vaginatum',
                                                                   "Kalmia angustifolia", 'Rhododendron groenlandicum'),
                                labels = c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum"))
ft_all$scientific_name <- as.character(ft_all$scientific_name)
levels.sc <- c("C. calyculata","E. vaginatum","K. angustifolia", "R. groenlandicum") # ordre utilisé pour les graphiques

ft_all$site_class = factor(ft_all$site_class, levels = c("GTV_SE", "GPB_canal", "MBP_open", "MBP_Intermediate", "MBP_High"), 
                           labels = c("GTV", "GPB", "MBP.No.N", "MBP.Mid.N", "MBP.High.N")) 
ft_all$site_class <- as.character(ft_all$site_class)
levels.site <- c("MBP.No.N", "GTV", "GPB", "MBP.Mid.N", "MBP.High.N") # ordre utilisé pour les graphiques

## prevent scientific notation
options(scipen = 999)


#=-=-=-=-=-=-=-=-=-=-=-=-
# Mixed effect coefficients ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# ## create data summary for computed MLMs (methode REML par défaut)
# 
# # créer les tables de résultats vides
# # res.RDM.plus.anova <- list()
# res.RDM.mult.anova <- list()
# 
# res.FXD.anova <- list()
# 
# # res.FXD.plus.elements <- list()
# res.FXD.mult.elements <- list()
# 
# # liste avec tous les traits fonctionnels : tester l'interaction (trait ~ espèce * site) pour chaque trait
# inFT_df <- as.list(colnames(ft_all)[5:22], dimnames = list(inFT_df))
# 
# # liste avec traits fonctionnels ou AUCUNE interaction (trait ~ espèce + site), deux traits : chla/masse, chla/aire
# # # inFT_df.plus <- as.list(colnames(ft_all)[c(15:16)], dimnames = list(inFT_df))
# # # ft.plus <- ft_all[c(1:8, 15:16)]
# # # # DROPPER INTERACTION (pour elements et aussi pour tester RDM)  : chla masse, chla aire
# 
# # liste avec traits fonctionnels  ??????                              ######(trait ~ espèce + site), deux traits : chla/masse, chla/aire
# # inFT_df <- as.list(colnames(ft_all)[c(5:22)], dimnames = list(inFT_df))
# # ft.mult <- ft_all[c(1:22)]
# # # GARDER INTERACTION (pour elements et aussi pour tester RDM) : tous
# 
# # outmat.plus <- matrix(data = NA, nrow = length(inFT_df.plus), ncol = 4, dimnames = list(colnames(ft.plus)[9:10]))
# # outmat.random <- matrix(data = NA, nrow = length(inFT_df), ncol = 1, dimnames = list(colnames(ft_all[5:22])))
# outmat.int <- matrix(data = NA, nrow = length(inFT_df), ncol = 5, dimnames = list(colnames(ft_all)[5:22]))
# 
# 
# # # créer les df
# # for (i in 1:length(inFT_df.plus)){
# #   inFT_df.plus[[i]] <- as.data.frame(ft.plus[,c(8, 5, 7, i+8)]) # no interaction
# # }
# # for (i in 1:length(inFT_df)){
# #   inFT_df[[i]] <- as.data.frame(ft.mult[,c(1, 3, 4, i+4)]) # interaction  # [, 4] <- TF, [, 1:2] site class et sc, [,3] <- plot field id
# # }
# for (i in 1:length(inFT_df)){
#   inFT_df[[i]] <- as.data.frame(ft_all[,c(1, 3, 4, i+4)]) # tous # [, 4] <- TF, [, 1:2] site class et sc, [,3] <- plot field id
# }
# 
# 
# # ## compute p-values of RDM, intercept (≠ zéro), site/species effect and interaction
# # for(i in 1:length(inFT_df)) {
# #   print(i)
# # 
# #   # formule
# #   form.plus <- as.formula(paste(colnames(inFT_df[[i]][4]), "~", paste(colnames(inFT_df[[i]][c(1:2)]), collapse = "+"), sep = ""))
# # 
# #   # calculer p-value EFFET RANDOM (pas interaction)
# #   res.RDM.plus.anova[[i]] = anova(lme(form.plus, na.action = na.omit, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  method = "REML", data = inFT_df[[i]]),
# #                                   gls(form.plus, na.action = na.omit, weights = varIdent(form = ~ 1 |site_class), method = "REML", data = inFT_df[[i]]))
# #   outmat.random[i, 4] <- round(res.RDM.plus.anova[[i]]$`p-value`[2], 4)
# # 
# #   # calculer les groupes entre ÉLÉMENTS
# #   res.FXD.plus.elements[[i]] = anova(lme(form.plus, na.action = na.omit, random = ~ 1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class),  data = inFT_df[[i]])) # method = "REML" par défaut
# #   outmat.random[i, 1:3] <- format(round(res.FXD.plus.elements[[i]]$`p-value`, 4), digits = 4)
# # 
# # }
# # outmat.random <- as.data.frame(outmat.random)
# # # ajouter noms de TF (ligne)
# # outmat.random <- rownames_to_column(outmat.random)
# # colnames(outmat.random) <- c("FT", "Intercept ≠ zero", "Site effect", "Species effect", "Random effect")
# 
# for(i in 1:length(inFT_df)) {
#   print(i)
# 
#   # formule
#   form.mult <- as.formula(paste(colnames(inFT_df[[i]][4]), "~", paste(colnames(inFT_df[[i]][c(1:2)]), collapse = "*"), sep = ""))
# 
#   # calculer p-value EFFET RANDOM (interaction)
#   res.RDM.mult.anova[[i]] = anova(lme(form.mult, na.action = na.omit, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  method = "REML", data = inFT_df[[i]]),
#                                   gls(form.mult, na.action = na.omit, method = "REML", data = inFT_df[[i]]))
#   #outmat.mult[i, 5] <- format(round(res.RDM.mult.anova[[i]]$`p-value`[2], 5), digits = 5)
#   outmat.int[i, 5] <- round(res.RDM.mult.anova[[i]]$`p-value` [2], 4)
# 
#   # calculer les groupes entre ÉLÉMENTS 
#   res.FXD.mult.elements[[i]] = anova(lme(form.mult, na.action = na.omit, random = ~ 1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class),  data = inFT_df[[i]])) # method = "REML" par défaut
#   outmat.int[i, 1:4] <- format(round(res.FXD.mult.elements[[i]]$`p-value`, 4), digits = 4)
# 
# }
# outmat.int <- as.data.frame(outmat.int, stringsAsFactors = F)
# # ajouter noms de TF (ligne)
# outmat.int <- rownames_to_column(outmat.int)
# colnames(outmat.int) <- c("FT", "Intercept ≠ zero", "Site effect", "Species effect", "Site : Species", "Random effect")
# 
# # for(i in 1:length(inFT_df)) {
# #   print(i)
# # 
# #   # formules
# #   form.plus.all <- as.formula(paste(colnames(inFT_df[[i]][4]), "~", paste(colnames(inFT_df[[i]][c(1:2)]), collapse = "+"), sep = ""))
# #   form.mult.all <- as.formula(paste(colnames(inFT_df[[i]][4]), "~", paste(colnames(inFT_df[[i]][c(1:2)]), collapse = "*"), sep = ""))
# # 
# #   # calculer p-value INTERACTIONS
# #   res.FXD.anova[[i]] = anova(lme(form.mult.all, na.action = na.omit, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), method = "ML", data = inFT_df[[i]]),
# #                              lme(form.plus.all, na.action = na.omit, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), method = "ML", data = inFT_df[[i]]))
# #   outmat[i, 1] <- format(round(res.FXD.anova[[i]]$`p-value`[2], 4), digits = 4)
# # }
# # outmat <- as.data.frame(outmat, stringsAsFactors = F)
# # outmat <- rownames_to_column(outmat)
# # colnames(outmat) <- c("FT", "Interaction")
# 
# # rassembler df ensemble
# outmat.all <- outmat.int %>% #outmat.plus %>%
#   # full_join(outmat.random) %>%
#   dplyr::select( "Functional trait" = "FT", "Intercept ≠ zero", "Site effect",  # change order
#                  "Species effect", "Site : Species", "Random effect")
# 
# # outmat.all[, 5:6] <- mutate_if(outmat.all[, 5:6], is.factor, as.character)
# # outmat.all[, 5:6] <- mutate_if(outmat.all[, 5:6], is.character, as.numeric)
# 
# ## rename FT
# outmat.all$`Functional trait` <- factor(outmat.all$`Functional trait`, levels = c("specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2", "leaf_dry_matter_content_mg_g", "leaf_water_content_mg_g",
#                                                                                   "equivalent_water_thickness_cm", "N_pourc", "C_pourc", "hemicellulose", "cellulose", "lignine",
#                                                                                   "soluble_contents", "chlA_mg_g", "chlA_mg_m2", "chlB_mg_g", "chlB_mg_m2", "carotenoides_mg_g", "C_N_ratio", "chlA_chlB"),
#                                         labels = c('SLA (m2 kg-1)', 'LMA (g m-2)', 'LDMC (mg g-1)', 'LWC (mg g-1)', 'EWT (cm)',  'N (%)', 'Total C (%)', 'Hemicellulose (%)', 'Cellulose (%)',
#                                                    'Lignin (%)', 'Soluble carbon (%)', 'Chl a (mg g-1)', 'Chl a (mg m2)', 'Chl b (mg g-1)', 'Chl b (mg m2)',
#                                                    'Carotenoids (mg g-1)', 'C:N ratio',  'Chl a:b ratio'))
# 
# ## add asterix if significant
# for (i in 1:length(colnames(outmat.all[, 2:6]))) {
#   outmat.all[, i+1] <- ifelse(outmat.all[, i+1] < 0.05,
#                               paste0(format(outmat.all[, i+1], digits = 3), " *"),
#                               format(outmat.all[, i+1], digits = 3))
#   outmat.all[, i+1] <- ifelse(outmat.all[, i+1] < 0.01,
#                               paste0(format(outmat.all[, i+1], digits = 3), "*"),
#                               format(outmat.all[, i+1], digits = 3))
#   outmat.all[, i+1] <- ifelse(outmat.all[, i+1] < 0.001,
#                               paste("< 0.001 ***"),
#                               format(outmat.all[, i+1], digits = 3))
# }
# # save(outmat.all, file = "lmm.all.values.RData")
# # write.csv(outmat.all, "LMM.results.FT.all.csv") # fait 7 octobre/ refait 31 oct -> enlevé Recalcitrant


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Calculs pour visualisation ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

## joindre valeurs réelles (J'ai calculé les modèles linéaires, mais je veux afficher les valeurs moyennes
# par espèce et/ou par site, qui sont à calcluer.)
# 2 effets ≠ interaction (moyenne/site, puis calclul indépendant de moyenne/espèce)
ft_site <- ft_all %>% 
  dplyr::select("site_class", "specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2", "equivalent_water_thickness_cm",
                "chlA_mg_g", "chlA_mg_m2", "chlB_mg_g") %>% 
  group_by(site_class) %>%
  summarise_all(funs(mean, std.error), na.rm = TRUE) %>%
  ungroup()
ft_site <- unite(ft_site, col = "id", c(site_class))

ft_sc <- ft_all %>% 
  dplyr::select("scientific_name", "specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2", "equivalent_water_thickness_cm",
                "chlA_mg_g", "chlA_mg_m2", "chlB_mg_g") %>% 
  group_by(scientific_name) %>% 
  summarise_all(funs(mean, std.error), na.rm = TRUE) %>%
  ungroup()
ft_sc <- unite(ft_sc, col = "id", c(scientific_name))


## interaction (moyenne/site et par espèce, voir "group_by") POUR AFFICHAGE ####
ft_all_df <- ft_all %>% 
  dplyr::select("site_class", "scientific_name", "leaf_dry_matter_content_mg_g","leaf_water_content_mg_g",
                "N_pourc", "C_N_ratio", "C_pourc", "hemicellulose", "cellulose",  
                "soluble_contents", "chlB_mg_m2","carotenoides_mg_g") %>% 
  group_by(site_class, scientific_name) %>% 
  summarise_all(funs(mean, std.error), na.rm = TRUE) %>% 
  ungroup()

## interaction (moyenne/site et par espèce, voir "group_by") POUR GROUPEMENTS
# séparer les espèces (on veut voir l'effet site pour chaque espèce, indépendemment l'une de l'autre)
ft_all.CC <- ft_all %>% 
  dplyr::filter(scientific_name == "C. calyculata") 
ft_all.RG <- ft_all %>% 
  dplyr::filter(scientific_name == "R. groenlandicum")
ft_all.KA <- ft_all %>% 
  dplyr::filter(scientific_name == "K. angustifolia")
ft_all.EV <- ft_all %>% 
  dplyr::filter(scientific_name == "E. vaginatum")

CC <- matrix(data = NA, nrow = 5, ncol = 20) # groupe associé (lignes) // 10 traits (colonnes) + "groupe"
RG <- matrix(data = NA, nrow = 5, ncol = 20) 
EV <- matrix(data = NA, nrow = 5, ncol = 20) 
KA <- matrix(data = NA, nrow = 5, ncol = 20) 

# LME PAR ESPÈCE pour trouver les groupements
{ ####
  ## CC
  # nom des colonnes **EN ORDRE**
  colnames(CC) <- c("groupes", "N_pourc", "groupes", "C_N_ratio", "groupes", "C_pourc", "groupes", "hemicellulose", "groupes", "cellulose",
                    "groupes", "soluble_contents", "groupes", "chlB_mg_m2","groupes", "carotenoides_mg_g",
                    "groupes", "leaf_dry_matter_content_mg_g", "groupes", "leaf_water_content_mg_g")
  # N_pourc
  N_pourc.lme = lme(N_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.CC)
  N_pourc.tuckey <- emmeans(N_pourc.lme, pairwise ~ site_class, adjust = "tukey") #scientific_name * 
  N_pourc.groupes <- cld(N_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 1] <- as.vector(N_pourc.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 2] <- N_pourc.groupes$.group
  
  # "C_N_ratio"
  C_N.lme = lme(C_N_ratio ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  C_N.tuckey <- emmeans(C_N.lme, pairwise ~ site_class, adjust = "tukey")                   
  C_N.groupes <- cld(C_N.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 3] <- as.vector(C_N.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 4] <- C_N.groupes$.group
  
  # "C_pourc"  
  C_pourc.lme = lme(C_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.CC)
  C_pourc.tuckey <- emmeans(C_pourc.lme, pairwise ~ site_class, adjust = "tukey")
  C_pourc.groupes <- cld(C_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 5] <- as.vector(C_pourc.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 6] <- C_pourc.groupes$.group
  
  # "hemicellulose"
  Hcell.lme = lme(hemicellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  Hcell.tuckey <- emmeans(Hcell.lme, pairwise ~ site_class, adjust = "tukey")
  Hcell.groupes <- cld(Hcell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 7] <- as.vector(Hcell.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 8] <- Hcell.groupes$.group
  
  # "cellulose" 
  cell.lme = lme(cellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  cell.tuckey <- emmeans(cell.lme, pairwise ~ site_class, adjust = "tukey")
  cell.groupes <- cld(cell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 9] <- as.vector(cell.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 10] <- cell.groupes$.group
  
  # "soluble_contents"
  solC.lme = lme(soluble_contents ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  solC.tuckey <- emmeans(solC.lme, pairwise ~ site_class, adjust = "tukey")
  solC.groupes <- cld(solC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 11] <- as.vector(solC.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 12] <- solC.groupes$.group
  
  # "chlB_mg_m2"
  chlBaire.lme = lme(chlB_mg_m2 ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  chlBaire.tuckey <- emmeans(chlBaire.lme, pairwise ~ site_class, adjust = "tukey")
  chlBaire.groupes <- cld(chlBaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 13] <- as.vector(chlBaire.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 14] <- chlBaire.groupes$.group
  
  # "carotenoides_mg_g"
  carot.lme = lme(carotenoides_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.CC)
  carot.tuckey <- emmeans(carot.lme, pairwise ~ site_class * site_class, adjust = "tukey")
  carot.groupes <- cld(carot.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 15] <- as.vector(carot.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 16] <- carot.groupes$.group
  
  # "leaf_dry_matter_content_mg_g"
  LDMC.lme <- lme(leaf_dry_matter_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.CC) 
  LDMC.tuckey <- emmeans(LDMC.lme, pairwise ~ site_class, adjust = "tukey")
  LDMC.groupes <- cld(LDMC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 17] <- as.vector(LDMC.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 18] <- LDMC.groupes$.group
  
  # "leaf_water_content_mg_g"
  LWC.lme <- lme(leaf_water_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.CC) 
  LWC.tuckey <- emmeans(LWC.lme, pairwise ~ site_class, adjust = "tukey")
  LWC.groupes <- cld(LWC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  CC[1:5, 19] <- as.vector(LWC.groupes$site_class) # ajouter nom des sites comparés
  CC[1:5, 20] <- LWC.groupes$.group
  
  write.csv(CC, "CC.group.csv") # fait 22 décembre
}
{ ####
  ## RG
  # nom des colonnes **EN ORDRE**
  colnames(RG) <-  c("groupes", "N_pourc", "groupes", "C_N_ratio", "groupes", "C_pourc", "groupes", "hemicellulose", "groupes", "cellulose",
                     "groupes", "soluble_contents", "groupes", "chlB_mg_m2","groupes", "carotenoides_mg_g",
                     "groupes", "leaf_dry_matter_content_mg_g", "groupes", "leaf_water_content_mg_g")
  # N_pourc
  N_pourc.lme = lme(N_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.RG)
  N_pourc.tuckey <- emmeans(N_pourc.lme, pairwise ~ site_class, adjust = "tukey") #scientific_name * 
  N_pourc.groupes <- cld(N_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 1] <- as.vector(N_pourc.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 2] <- N_pourc.groupes$.group
  
  # "C_N_ratio"
  C_N.lme = lme(C_N_ratio ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  C_N.tuckey <- emmeans(C_N.lme, pairwise ~ site_class, adjust = "tukey")                   
  C_N.groupes <- cld(C_N.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 3] <- as.vector(C_N.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 4] <- C_N.groupes$.group
  
  # "C_pourc"  
  C_pourc.lme = lme(C_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.RG)
  C_pourc.tuckey <- emmeans(C_pourc.lme, pairwise ~ site_class, adjust = "tukey")
  C_pourc.groupes <- cld(C_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 5] <- as.vector(C_pourc.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 6] <- C_pourc.groupes$.group
  
  # "hemicellulose"
  Hcell.lme = lme(hemicellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  Hcell.tuckey <- emmeans(Hcell.lme, pairwise ~ site_class, adjust = "tukey")
  Hcell.groupes <- cld(Hcell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 7] <- as.vector(Hcell.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 8] <- Hcell.groupes$.group
  
  # "cellulose" 
  cell.lme = lme(cellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  cell.tuckey <- emmeans(cell.lme, pairwise ~ site_class, adjust = "tukey")
  cell.groupes <- cld(cell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 9] <- as.vector(cell.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 10] <- cell.groupes$.group
  
  # "soluble_contents"
  solC.lme = lme(soluble_contents ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  solC.tuckey <- emmeans(solC.lme, pairwise ~ site_class, adjust = "tukey")
  solC.groupes <- cld(solC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 11] <- as.vector(solC.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 12] <- solC.groupes$.group
  
  # "chlB_mg_m2"
  chlBaire.lme = lme(chlB_mg_m2 ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  chlBaire.tuckey <- emmeans(chlBaire.lme, pairwise ~ site_class, adjust = "tukey")
  chlBaire.groupes <- cld(chlBaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 13] <- as.vector(chlBaire.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 14] <- chlBaire.groupes$.group
  
  # "carotenoides_mg_g"
  carot.lme = lme(carotenoides_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.RG)
  carot.tuckey <- emmeans(carot.lme, pairwise ~ site_class * site_class, adjust = "tukey")
  carot.groupes <- cld(carot.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 15] <- as.vector(carot.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 16] <- carot.groupes$.group
  
  # "leaf_dry_matter_content_mg_g"
  LDMC.lme <- lme(leaf_dry_matter_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.RG) 
  LDMC.tuckey <- emmeans(LDMC.lme, pairwise ~ site_class, adjust = "tukey")
  LDMC.groupes <- cld(LDMC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 17] <- as.vector(LDMC.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 18] <- LDMC.groupes$.group
  
  # "leaf_water_content_mg_g"
  LWC.lme <- lme(leaf_water_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.RG) 
  LWC.tuckey <- emmeans(LWC.lme, pairwise ~ site_class, adjust = "tukey")
  LWC.groupes <- cld(LWC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  RG[1:5, 19] <- as.vector(LWC.groupes$site_class) # ajouter nom des sites comparés
  RG[1:5, 20] <- LWC.groupes$.group
  
  write.csv(RG, "RG.group.csv") # fait 22 décembre
}
{ ####
  ## EV
  # nom des colonnes **EN ORDRE**
  colnames(EV) <-  c("groupes", "N_pourc", "groupes", "C_N_ratio", "groupes", "C_pourc", "groupes", "hemicellulose", "groupes", "cellulose",
                     "groupes", "soluble_contents", "groupes", "chlB_mg_m2","groupes", "carotenoides_mg_g",
                     "groupes", "leaf_dry_matter_content_mg_g", "groupes", "leaf_water_content_mg_g")
  # N_pourc
  N_pourc.lme = lme(N_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.EV)
  N_pourc.tuckey <- emmeans(N_pourc.lme, pairwise ~ site_class, adjust = "tukey") #scientific_name * 
  N_pourc.groupes <- cld(N_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 1] <- as.vector(N_pourc.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 2] <- N_pourc.groupes$.group
  
  # "C_N_ratio"
  C_N.lme = lme(C_N_ratio ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  C_N.tuckey <- emmeans(C_N.lme, pairwise ~ site_class, adjust = "tukey")                   
  C_N.groupes <- cld(C_N.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 3] <- as.vector(C_N.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 4] <- C_N.groupes$.group
  
  # "C_pourc"  
  C_pourc.lme = lme(C_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.EV)
  C_pourc.tuckey <- emmeans(C_pourc.lme, pairwise ~ site_class, adjust = "tukey")
  C_pourc.groupes <- cld(C_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 5] <- as.vector(C_pourc.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 6] <- C_pourc.groupes$.group
  
  # "hemicellulose"
  Hcell.lme = lme(hemicellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  Hcell.tuckey <- emmeans(Hcell.lme, pairwise ~ site_class, adjust = "tukey")
  Hcell.groupes <- cld(Hcell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 7] <- as.vector(Hcell.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 8] <- Hcell.groupes$.group
  
  # "cellulose" 
  cell.lme = lme(cellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  cell.tuckey <- emmeans(cell.lme, pairwise ~ site_class, adjust = "tukey")
  cell.groupes <- cld(cell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 9] <- as.vector(cell.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 10] <- cell.groupes$.group
  
  # "soluble_contents"
  solC.lme = lme(soluble_contents ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  solC.tuckey <- emmeans(solC.lme, pairwise ~ site_class, adjust = "tukey")
  solC.groupes <- cld(solC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 11] <- as.vector(solC.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 12] <- solC.groupes$.group
  
  # "chlB_mg_m2"
  chlBaire.lme = lme(chlB_mg_m2 ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  chlBaire.tuckey <- emmeans(chlBaire.lme, pairwise ~ site_class, adjust = "tukey")
  chlBaire.groupes <- cld(chlBaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 13] <- as.vector(chlBaire.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 14] <- chlBaire.groupes$.group
  
  # "carotenoides_mg_g"
  carot.lme = lme(carotenoides_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.EV)
  carot.tuckey <- emmeans(carot.lme, pairwise ~ site_class * site_class, adjust = "tukey")
  carot.groupes <- cld(carot.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 15] <- as.vector(carot.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 16] <- carot.groupes$.group
  
  # "leaf_dry_matter_content_mg_g"
  LDMC.lme <- lme(leaf_dry_matter_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.EV) 
  LDMC.tuckey <- emmeans(LDMC.lme, pairwise ~ site_class, adjust = "tukey")
  LDMC.groupes <- cld(LDMC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 17] <- as.vector(LDMC.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 18] <- LDMC.groupes$.group
  
  # "leaf_water_content_mg_g"
  LWC.lme <- lme(leaf_water_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.EV) 
  LWC.tuckey <- emmeans(LWC.lme, pairwise ~ site_class, adjust = "tukey")
  LWC.groupes <- cld(LWC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  EV[1:5, 19] <- as.vector(LWC.groupes$site_class) # ajouter nom des sites comparés
  EV[1:5, 20] <- LWC.groupes$.group
  
  write.csv(EV, "EV.group.csv") # fait 22 décembre
}
{ ####
  ## KA
  # nom des colonnes **EN ORDRE**
  colnames(KA) <- c("groupes", "N_pourc", "groupes", "C_N_ratio", "groupes", "C_pourc", "groupes", "hemicellulose", "groupes", "cellulose",
                    "groupes", "soluble_contents", "groupes", "chlB_mg_m2","groupes", "carotenoides_mg_g",
                    "groupes", "leaf_dry_matter_content_mg_g", "groupes", "leaf_water_content_mg_g")
  # N_pourc
  N_pourc.lme = lme(N_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.KA)
  N_pourc.tuckey <- emmeans(N_pourc.lme, pairwise ~ site_class, adjust = "tukey") #scientific_name * 
  N_pourc.groupes <- cld(N_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 1] <- as.vector(N_pourc.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 2] <- N_pourc.groupes$.group
  
  # "C_N_ratio"
  C_N.lme = lme(C_N_ratio ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.KA)
  C_N.tuckey <- emmeans(C_N.lme, pairwise ~ site_class, adjust = "tukey")                   
  C_N.groupes <- cld(C_N.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 3] <- as.vector(C_N.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 4] <- C_N.groupes$.group
  
  # "C_pourc"  
  C_pourc.lme = lme(C_pourc ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, method = "REML", data = ft_all.KA)
  C_pourc.tuckey <- emmeans(C_pourc.lme, pairwise ~ site_class, adjust = "tukey")
  C_pourc.groupes <- cld(C_pourc.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 5] <- as.vector(C_pourc.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 6] <- C_pourc.groupes$.group
  
  # "hemicellulose"
  Hcell.lme = lme(hemicellulose ~ site_class, random = ~1|plot_field_id,  na.action = na.omit, method = "REML", data = ft_all.KA)
  Hcell.tuckey <- emmeans(Hcell.lme, pairwise ~ site_class, adjust = "tukey")
  Hcell.groupes <- cld(Hcell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 7] <- as.vector(Hcell.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 8] <- Hcell.groupes$.group
  
  # "cellulose" 
  cell.lme = lme(cellulose ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.KA)
  cell.tuckey <- emmeans(cell.lme, pairwise ~ site_class, adjust = "tukey")
  cell.groupes <- cld(cell.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 9] <- as.vector(cell.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 10] <- cell.groupes$.group
  
  # "soluble_contents"
  solC.lme = lme(soluble_contents ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.KA)
  solC.tuckey <- emmeans(solC.lme, pairwise ~ site_class, adjust = "tukey")
  solC.groupes <- cld(solC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 11] <- as.vector(solC.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 12] <- solC.groupes$.group
  
  # "chlB_mg_m2"
  chlBaire.lme = lme(chlB_mg_m2 ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.KA)
  chlBaire.tuckey <- emmeans(chlBaire.lme, pairwise ~ site_class, adjust = "tukey")
  chlBaire.groupes <- cld(chlBaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 13] <- as.vector(chlBaire.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 14] <- chlBaire.groupes$.group
  
  # "carotenoides_mg_g"
  carot.lme = lme(carotenoides_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all.KA)
  carot.tuckey <- emmeans(carot.lme, pairwise ~ site_class, adjust = "tukey")
  carot.groupes <- cld(carot.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 15] <- as.vector(carot.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 16] <- carot.groupes$.group
  
  # "leaf_dry_matter_content_mg_g"
  LDMC.lme <- lme(leaf_dry_matter_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.KA) 
  LDMC.tuckey <- emmeans(LDMC.lme, pairwise ~ site_class, adjust = "tukey")
  LDMC.groupes <- cld(LDMC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 17] <- as.vector(LDMC.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 18] <- LDMC.groupes$.group
  
  # "leaf_water_content_mg_g"
  LWC.lme <- lme(leaf_water_content_mg_g ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all.KA) 
  LWC.tuckey <- emmeans(LWC.lme, pairwise ~ site_class, adjust = "tukey")
  LWC.groupes <- cld(LWC.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
  # compiler dans tableau de résultats
  KA[1:5, 19] <- as.vector(LWC.groupes$site_class) # ajouter nom des sites comparés
  KA[1:5, 20] <- LWC.groupes$.group
  
  write.csv(KA, "KA.group.csv") # fait 22 décembre
}

# SANS effet site (on élimine la valeur de site) ####
ft_all_df_NO.eff_site <- ft_all %>%  
  dplyr::select('scientific_name', "lignine") %>% 
  group_by(scientific_name) %>% 
  summarise_all(funs(mean, std.error), na.rm = TRUE) %>% 
  ungroup()

# Unir les colonnes sci.name&site pour joindre les tables de LMM et les tables de valeurs moyenne mesurées pour graphiques
ft_all_df <- unite(ft_all_df, col = "id", c(scientific_name, site_class))


#  Calcul des modèles sur les traits pour en extraire les groupements (affiché dans les graphiques) ####

# 2 effets, pas d'interaction, pour les traits fonctionnels suivants :
# "specific_leaf_area_m2_kg", "leaf_mass_per_area_g_m2", "equivalent_water_thickness_cm",
# "chlA_mg_g", "chlA_mg_m2", "chlB_mg_g"

# SLA - SITE
# afficher un graph
SLA <- ggplot(SLA.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = SLA.mult.site$specific_leaf_area_m2_kg_mean, color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(SLA.mult.site, mapping = aes(
    ymin = specific_leaf_area_m2_kg_mean - specific_leaf_area_m2_kg_std.error,
    ymax = specific_leaf_area_m2_kg_mean + specific_leaf_area_m2_kg_std.error), width = 0) +
  # ajouter à la main *** / mentionner dans le texte : "/espèce" / arranger éventuellement
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (specific_leaf_area_m2_kg_mean + specific_leaf_area_m2_kg_std.error + 0.3)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("SLA (m" ^ 2 ~"kg" ^ -1~ ")"), x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
SLA
#ggplot2::ggsave('~/Desktop/LMM_SLA.site.png', SLA, width = 2.7, height = 2.4)

# ESPÈCE
SLA.lme = lme(specific_leaf_area_m2_kg ~ scientific_name, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |scientific_name),  na.action = na.omit, method = "REML", data = ft_all)
SLA.tuckey <- emmeans(SLA.lme, pairwise ~ scientific_name, adjust = "tukey")
SLA.mult <- cld(SLA.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
SLA.mult.prov <- unite(SLA.mult, col = "id", c(scientific_name))
SLA.mult.sc = full_join(SLA.mult.prov, ft_sc[, c('id', "specific_leaf_area_m2_kg_mean", "specific_leaf_area_m2_kg_std.error")], by = "id")
SLA.mult.sc <-  separate(SLA.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

SLA <- ggplot(SLA.mult.sc, mapping = aes(x = scientific_name, y = SLA.mult.sc$specific_leaf_area_m2_kg_mean, color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(SLA.mult.sc, mapping = aes(
    ymin = specific_leaf_area_m2_kg_mean - specific_leaf_area_m2_kg_std.error,
    ymax = specific_leaf_area_m2_kg_mean + specific_leaf_area_m2_kg_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (specific_leaf_area_m2_kg_mean + specific_leaf_area_m2_kg_std.error + 0.3)),
            position = position_nudge(x = c(0.01,0,0,-0.07))) +  # EV, RG, CC, KA
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("SLA (m" ^ 2 ~"kg" ^ -1~ ")"), x = "Species") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines")) #, 
#strip.text.x = element_text(face = "italic"))
SLA
#ggplot2::ggsave('~/Desktop/LMM_SLA.sc.png', SLA, width = 2.7, height = 2.95)


# LMA  - SITE
LMA.lme <- lme(leaf_mass_per_area_g_m2 ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all)  
LMA.tuckey <- emmeans(LMA.lme, pairwise ~ site_class, adjust = "tukey")
LMA.mult <- cld(LMA.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
LMA.mult.prov <- unite(LMA.mult, col = "id", site_class)
LMA.mult.site = full_join(LMA.mult.prov, ft_site[, c('id', "leaf_mass_per_area_g_m2_mean", "leaf_mass_per_area_g_m2_std.error")], by = "id")
LMA.mult.site <-  separate(LMA.mult.site, id, "site_class", sep = "_", remove = TRUE)

LMA <- ggplot(LMA.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = leaf_mass_per_area_g_m2_mean, 
                                           color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(LMA.mult.site, mapping = aes(
    ymin = leaf_mass_per_area_g_m2_mean - leaf_mass_per_area_g_m2_std.error,
    ymax = leaf_mass_per_area_g_m2_mean + leaf_mass_per_area_g_m2_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (leaf_mass_per_area_g_m2_mean + leaf_mass_per_area_g_m2_std.error + 4)),
            position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("LMA (g" ~ "m" ^ -2 ~")"), x = "Site") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none",  
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
LMA
#ggplot2::ggsave('~/Desktop/LMM_LMA.site.png', LMA, width = 2.7, height = 2.4)
# #ggplot2::ggsave('~/Desktop/LMM_LÉGENDE.site.png', LMA, width = 2.7, height = 2.4)

# LMA  - ESPÈCE
LMA.lme <- lme(leaf_mass_per_area_g_m2 ~ scientific_name, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |scientific_name), na.action = na.omit, data = ft_all)  
LMA.tuckey <- emmeans(LMA.lme, pairwise ~ scientific_name, adjust = "tukey")
LMA.mult <- cld(LMA.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
LMA.mult.prov <- unite(LMA.mult, col = "id", scientific_name)
LMA.mult.sc = full_join(LMA.mult.prov, ft_sc[, c('id', "leaf_mass_per_area_g_m2_mean", "leaf_mass_per_area_g_m2_std.error")], by = "id")
LMA.mult.sc <-  separate(LMA.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

LMA <- ggplot(LMA.mult.sc, mapping = aes(x = scientific_name, y = LMA.mult.sc$leaf_mass_per_area_g_m2_mean, color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(LMA.mult.sc, mapping = aes(
    ymin = leaf_mass_per_area_g_m2_mean - leaf_mass_per_area_g_m2_std.error,
    ymax = leaf_mass_per_area_g_m2_mean + leaf_mass_per_area_g_m2_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (leaf_mass_per_area_g_m2_mean + leaf_mass_per_area_g_m2_std.error + 4)),
            position = position_nudge(x = c(0.02,0,-0.1,-0.1))) +  # KA, CC, RG, EV  - vers GAUCHE, + vers DROITE
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 16, 16)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("LMA (g" ~ "m" ^ -2 ~")"), x = "Species") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), legend.position = "none",  
        panel.border = element_rect(0), panel.spacing = unit(0, "lines")) #, 
# legend.text = element_text(face = "italic"))
LMA
#ggplot2::ggsave('~/Desktop/LMM_LMA.sc.png', LMA, width = 2.7, height = 2.65)
# #ggplot2::ggsave('~/Desktop/LMM_LÉGENDE.sc.png', LMA, width = 2.7, height = 2.4)


# EWT  - SITE
EWT.lme <- lme(equivalent_water_thickness_cm ~ site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all) 
EWT.tuckey <- emmeans(EWT.lme, pairwise ~ site_class, adjust = "tukey")
EWT.mult <- cld(EWT.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
EWT.mult.prov <- unite(EWT.mult, col = "id", site_class)
EWT.mult.site = full_join(EWT.mult.prov, ft_site[, c('id', "equivalent_water_thickness_cm_mean", "equivalent_water_thickness_cm_std.error")], by = "id")
EWT.mult.site <-  separate(EWT.mult.site, id, "site_class", sep = "_", remove = TRUE)

EWT <- ggplot(EWT.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = equivalent_water_thickness_cm_mean, 
                                           color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(EWT.mult.site, mapping = aes(
    ymin = equivalent_water_thickness_cm_mean - equivalent_water_thickness_cm_std.error,
    ymax = equivalent_water_thickness_cm_mean + equivalent_water_thickness_cm_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (equivalent_water_thickness_cm_mean + equivalent_water_thickness_cm_std.error + 0.0006)),
            position = position_nudge(x = c(-0.04,-0.04,-0.04,-0.04,-0.04))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = "ETW (cm)", x = "Site") + # ~ "Equivalent water thickness (ETW)" 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
EWT
#ggplot2::ggsave('~/Desktop/LMM_EWT.site.png', EWT, width = 2.7, height = 2.4)

# EWT  - ESPÈCE
EWT.lme <- lme(equivalent_water_thickness_cm ~ scientific_name, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |scientific_name), na.action = na.omit, data = ft_all) 
EWT.tuckey <- emmeans(EWT.lme, pairwise ~ scientific_name, adjust = "tukey")
EWT.mult <- cld(EWT.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
EWT.mult.prov <- unite(EWT.mult, col = "id", scientific_name)
EWT.mult.sc = full_join(EWT.mult.prov, ft_sc[, c('id', "equivalent_water_thickness_cm_mean", "equivalent_water_thickness_cm_std.error")], by = "id")
EWT.mult.sc <-  separate(EWT.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

EWT <- ggplot(EWT.mult.sc, mapping = aes(x = scientific_name, y = equivalent_water_thickness_cm_mean, 
                                         color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(EWT.mult.sc, mapping = aes(
    ymin = equivalent_water_thickness_cm_mean - equivalent_water_thickness_cm_std.error,
    ymax = equivalent_water_thickness_cm_mean + equivalent_water_thickness_cm_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (equivalent_water_thickness_cm_mean + equivalent_water_thickness_cm_std.error + 0.001)),
            position = position_nudge(x = c(0.03,0,0,-0.09))) +
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = "ETW (cm)", x = "Species") + # ~ "Equivalent water thickness (ETW)" 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines")) #, 
# strip.text.x = element_text(face = "italic"))
EWT
#ggplot2::ggsave('~/Desktop/LMM_EWT.sc.png', EWT, width = 2.7, height = 2.65)


# CHLA, masse  - SITE
chlA.lme = lme(chlA_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all)
chlA.tuckey <- emmeans(chlA.lme, pairwise ~ site_class, adjust = "tukey") 
chlA.mult <- cld(chlA.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlA.mult.prov <- unite(chlA.mult, col = "id", site_class)
chlA.mult.site = full_join(chlA.mult.prov, ft_site[, c('id', "chlA_mg_g_mean", "chlA_mg_g_std.error")], by = "id")
chlA.mult.site <-  separate(chlA.mult.site, id, "site_class", sep = "_", remove = TRUE)

chlA <- ggplot(chlA.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = chlA_mg_g_mean, 
                                             color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(chlA.mult.site, mapping = aes(
    ymin = chlA_mg_g_mean - chlA_mg_g_std.error,
    ymax = chlA_mg_g_mean + chlA_mg_g_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (chlA_mg_g_mean + chlA_mg_g_std.error + 0.5)),
            position = position_nudge(x = c(-0.04,-0.04,-0.04,-0.04,-0.04))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Chl" ~ italic("a") ~ "(mg g" ^"-1"~")"), x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlA
#ggplot2::ggsave('~/Desktop/LMM_chlA.site.png', chlA, width = 2.7, height = 2.4)

# CHLA, masse  - ESPÈCE
chlA.lme = lme(chlA_mg_g ~ scientific_name, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |scientific_name),  na.action = na.omit, method = "REML", data = ft_all)
chlA.tuckey <- emmeans(chlA.lme, pairwise ~ scientific_name, adjust = "tukey") 
chlA.mult <- cld(chlA.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlA.mult.prov <- unite(chlA.mult, col = "id", scientific_name)
chlA.mult.sc = full_join(chlA.mult.prov, ft_sc[, c('id', "chlA_mg_g_mean", "chlA_mg_g_std.error")], by = "id")
chlA.mult.sc <-  separate(chlA.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

chlA <- ggplot(chlA.mult.sc, mapping = aes(x = scientific_name, y = chlA_mg_g_mean, color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(chlA.mult.sc, mapping = aes(
    ymin = chlA_mg_g_mean - chlA_mg_g_std.error,
    ymax = chlA_mg_g_mean + chlA_mg_g_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (chlA_mg_g_mean + chlA_mg_g_std.error + 0.5)),
            position = position_nudge(x = c(0.03,0,-0.01,-0.06))) + # RG, CC, KA, EV, + DROITE
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("Chl" ~ italic("a") ~ "(mg g" ^"-1"~")"), x = "Species") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlA
#ggplot2::ggsave('~/Desktop/LMM_chlA.sc.png', chlA, width = 2.7, height = 2.65)


# CHLA, aire  - SITE
chlAaire.lme = lme(chlA_mg_m2 ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all)
chlAaire.tuckey <- emmeans(chlAaire.lme, pairwise ~ site_class, adjust = "tukey") 
chlAaire.mult <- cld(chlAaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlAaire.mult.prov <- unite(chlAaire.mult, col = "id", site_class)
chlAaire.mult.site = full_join(chlAaire.mult.prov, ft_site[, c('id', "chlA_mg_m2_mean", "chlA_mg_m2_std.error")], by = "id")
chlAaire.mult.site <-  separate(chlAaire.mult.site, id, "site_class", sep = "_", remove = TRUE)

chlAaire <- ggplot(chlAaire.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = chlA_mg_m2_mean, color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(chlAaire.mult.site, mapping = aes(
    ymin = chlA_mg_m2_mean - chlA_mg_m2_std.error,
    ymax = chlA_mg_m2_mean + chlA_mg_m2_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (chlA_mg_m2_mean + chlA_mg_m2_std.error + 14)),
            position = position_nudge(x = c(-0.04,-0.04,-0.04,-0.04,-0.04))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Chl" ~ italic("a") ~ "(mg m"^"-2"~")"), x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlAaire
#ggplot2::ggsave('~/Desktop/LMM_chlAaire.site.png', chlAaire, width = 2.7, height = 2.4)

# CHLA, aire  - ESPÈCE
chlAaire.lme = lme(chlA_mg_m2 ~ scientific_name, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |scientific_name),  na.action = na.omit, method = "REML", data = ft_all)
chlAaire.tuckey <- emmeans(chlAaire.lme, pairwise ~ scientific_name, adjust = "tukey") 
chlAaire.mult <- cld(chlAaire.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlAaire.mult.prov <- unite(chlAaire.mult, col = "id", scientific_name)
chlAaire.mult.sc = full_join(chlAaire.mult.prov, ft_sc[, c('id', "chlA_mg_m2_mean", "chlA_mg_m2_std.error")], by = "id")
chlAaire.mult.sc <-  separate(chlAaire.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

chlAaire <- ggplot(chlAaire.mult.sc, mapping = aes(x = scientific_name, y = chlA_mg_m2_mean, color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(chlAaire.mult.sc, mapping = aes(
    ymin = chlA_mg_m2_mean - chlA_mg_m2_std.error,
    ymax = chlA_mg_m2_mean + chlA_mg_m2_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (chlA_mg_m2_mean + chlA_mg_m2_std.error + 30)),
            position = position_nudge(x = c(0,0,0,-0.05))) +
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("Chl" ~ italic("a") ~ "(mg m" ^"-2"~")"), x = "Species") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlAaire
#ggplot2::ggsave('~/Desktop/LMM_chlAaire.sc.png', chlAaire, width = 2.7, height = 2.65)


# CHLB, masse  - SITE
chlB.lme = lme(chlB_mg_g ~ site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all)
chlB.tuckey <- emmeans(chlB.lme, pairwise ~ site_class, adjust = "tukey") # all : 0.0209
chlB.mult <- cld(chlB.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlB.mult.prov <- unite(chlB.mult, col = "id", site_class)
chlB.mult.site = full_join(chlB.mult.prov, ft_site[, c('id', "chlB_mg_g_mean", "chlB_mg_g_std.error")], by = "id")
chlB.mult.site <-  separate(chlB.mult.site, id, "site_class", sep = "_", remove = TRUE)

chlB <- ggplot(chlB.mult.site, mapping = aes(x = factor(site_class, levels = levels.site), y = chlB_mg_g_mean, 
                                             color = site_class)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(chlB.mult.site, mapping = aes(
    ymin = chlB_mg_g_mean - chlB_mg_g_std.error,
    ymax = chlB_mg_g_mean + chlB_mg_g_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (chlB_mg_g_mean + chlB_mg_g_std.error + 0.1)),
            position = position_nudge(x = c(-0.04,-0.04,-0.04,-0.04,-0.04))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Chl" ~ italic("b") ~ "(mg g" ^"-1"~")"), x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlB
#ggplot2::ggsave('~/Desktop/LMM_chlB.site.png', chlB, width = 2.7, height = 2.4)

# CHLB, masse  - ESPÈCE
chlB.lme = lme(chlB_mg_g ~ scientific_name, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |scientific_name),  na.action = na.omit, method = "REML", data = ft_all)
chlB.tuckey <- emmeans(chlB.lme, pairwise ~ scientific_name, adjust = "tukey")
chlB.mult <- cld(chlB.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
chlB.mult.prov <- unite(chlB.mult, col = "id", scientific_name)
chlB.mult.sc = full_join(chlB.mult.prov, ft_sc[, c('id', "chlB_mg_g_mean", "chlB_mg_g_std.error")], by = "id")
chlB.mult.sc <-  separate(chlB.mult.sc, id, "scientific_name", sep = "_", remove = TRUE)

chlB <- ggplot(chlB.mult.sc, mapping = aes(x = scientific_name, y = chlB_mg_g_mean, color = scientific_name)) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = scientific_name)) +
  geom_errorbar(chlB.mult.sc, mapping = aes(
    ymin = chlB_mg_g_mean - chlB_mg_g_std.error,
    ymax = chlB_mg_g_mean + chlB_mg_g_std.error), width = 0) +
  geom_text(mapping = aes(x = scientific_name, label = .group,
                          y = (chlB_mg_g_mean + chlB_mg_g_std.error + 0.12)),
            position = position_nudge(x = c(0.01,0,0,-0.09))) +
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("Chl" ~ italic("b") ~ "(mg g" ^"-1"~")"), x = "Species") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
chlB
#ggplot2::ggsave('~/Desktop/LMM_chlB.sc.png', chlB, width = 2.7, height = 2.65)


# 2 effets + interaction (*) pour les traits fonctionnels suivants :  ####
# "leaf_dry_matter_content_mg_g","leaf_water_content_mg_g",
# "N_pourc", "C_N_ratio", "C_pourc", "hemicellulose", "cellulose",  
# "soluble_contents", "chlB_mg_m2","carotenoides_mg_g"

# LDMC
LDMC <- ggplot(LDMC.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = leaf_dry_matter_content_mg_g_mean, 
                                        color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(LDMC.mult, mapping = aes(
    ymin = leaf_dry_matter_content_mg_g_mean - leaf_dry_matter_content_mg_g_std.error,
    ymax = leaf_dry_matter_content_mg_g_mean + leaf_dry_matter_content_mg_g_std.error), width = 0) +
  geom_text(mapping = aes(x = site_class, label = .group,
                          y = (leaf_dry_matter_content_mg_g_mean + leaf_dry_matter_content_mg_g_std.error + 5)),
            position = position_nudge(x = c(-0.1,-0.1,-0.1,-0.1,-0.06)), angle = 45) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = "LDMC (mg g" ^ -1 ~")", x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
LDMC
#ggplot2::ggsave('~/Desktop/LMM_LDMC.png', LDMC, width = 4.7, height = 2.4) #, units = "mm")


# LWC
LWC <- ggplot(LWC.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = leaf_water_content_mg_g_mean, 
                                      color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(LWC.mult, mapping = aes(
    ymin = leaf_water_content_mg_g_mean - leaf_water_content_mg_g_std.error,
    ymax = leaf_water_content_mg_g_mean + leaf_water_content_mg_g_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (leaf_water_content_mg_g_mean + leaf_water_content_mg_g_std.error + 5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = "LWC (mg g" ^ -1 ~")", x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
LWC
#ggplot2::ggsave('~/Desktop/LMM_LWC.png', LWC, width = 4.7, height = 2.4)


# N %
N_pourc <- ggplot(N_pourc.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = N_pourc_mean, color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(N_pourc.mult, mapping = aes(
    ymin = N_pourc_mean - N_pourc_std.error,
    ymax = N_pourc_mean + N_pourc_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (N_pourc_mean + N_pourc_std.error + 0.05)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = "N (%)", x = "Site") + # ~ "Nitrogen concentration"
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
N_pourc
#ggplot2::ggsave('~/Desktop/LMM_N_pourc.png', N_pourc, width = 4.7, height = 2.5)


# C %
C_pourc <- ggplot(ft_all_df, mapping = aes(x = factor(site_class, levels = levels.site), y = C_pourc_mean, 
                                           color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(ft_all_df, mapping = aes(
    ymin = C_pourc_mean - C_pourc_std.error,
    ymax = C_pourc_mean + C_pourc_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (C_pourc_mean + C_pourc_std.error + 2)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Total C (%)"), x = "Site") + # ~ "Lignin" 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
C_pourc
# ggplot2::ggsave('~/Desktop/LMM_C_pourc.png', C_pourc, width = 4.7, height = 2.4)

# C:N
C_N <- ggplot(C_N.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = C_N_ratio_mean, 
                                      color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(C_N.mult, mapping = aes(
    ymin = C_N_ratio_mean - C_N_ratio_std.error,
    ymax = C_N_ratio_mean + C_N_ratio_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (C_N_ratio_mean + C_N_ratio_std.error + 2)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = "C:N ratio", x = "Site") + # ~ "Nitrogen concentration"
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
C_N
#ggplot2::ggsave('~/Desktop/LMM_C_N.png', C_N, width = 4.7, height = 2.5)

# chlBaire
chlBaire <- ggplot(chlBaire.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = chlB_mg_m2_mean, 
                                                color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(chlBaire.mult, mapping = aes(
    ymin = chlB_mg_m2_mean - chlB_mg_m2_std.error,
    ymax = chlB_mg_m2_mean + chlB_mg_m2_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (chlB_mg_m2_mean + chlB_mg_m2_std.error + 5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Chl" ~ italic("b") ~ "(mg m" ^"-2"~")"),x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
# chlBaire
#ggplot2::ggsave('~/Desktop/LMM_chlBaire.png', chlBaire, width = 4.7, height = 2.5)


# carot
carot <- ggplot(carot.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = carotenoides_mg_g_mean, 
                                          color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(carot.mult, mapping = aes(
    ymin = carotenoides_mg_g_mean - carotenoides_mg_g_std.error,
    ymax = carotenoides_mg_g_mean + carotenoides_mg_g_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (carotenoides_mg_g_mean + carotenoides_mg_g_std.error + 0.5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Carotenoids (mg g" ^"-1"~")"), x = "Site") + # ~ "Total carotenoids" ~ 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
carot
#ggplot2::ggsave('~/Desktop/LMM_carot.png', carot, width = 4.7, height = 2.5)


# cell
cell <- ggplot(ft_all_df, mapping = aes(x = factor(site_class, levels = levels.site), y = cellulose_mean, 
                                        color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(ft_all_df, mapping = aes(
    ymin = cellulose_mean - cellulose_std.error,
    ymax = cellulose_mean + cellulose_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (cellulose_mean + cellulose_std.error + 5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Cellulose (%)"), x = "Site") +  
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
cell
ggplot2::ggsave('~/Desktop/LMM_cell.png', cell, width = 4.7, height = 2.4)


# Hcell
Hcell <- ggplot(ft_all_df, mapping = aes(x = factor(site_class, levels = levels.site), y = hemicellulose_mean, color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(ft_all_df, mapping = aes(
    ymin = hemicellulose_mean - hemicellulose_std.error,
    ymax = hemicellulose_mean + hemicellulose_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (hemicellulose_mean + hemicellulose_std.error + 5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Hemicellulose (%)"), x = "Site") + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
Hcell
ggplot2::ggsave('~/Desktop/LMM_Hcell.png', Hcell, width = 4.7, height = 2.4)


# solC
solC <- ggplot(ft_all_df, mapping = aes(x = factor(site_class, levels = levels.site), y = soluble_contents_mean, color = site_class)) +
  facet_grid(~ scientific_name) +
  geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
  geom_errorbar(ft_all_df, mapping = aes(
    ymin = soluble_contents_mean - soluble_contents_std.error,
    ymax = soluble_contents_mean + soluble_contents_std.error), width = 0) +
  # geom_text(mapping = aes(x = site_class, label = .group,
  #                         y = (soluble_contents_mean + soluble_contents_std.error + 5)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Site", values = c(15, 17, 16, 16, 16)) +
  scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
  labs(y = bquote("Soluble carbon (%)"), x = "Site") +  # ~ "Soluble carbons"
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"), 
        strip.text.x = element_text(face = "italic"))
# solC
# ggplot2::ggsave('~/Desktop/LMM_solC.png', solC, width = 4.7, height = 2.4)


# chlA_chlB
# chlA_chlB.lme = lme(chlA_chlB ~ scientific_name * site_class, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |site_class),  na.action = na.omit, method = "REML", data = ft_all)
# chlA_chlB.tuckey <- emmeans(chlA_chlB.lme, pairwise ~ scientific_name * site_class, adjust = "tukey")
# chlA_chlB.mult <- cld(chlA_chlB.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
# chlA_chlB.mult.prov <- unite(chlA_chlB.mult, col = "id", c(scientific_name, site_class))
# chlA_chlB.mult = full_join(chlA_chlB.mult.prov, ft_all_df[, c('id', "chlA_chlB_mean", "chlA_chlB_std.error")], by = "id")
# chlA_chlB.mult <-  separate(chlA_chlB.mult, id, c("scientific_name", "site_class"), sep = "_", remove = TRUE)
# 
# chlA_chlB <- ggplot(chlA_chlB.mult, mapping = aes(x = factor(site_class, levels = levels.site), y = chlA_chlB_mean,
#                                                   color = site_class)) +
#   facet_grid(~ scientific_name) +
#   geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = site_class)) +
#   geom_errorbar(chlA_chlB.mult, mapping = aes(
#     ymin = chlA_chlB_mean - chlA_chlB_std.error,
#     ymax = chlA_chlB_mean + chlA_chlB_std.error), width = 0) +
#   geom_text(mapping = aes(x = site_class, label = .group,
#                           y = (chlA_chlB_mean + chlA_chlB_std.error + 0.2)),
#             position = position_nudge(x = c(0,0,0,0,0))) +
#   scale_shape_manual(values = c(15, 17, 16, 16, 16)) +
#   scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#B10026", "#FC4E2A", "#FD8D3C")) +
#   #   labs(y = "Chlorophyll" ~  italic("a:b"), x = "Site") + # "Specific leaf area (chlA_chlB)" ~
#   theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.border = element_rect(0), panel.spacing = unit(0, "lines"),
#         strip.text.x = element_text(face = "italic"))
# # chlA_chlB
# #ggplot2::ggsave('~/Desktop/LMM_chlA_chlB.png', chlA_chlB, width = 4.7, height = 2.4)


# SANS effet site : pour le trait focntionnel suivant : lignine ####

# lignine
lign.lme = lme(lignine ~ scientific_name, random = ~1|plot_field_id, weights = varIdent(form = ~ 1 |scientific_name),  na.action = na.omit, method = "REML", data = ft_all)
lign.tuckey <- emmeans(lign.lme, pairwise ~ scientific_name, adjust = "tukey")
lign.mult <- cld(lign.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
lign.mult = full_join(lign.mult, ft_all_df_NO.eff_site)

lignine <- ggplot(lign.mult, aes(x = scientific_name, y = mean, color = scientific_name)) +
  geom_point(stat = "identity", position = position_dodge(width = .7), size = 3, aes(shape = scientific_name)) +
  geom_errorbar(lign.mult, mapping = aes(
    ymin = mean - std.error,
    ymax = mean + std.error), width = 0) +
  # geom_text(mapping = aes(x = scientific_name, label = .group,
  #                         y = (mean + std.error + 3)),
  #           position = position_nudge(x = c(0,0,0,0,0))) +
  scale_shape_manual(name = "Species", values = c(15, 17, 16, 18)) +
  scale_colour_manual(name = "Species", values = brewer.pal(5, "Dark2")[c(1:3, 5)]) +
  labs(y = bquote("Lignin (%)"), x = "Species") + # ~ "Lignin" 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"), 
        panel.border = element_rect(0), panel.spacing = unit(0, "lines"))
lignine
#ggplot2::ggsave('~/Desktop/LMM_lignine.png', lignine, width = 2.7, height = 2.95)


ggplot(EWT.mult.sc, mapping = aes(x = scientific_name, y = equivalent_water_thickness_cm_mean, 
                                  color = scientific_name)) +
  
  
  # AUTRES ####
# # phe
# phe.lme <- lme(Phenols_mg_g ~ scientific_name * site_class, random = ~1 | plot_field_id, weights = varIdent(form = ~ 1 |site_class), na.action = na.omit, data = ft_all) 
# phe.tuckey <- emmeans(phe.lme, pairwise ~ scientific_name * site_class, adjust = "tukey")
# phe.mult <- cld(phe.tuckey$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey")
# phe.mult.prov <- unite(phe.mult, col = "id", c(scientific_name, site_class))
# phe.mult = full_join(phe.mult.prov, ft_all_df[, c('id', "Phenols_mg_g_mean", "Phenols_mg_g_std.error")], by = "id")
# phe.mult <-  separate(phe.mult, id, c("scientific_name", "site_class"), sep = "_", remove = TRUE)
# 
# phe <- ggplot(ft_all_df, mapping = aes(y = mean, x = ft_all_df$site_class, color = site_class)) +
#   facet_grid(~ scientific_name) +
#   geom_point(stat = "identity", size = 3, position = position_dodge(width = .7), aes(shape = ft_all_df$site_class)) +
#   geom_errorbar(ft_all_df, mapping = aes(
#     ymin = mean - std.error,
#     ymax = mean + std.error), width = 0) +
#   # geom_text(phe.mult, mapping = aes(x = phe.mult$site_class, y = (phe.mult$upper.CL), label = phe.mult$.group)) +
#   scale_shape_manual(values = c( 15, 17, 16, 16, 16)) +
#   scale_color_manual(name = "Site", values = c(brewer.pal(9, "Greys")[c(7,7)], "#FD8D3C", "#FC4E2A", "#B10026")) +
#   labs(y = "Phenols (mg g" ^ -1 ~")", x = "Site") + # ~ "Leaf water content (phe)" 
#   theme_bw() +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.border = element_rect(0),
#         panel.spacing = unit(0, "lines"))
# phe
# #ggplot2::ggsave('~/Desktop/LMM_phe.png', phe, width = 4.7, height = 2.4)
