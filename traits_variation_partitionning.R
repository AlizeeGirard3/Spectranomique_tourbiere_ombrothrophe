# 15 oct

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
#                      Partitionnement de la variation - TRAITS                         #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# NOTES : 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(ggplot2)
library(tidyr)
library(spectrolab)
library(tibble)
library(vegan)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-
setwd("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels")

load("tous_analyses_vX2.RData")
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


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Manipulations préalables  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

## définir les variables explicatives en variables binaires, enlever col.redondantes
sc.prov <- model.matrix(~ ft_all$sc)  #as.matrix()
scientific_name <- sc.prov[, 2:4]

# st.clss.prov <- as.vector.factor(ft_all$)
st.clss.prov <- as.matrix(model.matrix(~ ft_all$site_class))
site_class <- st.clss.prov[, 2:5]

# distance entre les samples -> TF
## centrer-réduire + distance de Hellinger
ft <- vegan::decostand(ft_all[, 5:22], method = "standardize", na.rm = T) %>% 
  vegan::vegdist("euclidean", diag = T, upper = T) 


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Partitionnement de la variance -- TRAITS
#-=-=-=-=-=-=-=-=-=-=-=-=-

var.part.traits <- varpart(ft, site_class, scientific_name) # fait le 9 avril 2020
# 
# Partition of squared  distance in dbRDA 
# 
# Call: varpart(Y = ft, X = site_class, scientific_name)
# 
# Explanatory tables:
# X1:  site_class
# X2:  scientific_name 
# 
# No. of explanatory tables: 2 
# Total variation (SS): 1660.1 
# No. of observations: 94 
# 
# Partition table:
#                       Df  R.squared  Adj.R.squared Testable
# [a+b] = X1            4   0.04570       0.00281     TRUE
# [b+c] = X2            3   0.69920       0.68918     TRUE
# [a+b+c] = X1+X2       7   0.75131       0.73107     TRUE
# Individual fractions                                    
# [a] = X1|X2           4                 0.04189     TRUE
# [b]                   0                -0.03908    FALSE
# [c] = X2|X1           3                 0.72826     TRUE
# [d] = Residuals                         0.26893    FALSE
# ---
#   Use function ‘dbrda’ to test significance of fractions of interest
plot(var.part.traits, bg = 1:4, digits = 2, cutoff = 0.0001) # $adj.r.squared = 0.7326357


## vérifier résultats avec RDA (devrait donner la même chose)
rda.traits <- capscale(ft ~ site_class + scientific_name) # Inertia is variance, proportion is r.squared (in vegan)
RsquareAdj(rda.traits) # un peu différent... pas d'intération dans ma formule ?
anova(rda.traits)

## Use function ‘dbrda’ to test significance of fractions of interest. Semi-Partial correlation
cpsc.site <- dbrda(ft ~ site_class) # condition pas disponible dans capscale !
RsquareAdj(cpsc.site) #  pas significatif # $adj.r.squared = 0.002807663
anova(cpsc.site) #  P-val = 0.388

cpsc.sc <- dbrda(ft ~ scientific_name)
RsquareAdj(cpsc.sc) # $adj.r.squared = 0.6891779
anova(cpsc.sc) # P-val <0.001 ***

