# 14 juillet

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                        Données Environnementales, Compilation
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : 

# VOIR : ARCHIVE_compilation+analyses_environnement.R dans doc "environnement"
# script initialement écrit pour le cours Legendre, upgradé et updaté ici


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(plyr)
library(ggplot2)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/plot_toutes_var.RData") # plot_toutes_var
# issu des manips il y a jadis dans ce script


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Manipulations pour obtenir "plot_toutes_var" ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# # Nitrates et ammonium et azote total
# chimie_tourbe_prov <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/environnement/chimie_tourbe.csv", sep = ";")
# chimie_tourbe <- chimie_tourbe_prov %>%
#   dplyr::select(-c(class_N_dep, date)) %>% 
#   dplyr::rename(nitrate_mg_kg_ou_L = "t_mg_NO3.kg", ammonium_mg_kg_ou_L = "t_mg_NH4.kg")
# chimie_tourbe$nitrate_mg_kg_ou_L <- ifelse(chimie_tourbe$nitrate_mg_kg_ou_L < 0, 0, chimie_tourbe$nitrate_mg_kg_ou_L) # retirer valeur négative
# 
# chimie_eau_prov <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/environnement/chimie_eau.csv", nrows = 18) # exclure 2 dernières lignes
# chimie_eau <- chimie_eau_prov %>%
#   dplyr::select(-c(classe_N_dep, date)) %>%
#   na.omit() %>%
#   dplyr::rename(nitrate_mg_kg_ou_L = "e_NO3.L", ammonium_mg_kg_ou_L = "e_mg_NH4.L", N_tot_mg_L = "e_N_tot_mg.N.L" )
# 
# azote_tot_tourbe_prov <- read.csv("~/Documents/Maîtrise/DonnéesAnalyses/environnement/CN_tourbe_Ntot.csv")
# azote_tot_tourbe_prov2 <-  azote_tot_tourbe_prov %>%
#   dplyr::select(site, plot_id, type_echantillon, N_pourc = "N_pourcent")
# 
# # ajouter 18 lignes vides pour pouvoir effectuer ligne 66 (if_else)
# ajout <- as.data.frame(matrix(data = NA, nrow = nrow(azote_tot_tourbe_prov2), ncol = ncol(azote_tot_tourbe_prov2)))
# names(ajout) <- names(azote_tot_tourbe_prov2)
# azote_tot_tourbe <- azote_tot_tourbe_prov2 %>% rbind(ajout)
# 
# # joindre pour créer "azote"
# azote <- full_join(chimie_tourbe, chimie_eau)
# 
# #azote$N_tot_mg_mg_ou_L[1:18] <- azote_tot_tourbe$N_tot_mg_mg_ou_L[1:18] # ajouter les valeurs N_tot_mg_mg_ou_L du df azote_tot_tourbe pour échantillons "tourbe"
# azote$N_pourc <- ifelse(azote$type_echantillon == "tourbe", azote_tot_tourbe$N_pourc, NA)
# 
# # pH et conductivité
# load("~/Documents/Maîtrise/DonnéesAnalyses/environnement/ph_cond.RData") # ph_cond
# # issu de "calcul_cond_corrigée.R" dans dossier environnement
# 
# # profondeur nappe phréatique
# prof_nappe <- na.omit(read.csv("~/Documents/Maîtrise/DonnéesAnalyses/environnement/prof_nappe.csv"))  %>% # message d'erreur pour le na.omit, mais ça fonctionne quand même
#   filter(! site %in% c("MBP_5N", "MBP_10N", "MBP_20N")) %>%
#   dplyr::select(-c("site", "latitude", "longitude", "type_echantillon", "date"))
# 
# # joindre les df
# plot_toutes_var <- full_join(azote, ph_cond) %>%
#   full_join(prof_nappe)
# 
# 
# # ajout nom commun pour plot_id
# load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/field_id.RData") # field_id
# # issu de manips sur données "ARCHIVEplot_toute_var2" (extraction de "no_plot" et "plot_id")
# 
# plot_toutes_var <- full_join(plot_toutes_var, field_id) %>% 
#   dplyr::select("site":"plot_id", "no_plot", "type_echantillon":"prof_nappe_cm")

# save(plot_toutes_var, file = "/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/plot_toutes_var.RData") # plot_toutes_var


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Obtenir moyenne et std par site ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

moy_std_site <- plot_toutes_var %>%
  group_by(site, type_echantillon) %>%
  dplyr::summarize_at(vars("nitrate_mg_kg_ou_L", "ammonium_mg_kg_ou_L", "N_tot_mg_L", "N_pourc", "pH", "corr_EC_microS_cm", "prof_nappe_cm"), funs(mean, sd)) %>%
  ungroup()
# write.csv(moy_std_site, "verif.csv")

# mettre le résultats STD dans la case MEAN entre parenthèse
{
eau_toutes_var <- plot_toutes_var %>% 
  dplyr::filter(type_echantillon == "eau_tourbe")
tourbe_toutes_var <- plot_toutes_var %>% 
  dplyr::filter(type_echantillon == "tourbe")

# sur eau
mean <- eau_toutes_var %>% 
  group_by(site) %>% 
  dplyr::summarize_at(vars("nitrate_mg_kg_ou_L", "ammonium_mg_kg_ou_L", "N_tot_mg_L", "pH", "corr_EC_microS_cm", "prof_nappe_cm"), funs(mean)) %>% 
  ungroup()

std <- eau_toutes_var %>% 
  group_by(site) %>% 
  dplyr::summarize_at(vars("nitrate_mg_kg_ou_L", "ammonium_mg_kg_ou_L", "N_tot_mg_L", "pH", "corr_EC_microS_cm", "prof_nappe_cm"), funs(sd)) %>% 
  ungroup()

# transformer en vecteurs
mvect <- round(c(t(mean[, 2:7])), 2)
stdvect <- round(c(t(std[, 2:7])), 2)

# coller les valerus std entre parenthèse aux valeurs moyenne
msd <- paste(mvect, " (", stdvect, ")", sep = "")

# mettre le vecteur résultant sous forme de data frame
rn <- seq(1:3)
cn <- colnames(eau_toutes_var[c(1, 6:9, 11:13)])

eau_prov <- array(dim = c(8, 3), dimnames = list(cn, rn))
eau_prov[3:8,] <- msd
colID<- matrix(data = c("GPB_canal", "GTV_SE", "MBP_open", rep("eau_tourbe", 3)), 3, 2)
eau_prov[c(1,2),] <- t(colID)
eau_prov2 <- data.frame(t(eau_prov)) 
eau <- plyr::rename(eau_prov2, c(nitrate_mg_kg_ou_L = "nitrate_mg_L", ammonium_mg_kg_ou_L = "ammonium_mg_L", 
                                 pH = "peat_wtr_pH", corr_EC_microS_cm = "peat_wtr_EC"))

# sur tourbe
mean2 <- tourbe_toutes_var %>% 
  group_by(site) %>% 
  dplyr::summarize_at(vars("nitrate_mg_kg_ou_L", "ammonium_mg_kg_ou_L", "N_pourc", "pH", "corr_EC_microS_cm", "prof_nappe_cm"), funs(mean)) 

std2 <- tourbe_toutes_var %>% 
  group_by(site) %>% 
  dplyr::summarize_at(vars("nitrate_mg_kg_ou_L", "ammonium_mg_kg_ou_L", "N_pourc", "pH", "corr_EC_microS_cm", "prof_nappe_cm"), funs(sd)) 

# transformer en vecteurs
mean_vect2 <- round(c(t(mean2[, 2:7])), 2)
std_vect2 <- round(c(t(std2[, 2:7])), 2)

# coller les valerus std entre parenthèse aux valeurs moyenne
msd2 <- paste(mean_vect2, " (", std_vect2, ")", sep = "")

# mettre le vecteur résultant sous forme de data frame
rn <- seq(1:3)
cn2 <- colnames(tourbe_toutes_var[c(1, 6:8, 10:13)])

tourbe_prov <- array(dim = c(8, 3), dimnames = list(cn2, rn))
tourbe_prov[ 3:8,] <- msd2
colID2<- matrix(data = c("GPB_canal", "GTV_SE", "MBP_open", rep("tourbe", 3)), 3, 2)
tourbe_prov[c(1,2),] <- t(colID2)
tourbe_prov2 <- data.frame(t(tourbe_prov))
tourbe <- plyr::rename(tourbe_prov2, c(nitrate_mg_kg_ou_L = "nitrate_mg_kg", ammonium_mg_kg_ou_L = "ammonium_mg_kg", 
                                 pH = "peat_pH", corr_EC_microS_cm = "peat_EC"))

# créer colonne "n"
n <- matrix(data = rep("6", 3))

# coller les matrices ensemble
msd_toutes <- cbind(eau$site, n, eau[, -c(1,2)], tourbe[, 3:7]) %>% 
  plyr::rename(c("eau$site" = "Site", n = "n", nitrate_mg_L = "NO3- (mg L-1)", nitrate_mg_kg = "NO3- (mg kg-1)", 
                 ammonium_mg_L = "NH4+ (mg L-1)", ammonium_mg_kg = "NH4+ (mg kg-1)", 
                 N_tot_mg_L = "N (mg L-1)", N_pourc = "N (%)",
                peat_wtr_pH = "pH (water)", peat_pH = "pH (peat)", 
                peat_wtr_EC = "corrected EC (µS/cm; water)", peat_EC = "corrected EC (µS/cm; peat)",
                prof_nappe_cm = "Water table depth (cm)")) %>% 
  dplyr::select( "Site", "n", "NO3- (mg L-1)", "NO3- (mg kg-1)", "NH4+ (mg L-1)", "NH4+ (mg kg-1)", "N (mg L-1)", "N (%)", 
                 "pH (water)", "pH (peat)" ,"corrected EC (µS/cm; water)", "corrected EC (µS/cm; peat)",
                 "Water table depth (cm)")
msd_toutes$Site = factor(msd_toutes$Site, levels = c('MBP_open','GTV_SE',"GPB_canal"), labels = c("MBP","GTV","GPB")) 
}

# renvoyer en .csv
# write.csv(msd_toutes, file = "/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/tableau_var_env.csv")



#################################### NE PAS FAIRE CETTE PARTIE #################################### 




#-=-=-=-=-=-=-=-=-=-=-=-=-
# Visualiation  ----
#-=-=-=-=-=-=-=-=-=-=-=-=-


# SCRIPT ALEXIS
# WTloss.plot.prov.horiz.mesh.CN.ecm <- ggplot(WT.mult.prov.horiz.CN.ecm, aes(x = horizon, y = emmean-1, fill = mesh)) +
#   facet_grid(~soil.type.prov, labeller = labeller(soil.type.prov = c('AS'='Maple soil provenance', 'FG'='Beech soil provenance'))) +
#   geom_bar(stat = "identity", colour="black", position=position_dodge(width = .7), width = .6) +
#   geom_errorbar(aes(ymin = emmean-1-SE, ymax = emmean-1+SE), position=position_dodge(width=0.7), width = .2, size = .5) +
#   #geom_text(aes(y = emmean-1+SE, label = .group, hjust = .5, vjust = -.4), position=position_dodge(width=0.7), show.legend = FALSE) +
#   scale_fill_manual(name="Mycorrhizal\nhyphae\nexcluded", labels=c("Yes", "No"), values = c("red", "white")) +
#   labs(x="Horizons", y= "C:N ratio change in EcM forest") +
#   scale_y_continuous(labels = function(y) y + 1)
# # Color strip
# g <- ggplot_gtable(ggplot_build(WTloss.plot.prov.horiz.mesh.CN.ecm))
# strip_both <- which(grepl('strip-', g$layout$name))
# fills <- c("#b2df8a","#33a02c")
# k <- 1
# for (i in strip_both) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# library(grid)
# grid.draw(g, recording=TRUE)
# png("CN/CNchange.prov.horiz.mesh.se.png", width = 6, height = 4, units = "in", res = 300); print(grid.draw(g)); dev.off()



#
###
#####
#######
#########  RENDUE LÀ
#######
#####
###
#






# format long
boxplot <- gather(data = plot_toutes_var, key = "variable", value = "value", c(nitrate_mg_kg_ou_L:prof_nappe_cm))
boxplot$variable <- as.factor(boxplot$variable)
boxplot$variable = droplevels(boxplot)$variable

# séparer les types d'échantillon
# eau
boxplot_eau <- boxplot %>% 
  dplyr::filter(type_echantillon == "eau_tourbe")
boxplot_eau$variable = droplevels(boxplot_eau)$variable

# pour que ggboxplot ne les mette pas en ordre alphabetique
boxplot_eau$variable_ordered =  factor(boxplot_eau$variable, levels = c("ammonium_mg_kg_ou_L", "nitrate_mg_kg_ou_L", "N_tot_mg_L", "pH", "corr_EC_microS_cm", "prof_nappe_cm"))

# renommer les variables
var_lab_eau <- c(ammonium_mg_kg_ou_L = expression(paste("Ammonium (mg L" ^ "-1", ")")), nitrate_mg_kg_ou_L = 'Nitrate (mg/L)', N_tot_mg_L = 'Total nitrogen (mg/l)', pH = 'pH',
             corr_EC_microS_cm = 'Conductivity (µS/cm)', prof_nappe_cm = 'Water table depth (cm)' )

# tourbe
boxplot_tourbe <- boxplot %>% 
  dplyr::filter(type_echantillon == "tourbe")






# CA MARCHE PAS !
boxplot_tourbe$variable = droplevels(boxplot_tourbe$variable)







# pour que ggboxplot ne les mette pas en ordre alphabetique
boxplot_tourbe$variable_ordered =  factor(boxplot_tourbe$variable, levels = c("ammonium_mg_kg_ou_L", "nitrate_mg_kg_ou_L", "N_pourc", "pH", "corr_EC_microS_cm", "prof_nappe_cm"))

# renommer les variables
var_lab_tourbe <- c("nitrate_mg_kg_ou_L" = 'Nitrate (mg/kg)', "ammonium_mg_kg_ou_L" = expression(paste("Ammonium (mg kg" ^ "-1", ")")),
                    "N_pourc"= 'Total nitrogen (%)', "pH" = 'pH', "corr_EC_microS_cm" = 'Conductivity (µS/cm)', "prof_nappe_cm" = 'Water table depth (cm)')

# gosh = c(expression(paste("Ammonium (mg kg" ^ "-1", ")")),
#   expression(paste("Nitrate (mg kg" ^ "-1", ")")),
#   expression(paste("Total N (mg mg" ^ "-1", ")")),
#   'Total nitrogen (%)', 'pH',
#   expression(paste("Conductivity (µS cm" ^ "-1", ")")),
#   'Water table depth (cm)')


# visualisation
boxplot_env <- ggplot(boxplot_tourbe, aes(x = site, y = value, fill = site)) + 
  geom_boxplot() +
  facet_wrap(~ variable_ordered, scales = "free",  labeller = labeller(variable = var_lab_tourbe)) #labeller(variable = c(expression(paste("Ammonium (mg kg" ^ "-1", ")")),
  #                                                                             expression(paste("Nitrate (mg kg" ^ "-1", ")")),
  #                                                                             expression(paste("Total N (mg mg" ^ "-1", ")")),
  #                                                                             'Total nitrogen (%)', 'pH',
  #                                                                             expression(paste("Conductivity (µS cm" ^ "-1", ")")),
  #                                                                             'Water table depth (cm)')  )) # label_parsed) # labeller(variable_ordered = gosh))
  labs(y = '', x = '') +
  scale_fill_manual(name = 'Site',
                    values = c("tomato3", "sienna2","goldenrod2"),
                    labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  #labs(title = c(ammonium_mg_kg_ou_L = expression(paste("Ammonium (mg " ^ "-1", ")")))) +
  scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggtitle("Données environnementales sur échantillons de tourbe") # et eau de tourbe")
boxplot_env
#ggsave('boxplot_env.png', width = 7, height = 7)
#ggsave('boxplot_env.pdf', width = 10, height = 10)









# graph en français
var_lab_FR <- c(t_mg_NH4.kg = 'Ammonium (mg/kg)', t_mg_NO3.kg = 'Nitrate (mg/kg)', t_N_pourcent = 'Azote total (%)', t_pH = 'pH',
                t_cond_corr_mV = 'Conductivité (mV)', prof_nappe_cm = 'Profondeur de la nappe phréatique (cm)' ) 

boxplot_env_francais <-ggplot(boxplot, aes(x = site, y = valeur, fill = site)) + 
  geom_boxplot() +
  facet_wrap(~ variable_ordered, scales = "free",  labeller = as_labeller(var_lab_FR)) +
  labs(y = '', x = '') +
  scale_fill_manual(name = 'Site',
                    values = c("tomato3", "sienna2","goldenrod2"),
                    labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
boxplot_env_francais
#ggsave('boxplot_env_fr.png', width = 7, height = 7)
#ggsave('boxplot_env_fr.pdf', width = 10, height = 10)



#=======================================# ARCHIVES #=======================================#                                                       ----


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Analyses Legendre (ARCHIVES) ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# objet valeurs moyennes standardisées par site, 3 x 6
# Get enviromental matrix 2 with same number of rows as spectra
moy_st_site <- ddply(plot, .(site), numcolwise(mean)) %>% # si je veux de nouvelles colonnes, mettre plot à place de var_imp
  dplyr::select(site, t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, t_pH, t_cond_corr_mV, prof_nappe_cm) %>% # changer noms des colonnes ici
  remove_rownames %>% 
  column_to_rownames(var="site")  %>% # 1iere colonne = rownames
  decostand(method = "standardize") %>% # matrice centrée-réduite (et/ou scale = T dans ordinations) 
  cbind(c("GPB_canal", "GTV_SE", "MBP_open"))
#write.csv(moy_st_site, file = "tableau_var_env.csv")
colnames(moy_st_site)[7] <- "site"

# objet moyenne par site + sd : DEVOIR LEGENDRE
moy.sd_site <- plot %>% 
  dplyr::select(site, t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, 
                t_pH, t_cond_corr_mV, prof_nappe_cm) %>% 
  group_by(site) %>% 
  dplyr::summarize_all(funs(mean, sd))
#write.csv(moy.sd_site, file ="tableau_var_env.csv")

# Get environmental matrix X2 with same number of rows as leaf spectra matrix (Y2)
#load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/spectra.adonis.2.RData") # spectra.adonis
env_mat2 <- spectra.adonis.2 %>% 
  left_join(moy_st_site, by = "site") %>% 
  dplyr::select(t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, t_pH, t_cond_corr_mV, prof_nappe_cm) 
#save(env_mat2, file = "env_mat2.RData")

# Get environmental matrix X3 with same number of rows as leaf spectra matrix (Y3)
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/environnement/spectra.adonis.6.RData") # spectra.adonis

env_mat3 <- spectra.adonis.6 %>% 
  left_join(Nlevels, by = "site") %>% 
  dplyr::select(site, N_fert_level)
#save(env_mat3, file = "env_mat3.RData")


# ---
# DENDROGRAMMES 
# toutes variables environnementales ####
# dist_site <- vegdist(moy_st_site, "euclidean") # aller voir la théorie !!
# gr_site <- hclust(dist_site, method = "ward.D2") # analogue à ANOVA, tente de réduire variance intra
# 
# fig_gr_site <- ggdendrogram(gr_site) + 
#   labs(x = "bla", 
#        y = "distance",
#        title = "Groupement des sites selon variables environnementales") # les noms d'axes apparaissent pas ! je gosse pas a dessus lala 
# fig_gr_site
# # ggsave("fig_gr_site.png", width = 3, height = 3)
# 
# # interprétation résultats dendrogramm final
# # Le trois sites sont environs aussi loin les 2 des autres
# # utilisation de toutes les var environnementales (pH, cond, NH4, NO3, %N dans tourbe)

# # DENDROGRAMME variables AZOTE - WARD ####
# moy_st_site_N <- moy_st_site %>% 
#   select(t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent) # sans pH, cond
# dist_N <- vegdist(moy_st_site_N, "euclidean") # aller voir la théorie !!
# 
# gr_N <- hclust(dist_N, method = "ward.D2")
# fig_gr_N <- ggdendrogram(gr_N) + 
#   labs(x = "bla", 
#        y = "distance",
#        title = "Sites selon variables en lien avec N") # les noms d'axes apparaissent pas ! je gosse pas a dessus lala 
# fig_gr_N
# # NH4 a GPB est très élevé ! cela correspond certainement a des poussières de fertilisant utilisé dans l'agriculture alentour..
# # groupes se distinguent, GTV et MBP ensemble
# 
# # DENDROGRAMME variables AZOTE - UPGMC ####
# moy_st_site_N <- moy_st_site %>% 
#   select(t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent) # sans pH, cond
# dist_N <- vegdist(moy_st_site_N, "euclidean") # aller voir la théorie !!
# 
# gr_N <- hclust(dist_N, method = "centroid")
# fig_gr_N <- ggdendrogram(gr_N) + 
#   labs(x = "bla", 
#        y = "distance",
#        title = "Sites selon variables en lien avec N") # les noms d'axes apparaissent pas ! je gosse pas a dessus lala 
# fig_gr_N
# # NH4 a GPB est très élevé ! cela correspond certainement a des poussières de fertilisant utilisé dans l'agriculture alentour..
# # groupes se distinguent, GTV et MBP ensemble


# # ---
# # Dendrogramme ####
# # WARD #
# dist_plots_N <- val_site %>% 
#   dplyr::select(no_plot, t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, t_pH, t_cond_corr_mV) %>%
#   remove_rownames %>% 
#   column_to_rownames(var = "no_plot")  %>% # 1iere colonne = rownames
#   decostand(method = "standardize") %>% 
#   vegdist("euclidean", diag = T, upper = T) 
# 
# WARD_plots_N <- hclust(dist_plots_N, method = "ward.D2") # ward.D2 BON voir cours 6#97
# 
# (WARD <- ggdendrogram(WARD_plots_N))
# 
# # montre grande hétérogénéité !
# 
# # WPGMC #
# dist_plots_N <- val_site %>% 
#   dplyr::select(no_plot, t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, t_pH, t_cond_corr_mV) %>%
#   remove_rownames %>% 
#   column_to_rownames(var = "no_plot")  %>% # 1iere colonne = rownames
#   decostand(method = "standardize") %>% 
#   vegdist("euclidean", diag = T, upper = T) 
# 
# WPGMC_plots_N <- hclust(dist_plots_N, method = "median") 
# 
# (WPGMC <- ggdendrogram(UPGMC_plots_N, labels = T))




# Minimum spanning tree azote ####
span.dist <- spantree(dist_plots_N)
plot(span.dist, type = "t")


# ---
# Analyse composantes principales ####
# ACP ####

# enlever colonne site
val_site_rda <- val_site %>% 
  dplyr::select(t_mg_NO3.kg, t_mg_NH4.kg, t_N_pourcent, t_pH, t_cond_corr_mV) %>% 
  decostand(moy_site, method = "standardize") # matrice centrée-réduite (et/ou scale = T dans ordinations)

# faire pca
pca_site <- rda(val_site_rda) # une matrice, fait donc une PCA automatiquement, scale = T -> corrélation (centrage-réd)

# visualiser la rda 
dtp <- data.frame('site' = val_site$site, scores(pca_site, display = 'site')[1:18,]) 

pca.plot <- ggplot(dtp) +
  geom_point(aes(x = PC1, y = PC2, color = dtp$site)) +
  geom_segment(data = df2, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = "grey22", arrow = arrow(length = unit(0.01,"npc"))) +
  geom_text(data = df2, 
            aes(x = PC1, y = PC2, label = rownames(df2),
                hjust = 0.5*(1 - sign(PC1)), vjust = 0.5*(1 - sign(PC2))), 
            color = "grey22", size = 3) +
  xlab('PCA axis 1 (52.3 %)') + 
  ylab('PCA axis 2 (28.1 %)') +
  labs(title = 'PCA, cadrage de type 1(?) des sites échantillonnés et variables environnementales') +
  scale_colour_manual(name = 'Site', 
                      values = c("goldenrod2", "sienna2", "tomato3"),
                      labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  stat_ellipse(aes(x = PC1, y = PC2, color = dtp$site)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_fixed() +
  theme_minimal() 
pca.plot
#ggsave('pca.sites.png', width = 7, height = 7)
#ggsave('pca.sites.pdf', width = 7, height = 7)


# ---
# Kaiser-Guttmen et le modèle du baton brisé ####
# Quelles axes expliquent plus que le hasard ? (voir théorie BIO6077/Manuel de cours/page 133)

# sortir eigen_values du modèle d'ordination retenu
eigen_values <- pca_site$CA$eig

# appliquer critère Kaiser-Guttman
eigen_values[ eigen_values > mean(eigen_values)]

# broken stick model
n <- length(eigen_values)
br_st_mod <- data.frame(j = seq(1:n), p = 0)
br_st_mod$p[1] <- 1/n
for (i in 2:n) {
  br_st_mod$p[i] = br_st_mod$p[i-1] + (1/n + 1 - i)
}
br_st_mod$p <- 100*br_st_mod$p/n
br_st_mod

# afficher le tout dans un graph
bla <- ggplot(eigen_values) + 
  labs(x = "", y = "",
       main = "% variance") +
  geom_bar() +
  bla

# selon manuel :
par(mfrow = c(1,2))
barplot(eigen_values, main = "Eigenvalues", col = "bisque", las = 2)
abline(h = mean(eigen_values, col = "red")) # average eigenvalues !
legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
barplot(t(cbind(100*eigen_values/sum(eigen_values), br_st_mod$p[n:1])), beside = T, 
        main = "% variance", col = c("bisque", 2), las = 2)
legend("topright", c("% eigenvalue", "Broken stick model"), pch = 15, col = c("bisque", 2), bty = "n")

# pourquoi pas de rouge pr pc1 ?
# interpréter
# mettre cute sur ggplot
# mettre belles couleur OMG ARKCALISS

# tourbe et eau de tourbe n'ont pas le mm impact sur les sites !
# sites très différents *


# ---
#### Covariance des variables ####

# 3.1 
# utiliser tous les plots (réplicats) vs la moyenne... mais pas considérer les points individuels !

# Azote covariables ?
ggpairs(plot,
        cardinality_threshold = 20,
        columns = c("t_mg_NO3.kg","t_mg_NH4.kg", "e_NO3.L", "e_mg_NH4.L", "t_N_pourcent", "e_N_tot_mg.N.L"),
        columnLabels = c("Nitrate tourbe", "Ammonium tourbe", "Nitrate eau", "Ammonium eau","Azote tourbe",
                         "Azote eau"),
        mapping = aes(color = site),
        upper = list(continuous = wrap('cor', method = "spearman"))) # valeurs de corrélation paire par pair 
# plus d'ammonium que de nitrates, par conséquent logique si plus de corrélation entr Ntot et ammonium selon eau (0.585) ou tourbe (0.549)
# Nitrates (tourbe et eau) aussi bcp corrélé avec Nitrate eau et azote total (+ mobile ?) et varie de la mm façon mm si contribution moindre
# nitrate tourbe PAS corrélé avec azote tourbe (trop faible qté ?)
# voir cet article : Bridgham, S. D., Updegraff, K., & Pastor, J. (1998). Carbon, nitrogen, and phosphorus mineralization in northern wetlands. Ecology, 79(5), 1545-1561.


# Autres covariables ?
ggpairs(plot,
        columns = c("e_pH", "e_cond_corr_mV", "e_sal_psu", "e_dO_pourcent", "t_pH", "t_cond_corr_mV", "t_sal_psu", "t_dO_pourcent",  "prof_nappe_cm"),
        columnLabels = c("pH eau", "Conductivité eau", "Salinité eau","O2 dissout eau", "pH tourbe", "Conductivité tourbe", 
                         "Salinité tourbe","O2 dissout tourbe", "Profonddeur nappe"),
        mapping = aes(color = site),
        upper = list(continuous = wrap('cor', method = "spearman"))) # valeurs de corrélation paire par pair 
# pH de tourbe et eau (comme attendu) corrélé avec conductivité respective (moy = 0.6655)
# salinité inversement corrélée à conductivité respectivement pour mm type (-0.6..)
# je voulais voir si certaine variables étaient redondantes
# profondeur de nappe corrélée (>0.4 + ou -) à presque toutes les variables, hypopthèse : plus profond plus grand apport de
# matières plus ou moins décomposées ?

# 3.2
# moyennes par site
# Azote covariables ?
ggpairs(val_site,
        columns = c("t_mg_NO3.kg","t_mg_NH4.kg", "e_NO3.L", "e_mg_NH4.L", "t_N_pourcent", "e_N_tot_mg.N.L"),
        columnLabels = c("Nitrate tourbe", "Ammonium tourbe", "Nitrate eau", "Ammonium eau","Azote tourbe",
                         "Azote eau"),
        mapping = aes(color = rownames(val_site)),
        upper = list(continuous = wrap('cor', method = "spearman"))) # valeurs de corrélation paire par pair 
# 

# Autres covariables ?
ggpairs(val_site,
        columns = c("e_pH", "e_cond_corr_mV", "e_sal_psu", "e_dO_pourcent", "t_pH", "t_cond_corr_mV", "t_sal_psu", "t_dO_pourcent",  "prof_nappe_cm"),
        columnLabels = c("pH eau", "Conductivité eau", "Salinité eau","O2 dissout eau", "pH tourbe", "Conductivité tourbe", 
                         "Salinité tourbe","O2 dissout tourbe", "Profonddeur nappe"),
        mapping = aes(color = rownames(val_site)),
        upper = list(continuous = wrap('cor', method = "spearman")))

# on a des corrélations de +ou- 1 et 0.5, pas un très bonne idée du phénomène...


# ---
# PERMANOVA ####


# ---
# BOXPLOT ####
# format long
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/RDA/plot_toute_var.RData") # plot
boxplot <- gather(data = plot, key = "variable", value = "valeur", t_mg_NO3.kg:prof_nappe_cm) %>% 
  dplyr::filter(variable %in% c("t_mg_NO3.kg", "t_mg_NH4.kg","t_N_pourcent", # sélectionner intéressants et bon ordre
                                "t_pH","t_cond_corr_mV","prof_nappe_cm"))
boxplot$variable <- as.factor(boxplot$variable)
boxplot$variable = droplevels(boxplot)$variable
boxplot$variable_ordered = factor(boxplot$variable, levels=c("t_mg_NH4.kg", "t_mg_NO3.kg", "t_N_pourcent", "t_pH", "t_cond_corr_mV", "prof_nappe_cm"))

# renommer les variables
var_lab <- c(t_mg_NH4.kg = 'Ammonium (mg/kg)', t_mg_NO3.kg = 'Nitrate (mg/kg)', t_N_pourcent = 'Total nitrogen (%)', t_pH = 'pH',
             t_cond_corr_mV = 'Conductivity (mV)', prof_nappe_cm = 'Water table depth (cm)' ) 

boxplot_env <-ggplot(boxplot, aes(x = site, y = valeur, fill = site)) + 
  geom_boxplot() +
  facet_wrap(~ variable_ordered, scales = "free",  labeller = as_labeller(var_lab)) +
  labs(y = '', x = '') +
  scale_fill_manual(name = 'Site',
                    values = c("tomato3", "sienna2","goldenrod2"),
                    labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggtitle("Données environnementales (pas toutes incluses) sur échantillons de tourbe")
boxplot_env
#ggsave('boxplot_env.png', width = 7, height = 7)
#ggsave('boxplot_env.pdf', width = 10, height = 10)

# graph en français
var_lab_FR <- c(t_mg_NH4.kg = 'Ammonium (mg/kg)', t_mg_NO3.kg = 'Nitrate (mg/kg)', t_N_pourcent = 'Azote total (%)', t_pH = 'pH',
                t_cond_corr_mV = 'Conductivité (mV)', prof_nappe_cm = 'Profondeur de la nappe phréatique (cm)' ) 

boxplot_env_francais <-ggplot(boxplot, aes(x = site, y = valeur, fill = site)) + 
  geom_boxplot() +
  facet_wrap(~ variable_ordered, scales = "free",  labeller = as_labeller(var_lab_FR)) +
  labs(y = '', x = '') +
  scale_fill_manual(name = 'Site',
                    values = c("tomato3", "sienna2","goldenrod2"),
                    labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
boxplot_env_francais
#ggsave('boxplot_env_fr.png', width = 7, height = 7)
#ggsave('boxplot_env_fr.pdf', width = 10, height = 10)


# ---
# VIOLIN PLOTS ####
env_violin_errorbar <-ggplot(boxplot, aes(x = site, y = valeur, fill = site)) + 
  #geom_violin() +
  #geom_errorbar(aes(ymin = ,                  * comment on calcule ça ?
  #                  ymax = )) +
  facet_wrap(~ variable_ordered, scales = "free",  labeller = as_labeller(var_lab)) +
  labs(y = '', x = '') +
  scale_fill_manual(name = 'Site',
                    values = c("tomato3", "sienna2","goldenrod2"),
                    labels =c("Plée Bleue, Lévis", "Grande Tourbière, Villeroy", "Mer Bleue, Ontario")) +
  scale_x_discrete(breaks = NULL) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggtitle("Données environnementales (pas toutes incluses) sur échantillons de tourbe")
env_violin_errorbar

# ---
# THÉORIE POUR COMPRENDRE DONNÉES ####

# A) pourquoi plus d'ammonium (NH4) dans mon sol ? (paired covariance plot, valeur/plot/site)
# ammonium -> nitrates (nitrification, fait par bactéries, aérobie), pas (bcp) de nitrification si waterlogged soil!
# aussi si pH entre 5.5 et 7...
# donc normal si plus d'amonium !! / ammonium (NH4) normalement drectement dispo aux plantes
# (ammonium augmente l'acidification du sol)
# Nitrates directement disponible aux plantes !

# littérature mentionne aussi disponibilité en P et water table level comment affecte la dispo en N (Limpens, J., Berendse, F., & Klees, H. (2003). 
# N deposition affects N availability in [...] )
# donc si augm dispo de P, augm de dispo de N !
# eux : Dry and wet atmospheric deposition of nitrogen, phosphorus and silicon in an agricultural region
# mentionnent que met deposition of NHx important -> fort income d'ammonium par la pluie dans 
# des condition d'agriculture !

# ---
# Test de normalité ####
shapiro.test(var_importantes$t_mg_NO3.kg) # normal
shapiro.test(var_importantes$t_mg_NH4.kg) # pas normal
shapiro.test(var_importantes$t_N_pourcent) # normal
shapiro.test(var_importantes$t_pH) # pas normal
shapiro.test(var_importantes$t_cond_corr_mV) # pas normal
shapiro.test(var_importantes$prof_nappe_cm) # normal


# ---
# ARCHIVES ####

# obtenir une matrice avec les valeurs moyennes
toto <- var_importantes %>% 
  group_by(site) %>% 
  summarise(t_mg_NO3.kg = mean(t_mg_NO3.kg), 
            t_mg_NH4.kg = mean(t_mg_NH4.kg), 
            t_N_pourcent = mean(t_N_pourcent), 
            t_pH = mean(t_pH),
            t_cond_corr_mV = mean(t_cond_corr_mV),
            prof_nappe_cm = mean(prof_nappe_cm))
site <- toto %>% remove_rownames %>% column_to_rownames(var="site") # pour que la standardisation comprenne que 1iere colonne=rownames
