# 2 JUILLET 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                               PLS-DA TOUTES ESPÈCES                                  #  
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(plyr) # revalue()
library(dplyr)
library(spectrolab)
library(caret) # Classification And REgression Training
library(plotly) # plots
library (agricolae) # Tuckey test
library(reshape)
library(corrplot) # plot de corrélation
library(pls) # importance des bandes
library(ggplot2)
library(tidyverse) # rownames_to_column()


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# importer données spectrales
load("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_normalized_large.RData") # rt_sg_normalized_large
refl_norm <- rt_sg_normalized_large %>% 
  dplyr::filter(propriete == "reflectance") # issu du script "spectre_normalisation.R" dans document scripts-NETTOYAGE

# importer modèles directement
# --->>> modèle test n_comps <<<---
# mods_tous <- readRDS("~/DA/savedRDS/ttsp_50iter_normalized_PLSDA_mods.rds") # enregistré le ? / refait le 5 août avec vector-normalized

# --->>> modèle final <<<---
# fin_mods_tous <- readRDS("~/DA/savedRDS/ttsp_50iter_normalized_PLSDA_fin_mods.rds") # enregistré le ? / refait le 6 août avec vector-normalized

# --->>> confus final <<<---
# confus <- readRDS("~/DA/savedRDS/confus_final.rds") 


#=======================================# PLS ITERATIVE : CHOIX DE NCOMP  #=======================================#                                                      ----

#-=-=-=-=-=-=-=-=-=-=-=-=- 
# Coder la boucle sans préciser de ncomp                                                                                                               ----
#-=-=-=-=-=-=-=-=-=-=-=-=-     

# faire une table
dada <- as.data.frame(refl_norm)
table(dada$scientific_name)

# créer objet de classes qu'on veut prédire
classi <- dada$scientific_name # facteurs *

n_iter <- 50 # beacoup d'itérations, entre 50 et 100

# déterminer quels échantillons dans train et test [... plusieurs lignes]
rndm_id <- list() # specify no samples per species for training, alternative to % partitioning
set.seed(1840)
for (i in seq(n_iter)){ # fonction pour inventer une séquence de 1:50 / 50 fois
  rndm_id[[i]] <- with(dada, ave(1:nrow(dada), scientific_name, FUN = function(x) {sample.int(length(x))}))
}

n_iter <- 50
mods_tous <- list() 

# code de l'itération
for (i in seq(n_iter)){               
  inTrain <- createDataPartition(y = classi, p = .70, list = FALSE) # 70 % des espèces dans le groupe training
  print(i)
  flush.console() 
  set.seed(i)
  traini <- dada[inTrain, 6:2005] 
  testi <- dada[!(inTrain), 6:2005]
  trainclass <- classi[inTrain]; 
  testclass <- classi[!(inTrain)]
  plsFit <- train(traini, trainclass, method = "simpls", tuneLength = 10,
                  probMethod = "Bayes", trControl = trainControl(method = "boot")) 
  mods_tous[[i]] <- plsFit
}

# EXPLICATIONS ----
# la formule ci-haut demande 
# √ pour chaque (i) dans chaque fois qu'on va rouler la fxn (50 fois)
# √ inTrain <- createDataPartition(y = classi, p = , list = FALSE) (quick) / dit TRUE si valeur <= au chiffre mentionné
# √ imprime i (après l'itération)
# √ et invente des groupes de training et de test selon les valeurs de inTrain obtenues
# √ ensuite fait une liste (mods) de i et fait une plsFit dessus
# ----

# saveRDS(mods_tous, "savedRDS/ttsp_50iter_normalized_PLSDA_mods.rds") # long

head(mods_tous[[2]]$results)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Choix de ncomp                                                                          ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

n_iter = 50

# créer n_comps à partir des résultats du modèle mods_tous
n_comps <- vector(length = n_iter)
for (i in  1:n_iter){
  n_comps[i] <- mods_tous[[i]]$finalModel$ncomp 
}

# proposé ?
table(n_comps) 
# n_comps           PAS NORMALISÉ, 50 ITER
# 3  4  5  6  7  8 
# 1 11 20  9  6  3 
# n_comps             NORMALISÉ, 50 ITER, pas les mm résultats (mais mm conclusion)
# 3  4  5  6  7  8 
# 2 13 14  9 10  2  # 14x proposé de prendre 5 ncomps, 13x 14 ncomps

########################  selon un test de Tukey sur résultats kappa  ######################## 

# créer une table de données
kappas <- data.frame(n_comps = 1:10, matrix(NA, nrow = 10, ncol = length(mods_tous))) # n_comps vont de  1 à 10
for (i in 1:length(mods_tous)){  
  kappas[, i+1] <- mods_tous[[i]]$results$Kappa # rouler mods_tous[[i]]$results$Kappa = 10, nrow = 10
}

# créer une table de données *de la transpose* des résultats kappa
kapp <- as.data.frame(as.numeric(t(kappas[,-1]))) %>% # t() retourne la transpose (500 lignes, 50 modèles * 10 ncomp)
  cbind(rep(1:10, each = length(mods_tous))) # colle à la transpose des résultats kappa, le chiffre 1-10 respectivement pour chaque ncomp
colnames(kapp) <- c("Kappa", "n_comps") # noms de colonne

# transformer en facteurs
kapp$n_comps <- as.factor(kapp$n_comps)

# fitter un modèle linéaire sur lequel on fera test de Tuckey
modi <- lm(kapp$Kappa ~ n_comps, kapp) # on veut kappa en fxn du nombre de composantes

# test Tuckey
tuk <- HSD.test(modi, trt = "n_comps") # mutliple comparisons of treatments with Tuckey

tuk_data <- as.data.frame(tuk$groups) %>% # df
  mutate(var = as.numeric(row.names(tuk$groups)))  # ajouter colonne = nom de ligne (groupé selon groupe alors pas en ordre)

tuk_data <- tuk_data[order(tuk_data$var,decreasing = F),] 

# créer objet "letters" pour visualisation des groupes
letters <- as.character(tuk_data$groups)

# visualisation des résultats "Kappa"
ggplot(kapp, aes(x = n_comps, y = Kappa)) +
  geom_boxplot() +
  annotate(geom = "text", x = 1:10, y = rep(1, 10), label = letters) +
  theme_bw() +
  labs(x = "Number of components", y = "Kappa", 
       title = "test de Tukey sur résultats kappa pour 
       choix de n_comps; paramètres = 
       refl_sg_large_famille_sp_obj_normalized / NORMALISÉ / discr. espèces/ fait 5/08/19,
       70-30%, 50 iter, enregistré 5 août")
ggsave('tuckey_2août_PLSDA_normalized_ttsp.pdf', height = 10, width = 12)
# 4 ou 5 ncomps sont bons (4,5,6 sont dans le même groupe, premiers plus haut = 4 )


#=======================================# MODEL FINAL #=======================================#                    ----

# nombre de composants
n_comps <- 4 # nombre de composants choisi
fin_mods_tous <- list()
n_iter = 50 # 50 itérations 

for (i in 1:n_iter){ # same as 1:50 
  print(i)
  flush.console()
  set.seed(i)
  inTrain <- createDataPartition(y = classi, p = .70, list = FALSE) # 70% des 
  training <- dada[inTrain,  6:2005]
  testing <- dada[!(inTrain),  6:2005]
  trainclass <- as.factor(classi[inTrain]); testclass <- as.factor(classi[!(inTrain)])
  
  finalModel <- plsda(training, trainclass, ncomp = n_comps, probMethod = "Bayes", method = "simpls") 
  fin_mods_tous[[i]] <- finalModel  
} 

head(fin_mods_tous[[1]])

# saveRDS(fin_mods_tous, "~/Documents/Maîtrise/DonnéesAnalyses/PLS/DA/savedRDS/ttsp_50iter_normalized_PLSDA_fin_mods.rds")
# 4 juillet : 70%-30%, 100 itérations, 4 n_comps / 2 août : 70-30%, 50 itérations, 4 n_comps
# 5 août : 70-30%, 50 iter, 4 n_comps, NORMALISÉ


#=======================================# VISUALISATION #=======================================#                    ----

#-=-=-=-=-=-=-=-=-=-=-=-=-
# graphique scores et classement (3D)                                                                         ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

# choix des couleurs !
col_ord <- c('#313695','#4575b4','#74add1','#abd9e9',"darkorchid4",'yellow2','#a50026','#d73027','#f46d43',
             '#fdae61',"grey50","black")
#col_ordi <- rep(col_ord[1:11], 11) #as.numeric(table(out$species))
#col_ordi <- rev(col_ordi)
# i = 1
out <- list()

for (i in 1:n_iter) {
  out[[i]] <- as.matrix(fin_mods_tous[[i]]$scores[1:68, 1:ncol(fin_mods_tous[[i]]$scores)]) # nombre de lignes = 12 (50%)
  out[[i]] <- as.data.frame(out[[i]])
  colnames(out[[i]]) <- gsub(pattern = " ", "_", colnames(out[[i]])) # " " à remplacer par "-" dans les colonnes créées
  out[[i]]$species <- as.character(trainclass) 
  out[[i]] <- out[[i]][order(out[[i]]$species, decreasing = T),]
  col_ordi <- rep(col_ord[1:12], length.out = 17)#lenght.out = as.numeric(table(out[[i]]$species)))
  out[[i]]$col_ordi <- rev(col_ordi)
}


head(out[[1]])  

# EXPLICATIONS ----
# √ pour les 50 itérations, donnes-moi le score des lignes 1:12 et de toutes les colonnes 
# "scores describe the position of each sample in each determined latent variable and weights describe the contribution of each variable to each LV "
# √ remplace les espaces dans out[[i]]
# √ crée une colonne avec les noms d'espèce
# √ met les lignes en ordre croissant 
# √ met les couleurs précisées ci-haut, répète les couleurs autant de fois qu'il y a d'items par espèce (species) = 3 dans ce cas

# affichage                                                                       COMMENT INTERPRÉTER ?
plot_ly(out[[1]], x = ~ Comp_2, y = ~ Comp_3, z = ~ Comp_4, type = "scatter3d", mode = "markers", # z = ~Comp_3
        color = ~ species, marker = list(size = 9, opacity = 0.8), # marquer : aes() pour les points
        colors = col_ordi) 
# distinction pas claire entre les espèces dans cet espace... essayer différentes variable latentes pour voir
# ERIOPHORUM dans son petit coin tho

max(fin_mods_tous[[i]]$scores) # vérifier si score bizzare pour certains échantillons


#-=-=-=-=-=-=-=-=-=-=-=-=-
# code pour graphique de probabilité et de confusion                                                              ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-
probis <- list()
confus <- list()

# we want probabilities and predicted classes using only *testing*
# probis : nbr between 0:100 % 
# and predicting class membership

for (i in seq(n_iter)){ 
  print(i)
  flush.console()
  set.seed(i)

  inTrain <- createDataPartition(y = classi, p = .70, list = FALSE)
  
  testing <- dada[-inTrain, 6:2005 ] 
  testclass <- as.factor(classi[-inTrain])
  
  plsProbs <- predict(fin_mods_tous[[i]], newdata = testing, type = "prob") # compare les résultats aux probabilités obtenues (?) // predict() / newdata = specifying the first place to look for explanatory variables to be used for prediction.
  # plsProbs -> par ligne (dans testing) on a la probabilité que chaque échantillon appartient à chaque classe
  # vérifier qquns pour voir si ça fonctionne 
  plsClasses <- predict(fin_mods_tous[[i]], newdata = testing, type = "class") # compare les résultats aux facteurs "testing"
  confus[[i]] <- confusionMatrix(data = plsClasses, testclass) # {caret} reference ~ predicted (same levels*) données plsClsses
  
  probs <- as.data.frame(plsProbs)
  names(probs) <- sapply(strsplit(names(probs), split = ".n"),"[",1) # sapply() pour des vecteurs
  probs <- cbind(testclass, probs)
  probis[[i]] <- probs 
}
# saveRDS(confus, "savedRDS/confus_final.rds") # long, fait 29 août

head(confus[[1]]) # lists of every statistics we calculated 

confus[[1]]$overall

# NOTE :
# Overall Statistics
# confus[[1]]$overall # 100 itérations
# Accuracy          Kappa       AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 9.782609e-01   9.709962e-01   8.847282e-01   9.994498e-01   2.608696e-01   1.878677e-25            NaN 
# is the model good YES // is it trustable ? oui
# confus[[1]]$overall # 50 itérations
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 9.782609e-01   9.709962e-01   8.847282e-01   9.994498e-01   2.608696e-01   1.878677e-25            NaN 
# confus[[1]]$overall # 50 iterations, NORMALIZED SAME SAME
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue  McnemarPValue 
# 9.782609e-01   9.709962e-01   8.847282e-01   9.994498e-01   2.608696e-01   1.878677e-25           
# Accuracy          Kappa  AccuracyLower  AccuracyUpper   AccuracyNull AccuracyPValue 
# 9.230769e-01   8.972332e-01   7.486971e-01   9.905446e-01   2.692308e-01   3.760347e-12 

#-=-=-=-=-=-=-=-=-=-=-=-=- 
# affichage graphique de probabilité                                                                                                                                   ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

arr <- array(unlist(probis), dim = c(dim(probis[[1]]), n_iter)) # faire un array avec la liste, dont la dimension est de 46:5:50
prob_mean <- apply(arr, 1:2, mean) # faire la moyenne avec les dimensions 1 (échantillons testé) et 2 (5 variables latentes) (sur 3 dimmentions,50 = n_iter)
                                   # / chaque combinaison (i,j) pareille = moyenne, toutes listes confondues
                                   # exemple : z <- array(1:24, dim = 2:4) / zseq <- apply(z, 1:2, mean) / zseq
prob_mean <- as.data.frame(prob_mean) # on n'a plus que 2 dim (échantillon x son classement dans testclass, toutes itération confondues)
                                      # mais la moyenne de la colonne 1 a pas rapport

prob_mean$V1 <- probis[[2]]$testclass   # remplaces V1 qui correpond à comment ça a été classé avec le modèle #1
colnames(prob_mean) <- colnames(probis[[1]]) 

pp <- reshape::melt(prob_mean, id.vars = 'testclass')

# vérifier si facteurs sont les mêmes
levels(pp$variable)
levels(pp$testclass) # nope

# demander d'écrire les facteurs pareil
pp$variable <- revalue(pp$variable, # revalue() => pkg "plyr"
                       c("Eriophorum vag" = "Eriophorum vaginatum",
                         'Kalmia ' = 'Kalmia angustifolia',
                         'Rhodod' = 'Rhododendron groenlandicum',
                         'Chamaedap' = 'Chamaedaphne calyculata'))

# vérifier si facteurs sont les mêmes
levels(pp$variable)
levels(pp$testclass) # yes

# ajouter colonne 'position'
pp$position <- ifelse(pp$testclass == pp$variable, 2, 1) # fait une colonne dans laquelle valeur son définies par : test = si pp$testclass est égal à pp$variable, si oui position = 2, sinon position = 1

# mettre en ordre
pp$testclass <- factor(pp$testclass, levels = rev(levels(pp$testclass)))

# choix des couleurs !
coli <- c('#313695','#4575b4','#74add1','#abd9e9',"darkorchid4",'yellow2','#a50026','#d73027','#f46d43',
          '#fdae61',"grey50","black")

# affichage graphique de probabilité
#pdf("./R_output/probability_plot_leaves_20comps.pdf",width = 6,height = 4)
ggplot(pp, aes(x = testclass, y = value, fill = variable, group = position)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = alpha(coli, 1)) +
  theme(panel.background = element_blank()) + # legend.title = element_blank()) +
  labs(x = "Species", y = "", title = "Probabilité de classement ; paramètres = 
      refl_sg_large_famille_sp_obj_normalized / normalized/ fait 6/08/19,
      70-30%, 50 iter, 4 n_comps, enregistré 5 août") +
  theme_bw() +
  coord_flip()
#ggsave('probabilité_PLSDA_50iter_normalized_ttsp.pdf', height = 10, width = 12)
#dev.off()


#-=-=-=-=-=-=-=-=-=-=-=-=- 
# affichage graphique de confusion                                                                                                                                   ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

tabs <- list()

for(i in 1:length(confus)){
  tabs[[i]] <- confus[[i]]$table
}

tabsi <- Reduce('+', tabs)
tab_mean <- as.data.frame.matrix(tabsi/length(confus))
# write.csv(tab_mean,"/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/DA/R_output/PLSDA_confumean_espèces_4comps.csv")

sums <- colSums(tab_mean)
tabs_perc <- matrix(NA, length(sums),length(sums))
for (i in 1:length(sums)){
  tabs_perc[,i] <- (tab_mean[,i]/sums[i]*100)
}

colnames(tabs_perc) <- colnames(confus[[1]]$table)
rownames(tabs_perc) <- rownames(confus[[1]]$table)
# write.csv(tabs_perc,"/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/DA/R_output/PLSDA_confuperc_especes_4comps.csv")

col <- colorRampPalette(c("black", "black", "brown", "gold", "forestgreen")) 

corrplot::corrplot(tabs_perc, p.mat = tabs_perc, insig = "p-value", sig.level = -1, addCoef.col = 1,
         tl.srt = 70, col = col(20), cl.lim = c(0, 1), tl.col = 1, tl.offset =1.5, 
         cl.ratio = 0.2, cl.align.text = "l", cl.cex = 0.7, tl.cex = 0.6,
         mar = c(1,3,3,3), title = "Probabilité de classement ; paramètres = \n
         refl_sg_large_famille_sp_obj_normalised \n 
         fait 6/08/19, 70-30%, 50 iter, 4 n_comps, enregistré 6 août")
mtext("Prediction", 2, line = 0, cex = 1.2)
mtext("Reference", at = 2, line = -6, cex = 1.2)
dev.off()


#================================# Statistiques finales du modèle PLS-DA  #================================#                   ----

#-=-=-=-=-=-=-=-=-=-=-=-=- 
# justesse du modèle                                                          ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

# changer ncomps
accu <- numeric(length = length(mods_tous))
for (i in 1:length(mods_tous)){
  accu[i] <- mods_tous[[i]]$results$Accuracy[n_comps]
}

(accmean <- mean(accu))
(accsd <- sd(accu))

# justesse du modèle  
kappas <- numeric(length=length(mods_tous))
for (i in 1:length(mods_tous)){
  kappas[i] <- mods_tous[[i]]$results$Kappa[n_comps]
}

(kappmean <- mean(kappas))
(kappsd <- sd(kappas))


#-=-=-=-=-=-=-=-=-=-=-=-=- 
# importance des bandes spectrales                                                          ----        
#-=-=-=-=-=-=-=-=-=-=-=-=-

# --->>> modèle final <<<---
fin_mods_tous <- readRDS("savedRDS/ttsp_50iter_normalized_PLSDA_fin_mods.rds") # enregistré le ? / refait le 6 août avec vector-normalized

# créer une liste des rowsums of loadings (extracting the loadings)
lls <- list() 
for(i in 1:length(fin_mods_tous)){                # CI_DESSOUS : 2:n_comps PCQ norm_magn 1 ière colonne, ne doit pas être incluse !
  lls[[i]] <- abs(loadings(fin_mods_tous[[i]])[1:dim(loadings(fin_mods_tous[[1]]))[1],2:n_comps]) # loadings: evry comp AND evry wvl ** PAS NOR_MAGNITUDE
  sumis <- lapply(lls, rowSums) # j'ai choisi 2 n_comps (dimensions : number of wvlght) pour toute dimention donne 
} # lapply : appliquer à tout le modèle

mm <- apply(simplify2array(sumis), 1, mean) # mean of loadings (of the wavelenghts) / importance in differencing
ss <- apply(simplify2array(sumis), 1, sd) # standard deviation 

mm <- as.data.frame(mm)
ldngs <- cbind(mm,ss)
ldngs <- rownames_to_column(ldngs)
ldngs[,1] <- as.character(ldngs[,1])
colnames(ldngs) <- c("wvl", "mean", "sd")
# write.csv(ldngs, file = paste("/Users/Aliz/Documents/Maîtrise/DonnéesAnalyses/PLS/DA/", n_comps, "ncomps_ldngs.csv", sep = ""), row.names = FALSE)
# enregistrer pour pouvoir générer beau graph dans ggplot, voir document PLS, "script VIP&var.importance.R"

# afficher le graph d'importance des bandes
plot(1:ncol(testing), mm$mm, type = "n", bty = "l", ylab = "abs (loadings)", 
     xaxt = "n", xlab = "Wavelength (nm)", xaxt = "n", xlim = c(300,2500))