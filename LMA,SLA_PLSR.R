# 12 juillet 2019

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#                                   PLS-R structure
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NOTES : voir RESOURCES dans document "PLS" pour scripts d'Anna

# GROUPE TF "aire et masse" : "PLATEAU NIR" 700-1350 (Asner et al., 2011)


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Libraries ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

library(dplyr)
library(spectrolab)
library(caret) # Classification And REgression Training
library(reshape) # Classification And REgression Training
library(agricolae)
library(pls)

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Données ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# définir le répertoire
setwd("~/Documents/Maîtrise/DonnéesAnalyses/PLS/R")
# setwd("~/Desktop/Alizée")

load("~/Documents/Maîtrise/DonnéesAnalyses/PLS/rt_sg_large.RData") # rt_sg_large
# load("rt_sg_large.RData") # rt_sg_large
# obtenu via manips sur script "nettoyage_spectres_lissage_correction.R" dans le document "scripts-NETTOYAGE"

load("~/Documents/Maîtrise/DonnéesAnalyses/traits_fonctionnels/tous_analyses_vX2.RData") # ft_all_vX2
# load("tous_analyses_vX2.RData") # ft_all_vX2
ft_all <- ft_all_vX2
# obtenu via script "nettoyage-TF-enlever.ouliers" document scripts-NETTOYAGE


#=======================================# Manipulations préalables #=======================================#                    ----

#-=-=-=-=-=-=-=-=-=-=-=-=-
# obtenir df pour analyses ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

refl_ft <- rt_sg_large %>%  # autres colonnes intéressantes, voir colnames(rt_sg_large)
  select(c("700":"1350")) %>%  
  dplyr::rename_all(function(x) paste0("X", x)) %>%
  cbind(rt_sg_large[, c(1:2, 4:5)]) %>%
  filter(propriete == "reflectance") %>%
  full_join(ft_all) %>%
  dplyr::select("sample_id", "scientific_name", "site_class", 'cellulose', 'leaf_dry_matter_content_mg_g', 
                'N_pourc', 'chlA_mg_g', 'Tannins_mg_g', 'lignine', 'leaf_water_content_mg_g', 'C_N_ratio',
                "plot_id", 'soluble_contents', 'Phenols_mg_g', 'C_pourc', 'equivalent_water_thickness_cm', 
                'chlB_mg_g', 'leaf_mass_per_area_g_m2', "specific_leaf_area_m2_kg", 'carotenoides_mg_g', 
                "site_id", 'chlA_chlB', 'recalcitrant', "hemicellulose", 'chlA_mg_m2', 'chlB_mg_m2',
                "propriete", "X700":"X1350")
#  188 obs √ 94 c'est 95 - 1 sample rejeté !
# ATTENTION ! échantillon # 12186521 manque FT necess poudre :  N%, C_N_ratio, CARBONE fr, C%, ect
# l'enlever pour ces analyses et/ou ne pas me demander pourquoi yen manque 1 dans les résultats !


#-=-=-=-=-=-=-=-=-=-=-=-=-
# ENLEVER ÉCHANTILLON                                                                        ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# refl_ft_exper <- refl_ft_exper %>% 
#   dplyr::filter(!sample_id == 12186521) # manque les FT nécessitant le sample BROYÉ : N%, C_N_ratio, CARBONE fr, C%, ect


#                               ###############################                               #
#                                  #########################                                  #
#                                     ###################                                     #
#                                        #############                                        #
#                                           #######                                           #
#                                              #                                              #


#=======================================# *Iterative* PLS-R #=======================================#                                                       ----

#                  #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#
#                  #               PARTIE 1 : trouver le nombre de composants                 #                  #           -----
#                  #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#

# NOTE : rouler le modèle itératif tous traits des mm wvls ensemble
#        rouler LA PARTIE : "" sur les traits séparés ! faire nouvelle matrice à ce moment là


#-=-=-=-=-=-=-=-=-=-=-=-=-
# créer matrice d'analyse pour modèle itératif ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

# créer dataframe et selectionner les traits à rouler ensemble
aire_masse_df <- refl_ft %>%
  dplyr::select(sample_id, leaf_mass_per_area_g_m2, specific_leaf_area_m2_kg, "X700":"X1350")

# créer objets utilisés dans boucle 
inBands <- names(aire_masse_df)[4:654]
inVars <- names(aire_masse_df)[2:3]

set.seed(1840) # pour obtenir tjrs les mm données mm si itérations roulées plusieurs fois
for (i in 1:length(inVars)){
  inVar <- inVars[i]
  iForm <- paste(inBands, collapse = "+")
  iForm <- as.formula(paste(inVar, "~", iForm))
  # dir.create(paste("~/Desktop/", inVar, sep = ""), recursive = T)
  
  # EXPLICATIONS ----
  # la formule ci-haut demande 
  # √  pour toutes les variables dans inVars (voir ci-haut)
  # √  mets des + entre les wvl et fait une relation (~) avec la i variable
  # √  crée un répertoire sur le bureau dans lequel mettre les données
  # ----
  
  # boucle pour trouver le nombre de composants
  nCut <- floor(0.7 * nrow(aire_masse_df)) # 80% data for calibration
  nComps <- 15
  nsims <- 100                          
  outMat <- matrix(data = NA, nrow = nsims, ncol = nComps) #
  outMatRM <- matrix(data = NA, nrow = nsims,ncol = nComps) # RM = root mean ?
  
  for (nsim in seq(nsims)){
    print(nsim) # écrit tes rendu à quelle itération
    flush.console() # ?
    aire_masse_df$RAND <- order(runif(nrow(aire_masse_df))) # runif(n, min = , max = ) is used to generate n uniform (autant de fréquence de chaque tranche) random numbers lie in the interval (min, max).
    subData <- aire_masse_df[aire_masse_df["RAND"] < nCut, ] # choisi ligne < à ordre random, subset pour tester
    resNCOMP <- plsr(iForm, data = subData, ncomp = nComps, validation = "LOO", method = "oscorespls")
    resPRESS <- as.vector(resNCOMP$validation$PRESS) # créer objet resPRESS et lui donner "n" places
    outMat[nsim, seq(resNCOMP$validation$ncomp)] <- resPRESS # met objet resPRESS dans outMat pour chaque itération
    resRMSEP <- as.numeric(RMSEP(resNCOMP, estimate = "CV", intercept = F)$val) # RMSEP : mvr object
    outMatRM[nsim, seq(resNCOMP$validation$ncomp)] <- resRMSEP # à chaque itération, met les résultats de resRMSEP dans outMat
  }
  
  # EXPLICATIONS ----
  # la formule ci-haut demande 
  # √  prends 80% des données, 15 comps et fait 100 it. 
  # √  fait deux matrices "outMat" et "outMatRM", pour l'instant, avec des NA et tant de colonnes et lignes
  # √  pour chaque itération, crée un colone "RAND" qui correspond à un mélange de données (selon chiffres assignés random et ordonnés)
  # √  ensuite prends en X nombre (80%) pour tester (subData)
  # √  fait une plsr dessus qui va s'appeller resNCOMP
  # √  fait ensuit resPRESS : contient le chiffre qu'on a donné pour ncomps 
  # √  dans outMat, met -> nsim : no de ligne du no d'itération + seq(resNCOMP$validation$ncomp) : soit 15 colonnes 
  # √  pour chaque itération, met les résultats resPRESS à la ligne correspondant au no d'itération dans outMat
  # √  calcule RMSEP(), renvoie les $val et crée resRMSEP
  # √  pour chaque itération, met les résultats resRMSEP à la ligne correspondant au no d'itération dans outMatRM
  # ----
  
  # "PRESS" stat: test de Tukey et visualisation                      où PRESS : Predictive Residual Sum of. Squares
  pressDF <- as.data.frame(outMat) # crée pressDF vide
  names(pressDF) <- as.character(seq(nComps)) # créé ligne 239
  pressDFres <- melt(pressDF) # forme accetables pour faire un modèle linaire
  
  modi <- lm(value ~ variable, pressDFres) # modèle linéaire puis test Tuckey
  tuk1 <- HSD.test(modi, trt = "variable")                        
  tuk_dat1 <- as.data.frame(tuk1$groups)
  tuk_dat1$var <- as.numeric(row.names(tuk_dat1))
  tuk_dat1 <- tuk_dat1[order(tuk_dat1$var, decreasing = F),]
  letters <- as.character(tuk_dat1$groups)
  
  jpeg(paste("resultats/", "800-2400_", inVar, "/", inVar, "_abs_PRESS.jpg", sep = ""),
       width = 6, height = 5, units = "in", res = 200)
  par(bty = "l")
  boxplot(pressDFres$value ~ pressDFres$variable,
          xlab = "n Components", ylab = "PRESS", main = inVar)
  text(x = 1:max(as.numeric(pressDFres$variable)), 
       y = rep(max(pressDFres$value),15), letters)
  dev.off()
  
  # RMSEP : Test de Tukey et visualisation
  RMDF <- as.data.frame(outMatRM)
  names(RMDF) <- as.character(seq(nComps))
  RMDFres <- melt(RMDF)
  
  modi <- lm(value ~ variable, RMDFres) 
  tuk2 <- HSD.test(modi, trt = "variable")
  tuk_dat2 <- as.data.frame(tuk2$groups)
  tuk_dat2$var <- as.numeric(row.names(tuk_dat2))
  tuk_dat2 <- tuk_dat2[order(tuk_dat2$var,decreasing = F),]
  letters <- as.character(tuk_dat2$groups)
  
  jpeg(paste("resultats/", "800-2400_", inVar, "/", inVar,"_abs_RMSEP.jpg",sep = ""),
       width = 6, height = 5, units = "in", res = 200)
  par(bty = "l")
  boxplot(RMDFres$value ~ RMDFres$variable,
          xlab = "n Components", ylab = "RMSEP", main = inVar)
  text(x = 1:max(as.numeric(RMDFres$variable)), 
       y = rep(max(RMDFres$value), 15), letters)
  dev.off()
  
  # EXPLICATIONS ----
  # les formules ci-haut demandent 
  # √  de calculer les valeurs PRESS et RMSEP
  # √  et de renvoyer un boxplot dans le répertoire créé
  # ----
  
}

saveData <- aire_masse_df  



#                     #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#
#                     #                                  PARTIE 2-1                                #           ----- 
#                     #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#

# leaf_mass_per_area_g_m2


#-=-=-=-=-=-=-=-=-=-=-=-=-
# Fonction d'Anna                                                          ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

VIPjh = function(object, j, h) {    
  if (object$method !="oscorespls") stop("Only implemented for orthogonal scores algorithm. Refit with 'method = \"oscorespls\"'")
  if(nrow(object$Yloadings) > 1) stop("Only implemented for single-response models")
  b = c(object$Yloadings)[1:h]
  T = object$scores[,1:h, drop = FALSE]
  SS = b ^ 2 * colSums(T ^ 2)
  W = object$loading.weights[ ,1:h, drop = FALSE]
  Wnorm2 = colSums(W ^ 2)
  VIP = sqrt(nrow(W) * sum(SS * W[j,] ^ 2 / Wnorm2) / sum(SS)) 
  return(VIP)
}

# CHOIX DU TRAIT, un à la fois * 
(inVar <- names(aire_masse_df[2])) 

# subsetter saveData (issu de la partie 1)
aire_masse_df_ITER <- saveData[complete.cases(saveData[,inVar]),]

# choix de n comps selon partie 1
nCompss <- 6

# formule de la PLS-R
iForm <- paste(inBands, collapse = "+")
iForm <- as.formula(paste(inVar, "~", iForm))

set.seed(1840)
for (i in 1:length(nCompss)){
  nsims <- 100
  nComps <- nCompss[i]
  nCut <- floor(0.7*nrow(aire_masse_df_ITER))
  
  # créer objets vides
  coefMat <- matrix(data = NA, nrow = nsims, ncol = length(inBands) + 1)
  coefStd <- matrix(data = NA, nrow = nsims, ncol = length(inBands) + 1)
  vipMat <- matrix(data = NA, ncol = nsims, nrow = length(inBands))
  statMat <- matrix(data = NA, nrow = nsims, ncol = 6)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    aire_masse_df_ITER$RAND <- order(runif(nrow(aire_masse_df_ITER))) # same de partie 1 : juste créer chiffre random sur elsquels baser les subsets test et training
    subData <- aire_masse_df_ITER[aire_masse_df_ITER["RAND"] < nCut,]
    tstData <- aire_masse_df_ITER[aire_masse_df_ITER["RAND"] >= nCut,] # ">=" ne pas séparer
    
    resX <- plsr(iForm, data = subData, ncomp = nComps, method = "oscorespls") 
    resS <- plsr(iForm, data = subData, ncomp = nComps, method = "oscorespls", scale = T) # scaled : POURQUOI ? "importance des wvls"
    
    # Coefficients (brut et standardisé)
    coefs <- as.vector(coef(resX, ncomp = nComps, intercept = T))
    zcoef <- as.vector(coef(resS, ncomp = nComps, intercept = T))
    
    coefMat[nsim,] <- coefs # à chaque ligne = no itération, ajoute les coeficients dans les colonnes respectives
    coefStd[nsim,] <- zcoef # standardized coeffis for ***importance of wvls*** /// SAME action
    
    # VIP
    vip <- c() # crée un objet "vip"
    for (j in seq(length(inBands))){ # 2000 bandes
      vip <- c(vip, VIPjh(resS, j, nComps)) # crée vip qui est un regroupement des vip précédantes, et du résultat VIPjh
      # où j = no de bande et h = le no de variables latentes choisies
    }
    
    vipMat[,nsim] <- vip # met "vip" dans vipMat[toutes lignes, colonne corr au no de simulation]
    
    # EXPLICATIONS ----
    # la formule ci-haut demande 
    # √ pas rapport que ça se fasse là; aurait cert. pu se faire n'importe quand avant lign 473... (?)
    # √ calcule les vip pour chaque bande dans chaque variable latente (dans chq boucle itérative)
    # √ crée vipMat (réutilisé ensuite, ligne 473 **)
    # ----
    
    # Statistiques du modèle
    fitX <- as.vector(unlist(resX$fitted.values[, 1, nComps])) # toutes lignes, 1iere colonne, nComps auquel on est rendu (nComps correspond au NO DE LISTE)
    preX <- as.vector(unlist(predict(resX, ncomp = nComps, newdata = tstData))) # prédit les résultats fittés à partir du set "testing"
    fitBias <- mean(as.numeric(subData[, inVar]) - fitX) # moyenne de l'erreur de chaque ligne (échantillon)
    valBias <- mean(as.numeric(tstData[, inVar]) - preX) # SAME, mais avec les valeurs prédites
    subData_df <- data.frame(subData[,])
    tstData_df <- data.frame(tstData[,])
    fitR2 <- summary(lm(subData_df[, inVar] ~ fitX))$r.squared
    valR2 <- summary(lm(tstData_df[, inVar] ~ preX))$r.squared
    fitRMSE <- sqrt(mean((as.numeric(subData[,inVar]) - fitX) ^ 2))
    valRMSE <- sqrt(mean((as.numeric(tstData[,inVar]) - preX) ^ 2))
    outVec <- c(fitR2, fitRMSE, fitBias, valR2, valRMSE, valBias)
    statMat[nsim,] <- outVec
  }
  
  statMat <- as.data.frame(statMat)
  names(statMat) <- c("fitR2","fitRMSE","fitBias","valR2","valRMSE","valBias")
  write.csv(statMat, paste("resultats/calcul_NRMSE/", inVar, "_", nComps,"comps_abs_stats.csv", sep=""), row.names=FALSE)
  write.csv(statMat, paste("resultats/", "800-2400_", inVar, "/", inVar, "_", nComps,"comps_abs_stats.csv", sep=""), row.names=FALSE)
  
  # EXPLICATIONS ----
  # la formule ci-haut demande 
  # √ met dans "statMat" les valeurs "outVec" calculées, mais à la ligne corr. au no d'intération
  # √ valeurs outVec : fitR2, fitRMSE, fitBias, valR2, valRMSE et valBias,
  # √ calculées à partir de resX et de prédictions subséquentes 
  # √ enregistre statMat
  # ----
  
  # Coefficients du modèle
  coeffis <- data.frame(matrix(nrow = length(inBands) + 1, ncol = 3))
  names(coeffis) <- c("bands", "mean","stdv")
  
  # entrer les données nécessaires dans le df créé
  coeffis$bands <- c("Intercept",inBands) # crée un vecteur où lignes = ça
  coeffis$mean <- apply(coefMat, MARGIN = 2, FUN = mean) # coefMat : coefficients de la plsr lignes 377, 378
  coeffis$stdv <- apply(coefMat, MARGIN = 2, FUN = sd)
  
  # enregistrer
  # write.table(coeffis, paste("~/Desktop/", inVar, "/", inVar, "_", nComps, "comps_abs_coeffMEAN.csv", sep = ""), 
  #             sep = ",", col.names = T, row.names = F)
  
  
  # Prédictions du modèle
  specMat <- aire_masse_df_ITER[, inBands] # ne sort que les colonnes inBands d'aire_masse_df_ITER
  specMat <- cbind(rep(1, nrow(specMat)), specMat) # colles-y le nombre de colonne = nb échantillons + le specMat existant (itérations)
  specMat <- as.matrix(specMat)
  
  predMat <- specMat %*% t(coefMat) # {base} x %*% y = matrix mutliplication // multiplie matrice specMat par transpose de coefMat (voir Legendre)
  predMean <- apply(predMat, FUN = mean, MARGIN = 1) # 
  predStdv <- apply(predMat, FUN = sd, MAR = 1)
  
  preds <- aire_masse_df_ITER[, !(names(aire_masse_df_ITER) %in% inBands | names(aire_masse_df_ITER) == "RAND")] # 2 critères de sélection des colonnes d'aire_masse_df_ITER à mettre dans
  # preds : [ {base} "x | y" = "OR" ] / DONC pas les noms d'aire_masse_df_ITER 
  # qui sont dans inBands OU METTRE les colonnes names(aire_masse_df_ITER) == "RAND")...
  # pourquoi référer à des lignes ? (RAND)
  preds[, paste("predMean_", inVar, sep = "")] <- predMean # ajoute colonne predMean_Var à df preds
  preds[, paste("predStdv_", inVar, sep = "")] <- predStdv # same, pr chaque Var
  write.csv(preds,paste("resultats/pred_meas/", inVar, "_", nComps,"comps_abs_preds.csv", sep=""), row.names=FALSE)
  write.csv(preds, paste("resultats/", "800-2400_", inVar, "/", inVar, "_", nComps,"comps_abs_preds.csv", sep=""), row.names=FALSE)
  
  # Visualisation des prédictions 
  modCI <- quantile(statMat$fitR2, probs = c(0.05, 0.95)) # CI : confidence intervals
  
  formi <- as.formula(paste(paste(inVar, " ~ predMean_", inVar, sep = ""))) # pourquoi paste(paste()) ?
  lmformi <-   as.formula(paste(paste("predMean_", inVar, " ~ ", inVar, sep = "")))
  
  jpeg(paste("resultats/", "800-2400_", inVar,  "/", inVar, "_", nComps, "comps_abs_predplot.jpg", sep = ""),
       width = 6, height = 5, units = "in", res = 300)
  plot(formi, data = preds, pch = 16, cex = 0.8, ylab = "measured", xlab = "predicted",
       main = inVar, xlim = c(min(predMean - predStdv), max(predMean + predStdv)))
###  abline(lm(lmformi, data = preds))
  abline(a = 0, b = 1, lty = 2)
  # arrows(predMean, aire_masse_df_ITER[, inVar], predMean + predStdv, aire_masse_df_ITER[, inVar], angle = 90, length = 0.05, lwd = 0.8)
  # arrows(predMean, aire_masse_df_ITER[, inVar], predMean - predStdv, aire_masse_df_ITER[, inVar], angle = 90, length = 0.05, lwd = 0.8)
  legend("topleft", bty = "n", cex = 0.8,
         c(paste("R² = ", sprintf("%.2f", signif(mean(statMat$fitR2), 3)), " [", signif(modCI[1], 2), ",", signif(modCI[2], 2), "]", sep = ""), 
           paste("RMSEP =", sprintf("%.2f", signif(mean(statMat$fitRMSE), 3)), sep = " "),
           paste("ncomps =", nComps, sep = " ")))
  dev.off()
  
  # Visualisation des VIP ("vip aggrégé")
  vipAggr <- as.data.frame(t(apply(vipMat, MARGIN = 1, FUN = quantile, probs = c(0.05,0.5,0.95)))) # crée df de la transpose (dcast genre) des quantiles de la vipMat
  # où vipMat = vip pour chaque bande dans chaque variable latente (dans chq boucle itérative)
  
  vipAggr$mean_VIP <- apply(vipMat, MARGIN = 1, FUN = mean) # ajoute à vipAggr la moyenne des VIP
  vipAggr$stdv <- apply(vipMat, MARGIN = 1, FUN = sd) # same, std
  stderr <- function(x) sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))  # standard error of mean 
  vipAggr$std_err <- apply(vipMat, MARGIN = 1, FUN = stderr)
  vipAggr$band <- inBands # ajoute un colonne où ligne = identifiant de la bande correspondante
  # vipAggr utilisé dans visualisation, ligne 512
  vipAggr$band <- gsub("X", "", vipAggr$band)
  vipAggr$band <- as.numeric(vipAggr$band)
  write.csv(vipAggr,  file = paste("resultats/vip/", inVar, "_", "vipAggr.csv", sep = ""), row.names = FALSE)
  
  # Coefficients standardisés pour la visualisation 
  coeff_std <- data.frame(matrix(nrow = length(inBands) + 1, ncol = 3))
  names(coeff_std) <- c("bands", "mean", "stdv")
  
  # ajoute au df coeff_std les valeurs suivantes, tirées de la plsr
  coeff_std$bands <- c("Intercept", inBands) # id bande
  coeff_std$mean <- apply(coefStd, MARGIN = 2, FUN = mean) # coefStd créé ligne 385
  coeff_std$stdv <- apply(coefStd, MARGIN = 2, FUN = sd) # same
  
  # # Visualisation VIP et des coefficients standardisés 
  # jpeg(paste("resultats/", "800-2400_", inVar,  "/", inVar, "_", nComps, "comps_abs_varimp.jpg", sep = ""),
  #      width = 6, height = 7, units = "in", res = 300)
  # par(mfrow = c(2,1), mar = c(1.5,4,2.5,1.5), oma = c(3,0,0,0))
  # plot(coeff_std$mean[-1] ~ as.numeric(substr(coeff_std$bands[-1], 2, nchar(coeff_std$bands[-1]))),
  #      type = "p", pch = 19, xlab = "", ylab = "coeff_stdmean", main = paste(inVar, nComps, "comps", sep = "_"),
  #      ylim = c(-max(abs(coeff_std$mean[-1])), max(abs(coeff_std$mean[-1]))), bty = "l")
  # abline(h = 0)
  # points(abs(coeff_std$mean)[-1] ~ as.numeric(substr(coeff_std$bands[-1], 2, nchar(coeff_std$bands[-1]))),
  #        xlab = "wvl", ylab = "coeff_stdmean", col = 2, pch = 16, cex = 0.8)
  # 
  # lines(abs(coeff_std$mean)[-1] ~ as.numeric(substr(coeff_std$band[-1], 2, nchar(coeff_std$band[-1]))), col = 2)
  # 
  # plot(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))), vipAggr$mean_VIP, type = "l",
  #      xlab = "wvl", ylab = "VIP", bty ="l")
  # polygon(x = c(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))),
  #               rev(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))))),
  #         y = c(vipAggr$mean_VIP + vipAggr$stdv * 1.96, rev(vipAggr$mean_VIP - vipAggr$stdv * 1.96)),  
  #         col =  adjustcolor("red", alpha.f = 0.2), border = NA)
  # mtext("wavelength(nm)", 1, outer = T, line = 1)
  # dev.off()
}


#                     #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#
#                     #                                  PARTIE 2-2                               #           ----- 
#                     #----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#

# "specific_leaf_area_m2_kg"

#-=-=-=-=-=-=-=-=-=-=-=-=-
# Fonction d'Anna                                                          ----
#-=-=-=-=-=-=-=-=-=-=-=-=-

VIPjh = function(object, j, h) {    
  if (object$method !="oscorespls") stop("Only implemented for orthogonal scores algorithm. Refit with 'method = \"oscorespls\"'")
  if(nrow(object$Yloadings) > 1) stop("Only implemented for single-response models")
  b = c(object$Yloadings)[1:h]
  T = object$scores[,1:h, drop = FALSE]
  SS = b ^ 2 * colSums(T ^ 2)
  W = object$loading.weights[ ,1:h, drop = FALSE]
  Wnorm2 = colSums(W ^ 2)
  VIP = sqrt(nrow(W) * sum(SS * W[j,] ^ 2 / Wnorm2) / sum(SS)) 
  return(VIP)
}

# CHOIX DU TRAIT, un à la fois * 
(inVar <- names(aire_masse_df[3])) 

# subsetter saveData (issu de la partie 1)
aire_masse_df_ITER <- saveData[complete.cases(saveData[,inVar]),]

# choix de n comps selon partie 1
nCompss <- 6

# formule de la PLS-R
iForm <- paste(inBands, collapse = "+")
iForm <- as.formula(paste(inVar, "~", iForm))

set.seed(1840)
for (i in 1:length(nCompss)){
  nsims <- 100
  nComps <- nCompss[i]
  nCut <- floor(0.7*nrow(aire_masse_df_ITER))
  
  # créer objets vides
  coefMat <- matrix(data = NA, nrow = nsims, ncol = length(inBands) + 1)
  coefStd <- matrix(data = NA, nrow = nsims, ncol = length(inBands) + 1)
  vipMat <- matrix(data = NA, ncol = nsims, nrow = length(inBands))
  statMat <- matrix(data = NA, nrow = nsims, ncol = 6)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    aire_masse_df_ITER$RAND <- order(runif(nrow(aire_masse_df_ITER))) # same de partie 1 : juste créer chiffre random sur elsquels baser les subsets test et training
    subData <- aire_masse_df_ITER[aire_masse_df_ITER["RAND"] < nCut,]
    tstData <- aire_masse_df_ITER[aire_masse_df_ITER["RAND"] >= nCut,] # ">=" ne pas séparer
    
    resX <- plsr(iForm, data = subData, ncomp = nComps, method = "oscorespls") 
    resS <- plsr(iForm, data = subData, ncomp = nComps, method = "oscorespls", scale = T) # scaled : POURQUOI ? "importance des wvls"
    
    # Coefficients (brut et standardisé)
    coefs <- as.vector(coef(resX, ncomp = nComps, intercept = T))
    zcoef <- as.vector(coef(resS, ncomp = nComps, intercept = T))
    
    coefMat[nsim,] <- coefs # à chaque ligne = no itération, ajoute les coeficients dans les colonnes respectives
    coefStd[nsim,] <- zcoef # standardized coeffis for ***importance of wvls*** /// SAME action
    
    # VIP
    vip <- c() # crée un objet "vip"
    for (j in seq(length(inBands))){ # 2000 bandes
      vip <- c(vip, VIPjh(resS, j, nComps)) # crée vip qui est un regroupement des vip précédantes, et du résultat VIPjh
      # où j = no de bande et h = le no de variables latentes choisies
    }
    
    vipMat[,nsim] <- vip # met "vip" dans vipMat[toutes lignes, colonne corr au no de simulation]
    
    # EXPLICATIONS ----
    # la formule ci-haut demande 
    # √ pas rapport que ça se fasse là; aurait cert. pu se faire n'importe quand avant lign 473... (?)
    # √ calcule les vip pour chaque bande dans chaque variable latente (dans chq boucle itérative)
    # √ crée vipMat (réutilisé ensuite, ligne 473 **)
    # ----
    
    # Statistiques du modèle
    fitX <- as.vector(unlist(resX$fitted.values[, 1, nComps])) # toutes lignes, 1iere colonne, nComps auquel on est rendu (nComps correspond au NO DE LISTE)
    preX <- as.vector(unlist(predict(resX, ncomp = nComps, newdata = tstData))) # prédit les résultats fittés à partir du set "testing"
    fitBias <- mean(as.numeric(subData[, inVar]) - fitX) # moyenne de l'erreur de chaque ligne (échantillon)
    valBias <- mean(as.numeric(tstData[, inVar]) - preX) # SAME, mais avec les valeurs prédites
    subData_df <- data.frame(subData[,])
    tstData_df <- data.frame(tstData[,])
    fitR2 <- summary(lm(subData_df[, inVar] ~ fitX))$r.squared
    valR2 <- summary(lm(tstData_df[, inVar] ~ preX))$r.squared
    fitRMSE <- sqrt(mean((as.numeric(subData[,inVar]) - fitX) ^ 2))
    valRMSE <- sqrt(mean((as.numeric(tstData[,inVar]) - preX) ^ 2))
    outVec <- c(fitR2, fitRMSE, fitBias, valR2, valRMSE, valBias)
    statMat[nsim,] <- outVec
  }
  
  statMat <- as.data.frame(statMat)
  names(statMat) <- c("fitR2","fitRMSE","fitBias","valR2","valRMSE","valBias")
  write.csv(statMat, paste("resultats/calcul_NRMSE/", inVar, "_", nComps,"comps_abs_stats.csv", sep=""), row.names=FALSE)
  write.csv(statMat, paste("resultats/", "800-2400_", inVar, "/", inVar, "_", nComps,"comps_abs_stats.csv", sep=""), row.names=FALSE)
  
  # EXPLICATIONS ----
  # la formule ci-haut demande 
  # √ met dans "statMat" les valeurs "outVec" calculées, mais à la ligne corr. au no d'intération
  # √ valeurs outVec : fitR2, fitRMSE, fitBias, valR2, valRMSE et valBias,
  # √ calculées à partir de resX et de prédictions subséquentes 
  # √ enregistre statMat
  # ----
  
  # Coefficients du modèle
  coeffis <- data.frame(matrix(nrow = length(inBands) + 1, ncol = 3))
  names(coeffis) <- c("bands", "mean","stdv")
  
  # entrer les données nécessaires dans le df créé
  coeffis$bands <- c("Intercept",inBands) # crée un vecteur où lignes = ça
  coeffis$mean <- apply(coefMat, MARGIN = 2, FUN = mean) # coefMat : coefficients de la plsr lignes 377, 378
  coeffis$stdv <- apply(coefMat, MARGIN = 2, FUN = sd)
  
  # enregistrer
  # write.table(coeffis, paste("~/Desktop/", inVar, "/", inVar, "_", nComps, "comps_abs_coeffMEAN.csv", sep = ""), 
  #             sep = ",", col.names = T, row.names = F)
  
  
  # Prédictions du modèle
  specMat <- aire_masse_df_ITER[, inBands] # ne sort que les colonnes inBands d'aire_masse_df_ITER
  specMat <- cbind(rep(1, nrow(specMat)), specMat) # colles-y le nombre de colonne = nb échantillons + le specMat existant (itérations)
  specMat <- as.matrix(specMat)
  
  predMat <- specMat %*% t(coefMat) # {base} x %*% y = matrix mutliplication // multiplie matrice specMat par transpose de coefMat (voir Legendre)
  predMean <- apply(predMat, FUN = mean, MARGIN = 1) # 
  predStdv <- apply(predMat, FUN = sd, MAR = 1)
  
  preds <- aire_masse_df_ITER[, !(names(aire_masse_df_ITER) %in% inBands | names(aire_masse_df_ITER) == "RAND")] # 2 critères de sélection des colonnes d'aire_masse_df_ITER à mettre dans
  # preds : [ {base} "x | y" = "OR" ] / DONC pas les noms d'aire_masse_df_ITER 
  # qui sont dans inBands OU METTRE les colonnes names(aire_masse_df_ITER) == "RAND")...
  # pourquoi référer à des lignes ? (RAND)
  preds[, paste("predMean_", inVar, sep = "")] <- predMean # ajoute colonne predMean_Var à df preds
  preds[, paste("predStdv_", inVar, sep = "")] <- predStdv # same, pr chaque Var
  write.csv(preds,paste("resultats/pred_meas/", inVar, "_", nComps,"comps_abs_preds.csv", sep=""), row.names=FALSE)
  write.csv(preds, paste("resultats/", "800-2400_", inVar, "/", inVar, "_", nComps,"comps_abs_preds.csv", sep=""), row.names=FALSE)
  
  # Visualisation des prédictions 
  modCI <- quantile(statMat$fitR2, probs = c(0.05, 0.95)) # CI : confidence intervals
  
  formi <- as.formula(paste(paste(inVar, " ~ predMean_", inVar, sep = ""))) # pourquoi paste(paste()) ?
  lmformi <-   as.formula(paste(paste("predMean_", inVar, " ~ ", inVar, sep = "")))
  
  jpeg(paste("resultats/", "800-2400_", inVar,  "/", inVar, "_", nComps, "comps_abs_predplot.jpg", sep = ""),
       width = 6, height = 5, units = "in", res = 300)
  plot(formi, data = preds, pch = 16, cex = 0.8, ylab = "measured", xlab = "predicted",
       main = inVar, xlim = c(min(predMean - predStdv), max(predMean + predStdv)))
  ###  abline(lm(lmformi, data = preds))
  abline(a = 0, b = 1, lty = 2)
  # arrows(predMean, aire_masse_df_ITER[, inVar], predMean + predStdv, aire_masse_df_ITER[, inVar], angle = 90, length = 0.05, lwd = 0.8)
  # arrows(predMean, aire_masse_df_ITER[, inVar], predMean - predStdv, aire_masse_df_ITER[, inVar], angle = 90, length = 0.05, lwd = 0.8)
  legend("topleft", bty = "n", cex = 0.8,
         c(paste("R² = ", sprintf("%.2f", signif(mean(statMat$fitR2), 3)), " [", signif(modCI[1], 2), ",", signif(modCI[2], 2), "]", sep = ""), 
           paste("RMSEP =", sprintf("%.2f", signif(mean(statMat$fitRMSE), 3)), sep = " "),
           paste("ncomps =", nComps, sep = " ")))
  dev.off()
  
  # Visualisation des VIP ("vip aggrégé")
  vipAggr <- as.data.frame(t(apply(vipMat, MARGIN = 1, FUN = quantile, probs = c(0.05,0.5,0.95)))) # crée df de la transpose (dcast genre) des quantiles de la vipMat
  # où vipMat = vip pour chaque bande dans chaque variable latente (dans chq boucle itérative)
  
  vipAggr$mean_VIP <- apply(vipMat, MARGIN = 1, FUN = mean) # ajoute à vipAggr la moyenne des VIP
  vipAggr$stdv <- apply(vipMat, MARGIN = 1, FUN = sd) # same, std
  stderr <- function(x) sqrt(var(x, na.rm = TRUE)/length(na.omit(x)))  # standard error of mean 
  vipAggr$std_err <- apply(vipMat, MARGIN = 1, FUN = stderr)
  vipAggr$band <- inBands # ajoute un colonne où ligne = identifiant de la bande correspondante
  # vipAggr utilisé dans visualisation, ligne 512
  vipAggr$band <- gsub("X", "", vipAggr$band)
  vipAggr$band <- as.numeric(vipAggr$band)
  write.csv(vipAggr,  file = paste("resultats/vip/", inVar, "_", "vipAggr.csv", sep = ""), row.names = FALSE)
  
  # Coefficients standardisés pour la visualisation 
  coeff_std <- data.frame(matrix(nrow = length(inBands) + 1, ncol = 3))
  names(coeff_std) <- c("bands", "mean", "stdv")
  
  # ajoute au df coeff_std les valeurs suivantes, tirées de la plsr
  coeff_std$bands <- c("Intercept", inBands) # id bande
  coeff_std$mean <- apply(coefStd, MARGIN = 2, FUN = mean) # coefStd créé ligne 385
  coeff_std$stdv <- apply(coefStd, MARGIN = 2, FUN = sd) # same
  
  # # Visualisation VIP et des coefficients standardisés 
  # jpeg(paste("resultats/", "800-2400_", inVar,  "/", inVar, "_", nComps, "comps_abs_varimp.jpg", sep = ""),
  #      width = 6, height = 7, units = "in", res = 300)
  # par(mfrow = c(2,1), mar = c(1.5,4,2.5,1.5), oma = c(3,0,0,0))
  # plot(coeff_std$mean[-1] ~ as.numeric(substr(coeff_std$bands[-1], 2, nchar(coeff_std$bands[-1]))),
  #      type = "p", pch = 19, xlab = "", ylab = "coeff_stdmean", main = paste(inVar, nComps, "comps", sep = "_"),
  #      ylim = c(-max(abs(coeff_std$mean[-1])), max(abs(coeff_std$mean[-1]))), bty = "l")
  # abline(h = 0)
  # points(abs(coeff_std$mean)[-1] ~ as.numeric(substr(coeff_std$bands[-1], 2, nchar(coeff_std$bands[-1]))),
  #        xlab = "wvl", ylab = "coeff_stdmean", col = 2, pch = 16, cex = 0.8)
  # 
  # lines(abs(coeff_std$mean)[-1] ~ as.numeric(substr(coeff_std$band[-1], 2, nchar(coeff_std$band[-1]))), col = 2)
  # 
  # plot(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))), vipAggr$mean_VIP, type = "l",
  #      xlab = "wvl", ylab = "VIP", bty ="l")
  # polygon(x = c(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))),
  #               rev(as.numeric(substr(vipAggr$band, 2, nchar(vipAggr$band))))),
  #         y = c(vipAggr$mean_VIP + vipAggr$stdv * 1.96, rev(vipAggr$mean_VIP - vipAggr$stdv * 1.96)),  
  #         col =  adjustcolor("red", alpha.f = 0.2), border = NA)
  # mtext("wavelength(nm)", 1, outer = T, line = 1)
  # dev.off()
}


