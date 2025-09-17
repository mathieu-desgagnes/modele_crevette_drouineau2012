## ####
##
## Application de l'approche de Hilaire Drouineau 2012 au stocks de crevette
## Par Mathieu Desgagnés, printemps 2025
##
## ####

require('RTMB')

sourceGenerale <- file.path('//dcqcimlna01a','Projets','AP_INV','modele_crevette_drouineau2012')
anneesFittees <- 1990:2023

## charger les données
calculerData <- function(annees=1990:2023, new_strata=0){
    ## lecture des fichiers d'observation
    fl_abond <- read.csv2(file=file.path(sourceGenerale, 'data', 'fl_zone_sexe2_abd.csv'), dec='.')
    ##
    data <- list()
    data$anneesFittees <- annees
    data$nbAge <- 4 #age max de mâles, changement de sexe obligatoire ensuite
    data$midTaille <- unique(fl_abond$lc)
    data$nbMoisParSaison <- c(2,3,3,4)            #printemps=c(avril,mai), ete=c(juin:aout), automne=c(sept:nov), hiver=c(dec:mars)
    ##
    data$fl_abond_male <- subset(fl_abond, sexe2=='ma' & annee%in%data$anneesFittees & newstrata==new_strata, select=c('annee','zone','lc','abd_lc','newstrata'))
    ##
    ## 4 saison ("pas de temp") débutant au printemps
    data$nbPasDeTemps <- length(data$anneesFittees)*4
    data$effort <- rep(NA, data$nbPasDeTemps)
    ##
    data
}
data <- calculerData(annees=1990:2023)

calculerParam <- function(data){
    param <- list()
    param$log_valLinf <- log(26)
    param$log_valK <- log(0.38)
    param$trans_mu1 <- rep(qlogis((10-4)/16), length(data$anneesFittees)+data$nbAge)
    param$log_tailleCV <- log(0.1)
    param$log_valRsex <- log(10)
    param$log_l50sex <- rep(log(20), length(data$anneesFittees))
    param$log_M <- log(0.5)
    param$log_valRselComm <- log(10)
    param$log_vall50selComm <- log(20)
    param$log_qComm <- rep(log(0.1), length(data$anneesFittees)*4)
    param$log_valTarget <- log(1.5)
    param$log_recrutement <- rep(log(1e8), length(data$anneesFittees))
    param$log_Nmale0 <- rep(log(1e7), data$nbAge-1)
    param$log_Nprimi0 <- log(1e7)
    param$log_Nmulti0 <- log(1e7)
    param
}
param <- calculerParam(data=data)

## graphiques pour une zone (selectionner la zone avant d'appeler la fonction)
struct_taille <- function(val){
    par(mfrow=c(4,8))
    for(j in unique(val$annee)){
        temp <- subset(val, annee==j)
        ##
        plot(temp$lc, temp$abd_lc,
             type = "n",  # N'affiche rien pour l'instant
             xlim = range(temp$lc),
             ylim = range(temp$abd_lc),
             xlab = "Longueur de carapace (LC)",
             ylab = "Abondance",
             main = "Abondance par longueur de carapace")
        mid_classe <- unique(temp$lc)
        lim_classe <- c(0, tail(mid_classe,-1)-diff(mid_classe)/2, tail(mid_classe,1)+tail(diff(mid_classe)/2,1))
        lines(c(0,0), c(0,temp[temp$lc==mid_classe[1],'abd_lc']))
        for(lc in head(seq_along(mid_classe),-1)){
            lines(lim_classe[lc+c(0,1)], rep(temp[temp$lc==mid_classe[lc],'abd_lc'],2))
            lines(rep(lim_classe[lc+1],2), temp[temp$lc%in%mid_classe[lc+c(0,1)],'abd_lc'])
        }
    }
}

## modèle d'analyse modale
fnll <- function(param, fit=TRUE){
    getAll(param, data)
    nbAn <- length(anneesFittees)

    ## croissance
    valLinf <- exp(log_valLinf)
    valK <- exp(log_valK)
    mu1 <- 4 + 16*plogis(trans_mu1)
    tailleCV <- exp(log_tailleCV)

    ## changement de sexe
    valRsex <- exp(log_valRsex)
    l50sex <- exp(log_l50sex) #marche aléatoire

    ## mortalité naturelle et survie
    M <- exp(log_M) #marche aléatoire

    ## sélectivité commercial (male uniquement)
    valRselComm <- exp(log_valRselComm)
    vall50selComm <- exp(log_vall50selComm)
    selCommMale <- 1 / (1 + exp(-2*log(3)/valRselComm * (1:length(midTaille) - vall50selComm)))

    ## mortalité par la pêche
    ##

    ##
    ## attention: ici, la capturabilité commerciale est par pas de temps, et non année
    ##

    qComm <- exp(log_qComm) #marche aléatoire
    Fmale <- array(NA, dim=c(length(midTaille), nbPasDeTemps))
    for(i.long in 1:length(midTaille)){
        Fmale[i.long,] <- qComm * selCommMale[i.long] * effort
    }
    ## Fprimi <- qComm * effort
    ## valTarget <- exp(log_valTarget)
    ## for(i.step in 1:(nbPasDeTemps/4)){
    ##     Fmulti[4*(i.step-1)+1] <- qComm[i.step] * valTarget * effort[i.step]
    ##     Fmulti[4*(i.step-1)+2] <- qComm[i.step] * effort[i.step]
    ##     Fmulti[4*(i.step-1)+3] <- qComm[i.step] * effort[i.step]
    ##     Fmulti[4*(i.step-1)+4] <- qComm[i.step] * effort[i.step]
    ## }

    ## survie
    ##
    ## attention: dans la survie, considérer le nombre de mois réel

    SrMale <- array(NA, dim=c(length(midTaille), nbPasDeTemps))
    for(i.step in 1:nbPasDeTemps){
        SrMale[,i.step] <- exp(-(M + Fmale[,i.step])*1/4)
    }
    ## SrPrimi <- exp(-(M + Fprimi)*1/4)
    ## SrMulti <- exp(-(M + FMulti)*1/4)

    ## changement de sex
    propCsex <- array(0, dim=c(length(midTaille), nbAn)
    for(i.long in 1:length(midTaille)){
        propCsex[i.long,] <- 1 / (1 + exp(-2*log(3)/valRsex * (i.long - l50sex)))
    }


    ## longueur à l'age, pour chaque cohorte et pas de temps
    laa <- array(NA, dim=c(nbAn+nbAge, nbAge*4),
                 dimnames=list(cohorte=c(min(data$anneesFittees)+c(-nbAge:-1),data$anneesFittees), pasDeTemp=sprintf("%d%s", rep(1:nbAge, each=4), rep(c('A','B','C','D'), nbAge))))
    laa[,1] <- mu1
    for(i.step in 2:(nbAge*4)){
        laa[,i.step] <- laa[,i.step-1] + (valLinf-laa[,1]) * (1-exp(-valK*(nbMoisParSaison[(i.step-2)%%4+1]/12)))
    }

    NmaleAge <- array(NA, dim=c(nbAge, nbPasDeTemps), dimnames=list(age=1:nbAge, temps=1:nbPasDeTemps)) #au début du pas de temps
    NmaleAgeTaille <- array(NA, dim=c(nbAge, nbPasDeTemps, length(midTaille)), dimnames=list(age=1:nbAge, temps=1:nbPasDeTemps, taille=midTaille)) #une fois la croissance et la mortalité appliquée
    ## Nprimi <- rep(NA, nbPasDeTemps)
    ## Nmulti <- rep(NA, nbPasDeTemps)
    ## recrutement (au printemps)
    NmaleAge[1,1] <- exp(log_recrutement[1])
    for(i.step in 1){
        NmaleAgeTaille[1,i.step,] <- (pnorm(midTaille+0.25, mean=laa[1,i.step], sd=tailleCV*laa[1,i.step]) -
                                      pnorm(midTaille-0.25, mean=laa[1,i.step], sd=tailleCV*laa[1,i.step])) *
            NmaleAge[1,i.step]
    }
    ## année 0 (pas de temps 1, printemps)
    NmaleAge[2:nbAge,1] <- exp(log_Nmale0)
    for(i.age in 2:nbAge){ #il y a des cohorte avant la première année pour tenir compte des age >1 la première année
        NmaleAgeTaille[i.age,1,] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
                                     pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
            NmaleAge[i.age,1]
    }
    ## Nprimi[1] <- exp(log_Nprimi0)
    ## Nmulti[1] <- exp(log_Nmulti0)

    ## ensuite pour chaque année
    for(i.an in 1:(nbPasDeTemps/4)){
        ## saison 2 (été)
        i.step <- (i.an-1)*4+2
        for(i.age in 1:nbAge){ #
            ## croissance et redistribution des individus en classe de taille
            NmaleAgeTaille[i.age,i.step,] <- (pnorm(midTaille+0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+2], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+2]) -
                                              pnorm(midTaille-0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+2], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+2])) *
                NmaleAge[i.age,i.step-1]
            ##
            ## survie à la longueur (sur crevette après croissance)
            NmaleAgeTaille[i.age,i.step,] <- NmaleAgeTaille[i.age,i.step-1,] * SrMale[,i.step-1]
            ##
            ## survie à l'age
            NmaleAge[i.age,i.step] <- sum(NmaleAgeTaille[i.age,i.step,])
            ##
        }
        ## Nprimi[i.step] <- Nprimi[i.step-1] * SrPrimi[i.step-1]
        ## Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        ##
        ## saison 3 (automne, primi deviennent multi)
        i.step <- (i.an-1)*4+3
        for(i.age in 1:nbAge){ #
            ## croissance et redistribution des individus en classe de taille
            NmaleAgeTaille[i.age,i.step,] <- (pnorm(midTaille+0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+3], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+3]) -
                                              pnorm(midTaille-0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+3], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+3])) *
                NmaleAge[i.age,i.step-1]
            ##
            ## survie à la longueur (sur crevette après croissance)
            NmaleAgeTaille[i.age,i.step,] <- NmaleAgeTaille[i.age,i.step-1,] * SrMale[,i.step-1]
            ##
            ## survie à l'age
            NmaleAge[i.age,i.step] <- sum(NmaleAgeTaille[i.age,i.step,])
            ##
        }
        ## Nprimi[i.step] <- 0
        ## Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1] + Nprimi[i.step-1] * SrPrimi[i.step-1]
        ##
        ## saison 4 (hiver, aucune primi)
        i.step <- (i.an-1)*4+4
        for(i.age in 2:(nbAge)){ #
            ## croissance et redistribution des individus en classe de taille
            NmaleAgeTaille[i.age,i.step,] <- (pnorm(midTaille+0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+4], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+4]) -
                                              pnorm(midTaille-0.25, mean=laa[i.an+(4-i.age),(i.age-1)*4+4], sd=tailleCV*laa[i.an+(4-i.age),(i.age-1)*4+4])) *
                NmaleAge[i.age,i.step-1]
            ##
            ## survie à la longueur (sur crevette après croissance)
            NmaleAgeTaille[i.age,i.step,] <- NmaleAgeTaille[i.age,i.step-1,] * SrMale[,i.step-1]
            ##
            ## survie à l'age
            NmaleAge[i.age,i.step] <- sum(NmaleAgeTaille[i.age,i.step,])
            ##
        }
        Nprimi[i.step] <- 0
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        ##
        ## saison 1 (printemps: males avancent en age et changent de sexe, nouvelle primi)
        i.step <- (i.an-1)*4+5
        for(i.age in 1:nbAge){ #
            ## croissance et redistribution des individus en classe de taille
            NmaleAgeTaille[i.age,i.step,] <- (pnorm(midTaille+0.25, mean=laa[i.an+(4-i.age)+1,(i.age-1)*4+1], sd=tailleCV*laa[i.an+(4-i.age)+1,(i.age-1)*4+1]) -
                                              pnorm(midTaille-0.25, mean=laa[i.an+(4-i.age)+1,(i.age-1)*4+1], sd=tailleCV*laa[i.an+(4-i.age)+1,(i.age-1)*4+1])) *
                NmaleAge[i.age,i.step-1]
            ##

            ## NmaleAge[i.age,i.step] <- NmaleAge[i.age-1,i.step-1] * SrMale[i.step-1]
            ## NmaleAgeTaille[i.age,i.step,] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
            ##                                   pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
            ##     NmaleAge[i.age,i.step]
        }
        Nprimi[i.step] <- 0

        ## attention ici, vérifier que la structure de taille est adéquate pour le pas de temps...
        ##

        for(i.age in 1:nbAge){ #
            nbChangeSex.temp <- propCsex[,i.an] * NmaleAgeTaille[i.age,i.step-1,]
            NmaleAgeTaille[i.age+1,i.step,] <- NmaleAgeTaille[i.age,i.step-1,] - nbChangeSex.temp



            NmaleAge[i.age,i.step] <- NmaleAge[i.age,i.step] - nbChangeSex.temp
            Nprimi[i.step] <- Nprimi[i.step] + nbChangeSex.temp
        }
        Nprimi[i.step] <- Nprimi[i.step] + NmaleAge[nbAge,i.step-1] * SrMale[i.step-1]
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        NmaleAge[1,i.step] <- exp(log_recrutement[i.an])
        ##
    }

    ## sélectivité relevé male
    valRselReleve <- exp(log_valRselReleve)
    vall50selReleve <- exp(log_vall50selReleve)
    selReleveMale <- 1 / (1 + exp(-2*log(3)/valRselReleve * (1:length(midTaille) - vall50selReleve)))

    ## likelihood sur les structure de taille relevé

    REPORT(nll)
    ##
    nll
}


## modèle
fnll <- function(param, fit=TRUE){
    getAll(param, data)

    ## croissance
    valLinf <- exp(log_valLinf)
    valK <- exp(log_valK)
    mu1 <- 0.04 + 0.16*plogis(trans_mu1)
    tailleCV <- exp(log_tailleCV)

    ## changement de sexe
    valRsex <- exp(log_valRsex)
    l50sex <- exp(log_l50sex) #marche aléatoire

    ## mortalité naturelle et survie
    M <- exp(log_M) #marche aléatoire

    ## sélectivité commercial (male uniquement)
    valRselComm <- exp(log_valRselComm)
    vall50selComm <- exp(log_vall50selComm)
    selCommMale <- 1 / (1 + exp(-2*log(3)/valRselComm * (1:length(midTaille) - vall50selComm)))

    ## mortalité par la pêche
    qComm <- exp(log_qComm) #marche aléatoire
    for(i.long in 1:length(midTaille)){
        Fmale[i.long,] <- qComm * selMale[i.long] * effort
    }
    Fprimi <- qComm * effort
    valTarget <- exp(log_valTarget)
    for(i.step in 1:(nbPasDeTemps/4)){
        Fmulti[4*(i.step-1)+1] <- qComm[i.step] * valTarget * effort[i.step]
        Fmulti[4*(i.step-1)+2] <- qComm[i.step] * effort[i.step]
        Fmulti[4*(i.step-1)+3] <- qComm[i.step] * effort[i.step]
        Fmulti[4*(i.step-1)+4] <- qComm[i.step] * effort[i.step]
    }

    ## survie
    for(i.step in 1:nbPasDeTemps){
        SrMale[i.step] <- exp(-(M + Fmale[,i.step])*1/4)
    }
    SrPrimi <- exp(-(M + Fprimi)*1/4)
    SrMulti <- exp(-(M + FMulti)*1/4)

    laa[,1] <- mu1
    for(i.step in 2:nbPasDeTemps){
        laa[,i.step] <- laa[,i.step-1] + (valLinf-laa[,1]) * (1-exp(-valK*(1/4)))
    }

    propCsex <- array(0, dim=c(length(midTaille), nbPasDeTemps))
    for(i.long in 1:length(midTaille)){
        propCsex[i.long,] <- 1 / (1 + exp(-2*log(3)/valRsex * (i.long - l50sex)))
    }


    NmaleAge <- array(0, dim=c(nbAge, nbPasDeTemps))
    NmaleAgeTaille <- array(0, dim=c(nbAge, length(midTaille), nbPasDeTemps))
    Nprimi <- rep(0, nbPasDeTemps)
    Nmulti <- rep(0, nbPasDeTemps)
    ## recrutement
    NmaleAge[1,] <- exp(log_recrutement)
    for(i.step in 1:nbPasDeTemps){
        NmaleAgeTaille[1,,i.step] <- (pnorm(midTaille+0.25, mean=laa[1,i.step], sd=tailleCV*laa[1,i.step]) -
                                      pnorm(midTaille-0.25, mean=laa[1,i.step], sd=tailleCV*laa[1,i.step])) *
            NmaleAge[1,i.step]
    }
    ## année 0 (pas de temps 1, printemps)
    NmaleAge[2:nbAge,1] <- exp(log_Nmale0)
    for(i.age in 2:nbAge){
        NmaleAgeTaille[i.age,,1] <- (pnorm(midTaille+0.25, mean=laa[i.age,1], sd=tailleCV*laa[i.age,1]) -
                                     pnorm(midTaille-0.25, mean=laa[i.age,1], sd=tailleCV*laa[i.age,1])) *
            NmaleAge[i.age,1]
    }
    Nprimi[1] <- exp(log_Nprimi0)
    Nmulti[1] <- exp(log_Nmulti0)
    ##
    for(i.an in 1:(nbPasDeTemps/4)){
        ## saison 2 (été)
        i.step <- (i.an-1)*4+2
        for(i.age in 2:nbAge){ #
            NmaleAge[i.age,i.step] <- NmaleAge[i.age,i.step-1] * SrMale[i.step-1]
            NmaleAgeTaille[i.age,,i.step] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
                                              pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
                NmaleAge[i.age,i.step]
        }
        Nprimi[i.step] <- Nprimi[i.step-1] * SrPrimi[i.step-1]
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        ##
        ## saison 3 (automne, primi deviennent multi)
        i.step <- (i.an-1)*4+3
        for(i.age in 2:(nbAge)){ #
            NmaleAge[i.age,i.step] <- NmaleAge[i.age,i.step-1] * SrMale[i.step-1]
            NmaleAgeTaille[i.age,,i.step] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
                                              pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
                NmaleAge[i.age,i.step]
        }
        Nprimi[i.step] <- 0
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1] + Nprimi[i.step-1] * SrPrimi[i.step-1]
        ##
        ## saison 4 (hiver, aucune primi)
        i.step <- (i.an-1)*4+4
        for(i.age in 2:(nbAge)){ #
            NmaleAge[i.age,i.step] <- NmaleAge[i.age,i.step-1] * SrMale[i.step-1]
            NmaleAgeTaille[i.age,,i.step] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
                                              pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
                NmaleAge[i.age,i.step]
        }
        Nprimi[i.step] <- 0
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        ##
        ## saison 1 (printemps: males avancent en age et changent de sexe, nouvelle primi)
        i.step <- (i.an-1)*4+5
        for(i.age in 2:nbAge){ #
            NmaleAge[i.age,i.step] <- NmaleAge[i.age-1,i.step-1] * SrMale[i.step-1]
            NmaleAgeTaille[i.age,,i.step] <- (pnorm(midTaille+0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1]) -
                                              pnorm(midTaille-0.25, mean=laa[1+(4-i.age),(i.age-1)*4+1], sd=tailleCV*laa[1+(4-i.age),(i.age-1)*4+1])) *
                NmaleAge[i.age,i.step]
        }
        Nprimi[i.step] <- 0
        for(i.age in 2:nbAge){ #
            nbChangeSex.temp <- propCsex[,i.step-1] * NmaleAgeTaille[i.age-1,,i.step-1]
            NmaleAge[i.age,i.step] <- NmaleAge[i.age,i.step] - nbChangeSex.temp
            Nprimi[i.step] <- Nprimi[i.step] + nbChangeSex.temp
        }
        Nprimi[i.step] <- Nprimi[i.step] + NmaleAge[nbAge,i.step-1] * SrMale[i.step-1]
        Nmulti[i.step] <- Nmulti[i.step-1] * SrMulti[i.step-1]
        ##
    }

    ## captures (Baranov)
    CmaleAgeTaille <- array(0, dim=c(nbAge, length(midTaille), nbPasDeTemps))
    CmaleAge <- array(0, dim=c(nbAge, nbPasDeTemps))
    Cmale <- rep(0, nbPasDeTemps)
    Cprimi <- rep(0, nbPasDeTemps)
    Cmulti <- rep(0, nbPasDeTemps)
    for(i.step in 1:nbPasDeTemps){
        for(i.age in 1:nbAge){
            CmaleAgeTaille[i.age,,i.step] <- NmaleAgeTaille[i.age,,i.step] *
                Fmale[,i.step]/(M[i.step]+Fmale[,i.step]) * (1-exp(-(M[i.step]+Fmale[,i.step])*1/4))
            CmaleAge[i.age,i.step] <- sum(CmaleAgeTaille[i.age,,i.step])
        }
        Cmale[i.step] <- sum(CmaleAge[,i.step])
        Cprimi[i.step] <- Nprimi[i.step] * Fprimi[i.step]/(M[i.step]+Fprimi[i.step]) * (1-exp(-(M[i.step]+Fprimi[i.step])*1/4))
        Cmulti[i.step] <- Nmulti[i.step] * Fmulti[i.step]/(M[i.step]+Fmulti[i.step]) * (1-exp(-(M[i.step]+Fmulti[i.step])*1/4))
    }

    ## taux d'exploitation... nécessaire?

    ## sélectivité relevé male
    valRselReleve <- exp(log_valRselReleve)
    vall50selReleve <- exp(log_vall50selReleve)
    selReleveMale <- 1 / (1 + exp(-2*log(3)/valRselReleve * (1:length(midTaille) - vall50selReleve)))

    ## indice abondance
    qRel <- exp(log_qRel)
    iaMale <- rep(0, nbPasDeTemps)
    iaMaleAge <- array(0, dim=c(nbAge, nbPasDeTemps))
    iaMaleAgeTaille <- array(0, dim=c(nbAge, length(midTaille), nbPasDeTemps))
    for(i.step in 1:nbPasDeTemps){
        for(i.age in 1:nbAge){
            iaMaleAgeTaille[i.age,,i.step] <- selReleveMale * qRel * NmaleAgeTaille[i.age,,i.step]
            iaMaleAge[i.age,i.step] <- sum(iaMaleAgeTaille[i.age,,i.step])
        }
        iaMale[i.step] <- sum(iaMaleAge[,i.step])
    }
    iaFem <- qRel * [Nprimi + Nmulti]

    ## likelihood sur les structure de taille commerciales et relevé
    ## likelihood sur total catch par classe
    ## likelihood sur indice d'abondance par classe

    REPORT(nll)
    ##
    nll
}

randomVal <- c('logRpred','logBpred')
## randomVal <- c('logRpred','logBpred','transTauxExp')
##
mapVal <- list(log_sigma_C=factor(NA), log_sigma_longAge=factor(NA))
obj <- MakeADFun(fnll, dd_param, random=randomVal, map=mapVal)

##
## fit <- nlminb(obj$par, obj$fn, obj$gr)
fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max=100000000, iter.max=100000000))
fit
##
sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")
plr <- as.list(sdr, "Est", report=TRUE)
plrsd <- as.list(sdr, "Std", report=TRUE)
