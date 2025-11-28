#' Préparer les données et pramètres invariables utilisés dans le modèle
#'
#' @param annees la séquence d'années à ajuster par le modèle
#' @param new_strata
#'
#' @returns une `list()` des données et paramètres invariables
#' @export
#'
#' @examples
calculerData <- function(annees=1990:2007, new_strata=0){
  ## lecture des fichiers d'observation
  ##
  data <- list()
  data$anneesFittees <- annees; names(data$anneesFittees) <- data$anneesFittees
  data$nbAge <- 4 #age max de mâles, changement de sexe obligatoire ensuite
  data$saisons <- as.data.frame(list(nom=c('a','b','c','d'), nbMois=c(2,3,4,3))); dimnames(data$saisons)[[1]] <- data$saisons$nom
  data$midTaille <- seq(0, 29.5, by=0.5)+0.25; names(data$midTaille) <- data$midTaille #classe de tailles de 0.5cm
  ##
  ##
  ## 4 saison ("pas de temp") débutant au printemps
  data$nbPasDeTemps <- length(data$anneesFittees)*4
  temp <- expand.grid(annee=anneesFittees,saison=data$saisons$nom)
  temp <- temp[order(temp$annee, temp$saison),]
  temp$nom <- paste0(temp$annee, temp$saison)
  data$nomPasDeTemps <- temp
  ## effort commercial
  effort.init <- read.csv2(file=file.path('data', origineDonnees, 'effort_comm.csv'))
  effort.init$val <- effort.init$effort; effort.init$effort <- NULL
  effort.init$saison <- c('a','b','c','d')[effort.init$saison]
  data$effort_comm <- effort.init
  ##
  fl_abond.init <- read.csv2(file=file.path('data', origineDonnees, 'dfl_rel.csv'))
  fl_abond.male <- fl_abond.init[,1:2]; fl_abond.male$mids <- as.numeric(substr(dimnames(fl_abond.init)[[2]][3],2,10)); fl_abond.male$val <- fl_abond.init[,3]
  for(i.cl in 4:(ncol(fl_abond.init)-3)){
    temp <- fl_abond.init[,1:2]
    temp$mids <- as.numeric(substr(dimnames(fl_abond.init)[[2]][i.cl],2,10))
    temp$val <- fl_abond.init[,i.cl]
    fl_abond.male <- rbind(fl_abond.male, temp)
  }
  data$fl_rel.male <- fl_abond.male
  ##
  fl_abond.femelle <- fl_abond.init[,1:2]; fl_abond.femelle$stade <- 'femelle'; fl_abond.femelle$val <- fl_abond.init[,'femelles']
  for(i.cl in c('primi','multi')){
    temp <- fl_abond.init[,1:2]
    temp$stade <- i.cl
    temp$val <- fl_abond.init[,i.cl]
    fl_abond.femelle <- rbind(fl_abond.femelle, temp)
  }
  data$fl_rel.femelle <- fl_abond.femelle
  ##
  data
}

