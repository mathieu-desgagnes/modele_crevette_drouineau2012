#' Formater les paramètre à ajuster selon les données considérées
#'
#' @param data une `list()` des données à ajuster par le modèle, notamment les années ajustées
#'
#' @returns une `list()` des paramètre qui seront optimisés
#' @export
#'
#' @examples
calculerParam <- function(data){
  param <- list()
  param$log_valLinf <- log(26)
  param$log_valK <- log(0.38)
  param$trans_mu1 <- rep(qlogis((10-4)/16), length(data$anneesFittees)+data$nbAge-1); names(param$trans_mu1) <- (min(data$anneesFittees)-data$nbAge+1):max(data$anneesFittees)
  param$log_tailleCV <- log(0.1)
  param$log_valRsex <- log(10)
  param$log_l50sex <- rep(log(20), length(data$anneesFittees)); names(param$log_l50sex) <- data$anneesFittees
  param$log_M <- log(0.5)
  param$log_valRselComm <- log(10)
  param$log_vall50selComm <- log(20)
  param$log_qComm <- rep(log(0.00001), nrow(data$nomPasDeTemps)); names(param$log_qComm) <- data$nomPasDeTemps$nom
  param$log_valTarget <- log(1.5)
  param$log_recrutement <- rep(log(1e8), length(data$anneesFittees)); names(param$log_recrutement) <- data$anneesFittees
  param$log_Nmale0 <- rep(log(1e7), data$nbAge-1)
  param$log_Nprimi0 <- log(1e7)
  param$log_Nmulti0 <- log(1e7)
  param
}
