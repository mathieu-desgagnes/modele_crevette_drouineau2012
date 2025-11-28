#' graphiques pour une zone (selectionner la zone avant d'appeler la fonction)
#'
#' @param val
#'
#' @returns
#' @export
#'
#' @examples
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
