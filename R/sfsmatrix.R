#' Calculate the sfs matrix of frequencies
#'
#' @param vcf A vcf file
#'
#' @return The sfs matrix
#' @export
sfsmatrix <- function(vcf) {

  #Here we obtain the sfs as a matrix of frequencies from the rsfs matrix calculated from a vcf file
  #The row values are added together to know the amount of altered polymorphisms. That is converted to a contingency table
  #and then as a matrix

  sfs <- as.matrix(table(rowSums(rsfsmatrix(vcf), na.rm = TRUE), exclude = 0))

  #The column names are changed to add an f so they look like: fi for each frequency

  colnames(sfs) <- 'absolute frequency'
  for(i in 1:length(rownames(sfs))) {
    rownames(sfs)[i] <- paste0('f', rownames(sfs)[i])
  }

return(sfs)              #Returns the sfs matrix of frequencies
}
