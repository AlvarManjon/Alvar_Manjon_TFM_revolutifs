#' Calculate the SFS matrix of frequencies or its contingency table
#'
#' @param vcf A vcf file
#' @param cont_tab A logical variable. If TRUE the function returns the contingency table of SFS as a matrix. If FALSE it returns the SFS matrix.
#'
#' @return The SFS matrix or its contingency table
#' @export
sfsmatrix <- function(vcf, cont_tab = FALSE) {

  if (!is.logical(cont_tab) || length(cont_tab) != 1) {                      #Check the con_tab argument is just one logical value
    stop("Argument 'cont_tab' has to be a logical value (TRUE or FALSE)")    #Define the stop message if it is not logical
  }

  #Here we obtain the sfs as a matrix of frequencies from the genotype matrix calculated from a vcf file
  #We use the helper genotypematrix() over vcf
  #The row values are added together to know the amount of altered polymorphisms. That is converted to a contingency table
  #and then as a matrix

  sfs <- rowSums(genotypematrix(vcf), na.rm = TRUE)    #First create a vector of the amount of alternative alleles for each locus
  freqfac <- c(1:max(notmissing_vector(vcf)))          #Create a vector of all the possible altered polymorphism frequencies
                                                       #using the maximum of the helper notmissing() over the vcf file
  sfs_wfactors <- factor(sfs, levels = freqfac)        #We use the values of the last vector as factors for the alternative allele vector
  xi <- as.matrix(table(sfs_wfactors, exclude = 0))    #Create contingency table and transform it to a matrix

  #The column names are changed to add an f so they look like: fi for each frequency

  colnames(xi) <- 'absolute frequency'
  for(i in 1:length(rownames(xi))) {
    rownames(xi)[i] <- paste0('f', rownames(xi)[i])
  }

  #Finally, depending on the argument cont_tab, we return sfs as a matrix or the contingency table xi

  if (cont_tab == FALSE) {
    return(as.matrix(sfs))         #Returns the sfs as a matrix
  }
  else {
    return(xi)                     #Returns the contingency table of the sfs spectrum as a matrix
  }

}
