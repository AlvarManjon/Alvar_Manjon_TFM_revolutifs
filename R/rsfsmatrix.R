#This is the helper genotypematrix() to build the genotype matrix from the vcf file.
#It returns the numeric matrix of the amount of alternative alleles per sample.
genotypematrix <- function(vcf) {

  #begining of the function
  #Create a matrix with just the columns that contain the sample genotype data

  df <- vcf@gt[,-1]

  #Create a vector with the name of the samples from the column names of df

  samples <- colnames(df)

  #Iterate the matrix to transform each result into the number of alternative alleles for that sample

  for(row in 1:nrow(df)) {
    for(col in 1:ncol(df)) {
      df[row,col] <- stringr::str_count(stringr::str_split(df[row,col], ":", n=2, simplify = TRUE)[1], pattern = '1')
    }
  }

  #Here we obtain the matrix

  genomatrix <- matrix(as.integer(df),ncol = ncol(df))     #Convert the strings to integer into a new matrix
  colnames(genomatrix) <- samples                          #Add the samples vector to the column names

  return(genomatrix)                                       #The matrix is returned

}



#This is the helper notmissing_matrix() to build the sample size of non missing data from the vcf file.
#It returns the numeric matrix of the sample size of non missing data for each locus and for each subsample.
notmissing_matrix <- function(vcf) {

  #begining of the function
  #Create a matrix with just the columns that contain the sample genotype data

  df <- vcf@gt[,-1]

  #Iterate the matrix to transform each result into the number of not missing alleles for that sample (we take into account at least 9 alternative alleles)

  for(row in 1:nrow(df)) {
    for(col in 1:ncol(df)) {
      df[row,col] <- stringr::str_count(stringr::str_split(df[row,col], ":", n=2, simplify = TRUE)[1], pattern = '(0|1|2|3|4|5|6|7|8|9)')
    }
  }

  #Here we obtain the vector

  notmissing <- matrix(as.integer(df),ncol = ncol(df))     #Convert the strings to integer into a new matrix

  return(notmissing)                                       #The matrix is returned

}



#This is the helper notmissing_vector() to build the sample size of non missing data from the vcf file.
#It returns the numeric vector of the sample size of non missing data for each locus.
#This helper is a complement for the notmissing_matrix() helper
notmissing_vector <- function(vcf) {

  #We start from the matrix returned by the helper notmissing_matrix() applied over the vcf file

  nx <- rowSums(notmissing_matrix(vcf), na.rm = TRUE)      #Add every value of each row into a vector

  #This vector includes the sample size of the not missing data for each locus

  return(nx)                                               #The vector is returned

}



#' Calculate the rsfs matrix for every sample and every polymorphic site
#'
#' @param vcf A vcf file
#'
#' @return The rsfs matrix
#' @export
rsfsmatrix <- function(vcf) {

  #First we create the matrix of alternative alleles with the genotypematrix() helper

  rsfs <- genotypematrix(vcf)

  #We substitute every cell which has at least one alternative allele with the sum of the row
  #For this we iterate the matrix once again

  for(row in 1:nrow(rsfs)) {
    sumrow <- sum(rsfs[row,], na.rm = TRUE)          #Save the sum of the row before iterating the row
    for(col in 1:ncol(rsfs)) {
      if (sum(rsfs[row,col], na.rm = TRUE) > 0) {    #If the cell is not 0 or missing value
        rsfs[row,col] <- sumrow                      #Transform the cell into the sum of the row
      }
    }
  }

  return(rsfs)                                       #Return the rsfs matrix

}
