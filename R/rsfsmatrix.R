#' @export
rsfsmatrix <- function(vcf) {

  #beginin of function
  #Create a matrix with just the columns that contain the sample data

  df <- vcf@gt[,-1]

  #Create a vector with the name of the samples from the column names of df

  samples <- colnames(df)

  #Iterate the matrix to transform each result into the number of alternative alleles for that sample

  for(row in 1:nrow(df)) {
    for(col in 1:ncol(df)) {
      df[row,col] <- as.numeric(stringr::str_count(stringr::str_split(df[row,col], ":", n=2, simplify = TRUE)[1], pattern = '1'))
    }
  }

  #Here we obtain the rsfs as a matrix

  rsfs <- matrix(as.integer(df),ncol = ncol(df))     #Convert the strings to integer into a new matrix
  colnames(rsfs) <- samples                          #Add the samples vector to the column names
  return(rsfs)                                       #The rsfs
}
