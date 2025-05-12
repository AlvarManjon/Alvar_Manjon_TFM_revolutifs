#This is the helper omega() to build the omega weight
#
omega <- function(i, n, omega_type = "S", folded = FALSE) {

  if (omega_type == "xi1") {
    i[i != 1] <- 0
    return(i)
  }

  if (omega_type == "S") {
    return(1/i)
  }

  if (omega_type == "pi") {
    return(n - i)
  }

  if (omega_type == "H") {
    return(i)
  }

}



#This is the helper sum_omega() to sum a vector of omega weights
#
sum_omega <- function(n, omega_type = "S", folded = FALSE) {

  for (x in 1:length(n)) {

    vector <- c(1:(n[x]-1))
    n[x] <- sum(omega(vector, n[x], omega_type = omega_type, folded = folded))
  }

  return(n)

}



#' Compute theta for a complete sample of genetic sequences using the Site Frequency Spectrum (SFS).
#'
#' @param vcf A vcf file.
#' @param omega_type The model to follow for omega ("S", "pi", "xi1", "H"). "S" is set as default.
#' @param folded A logical variable. If true the sfs is considered as folded (no knowledge over which is the alternative variant).
#'
#' @returns The numerical value of theta
#' @export
theta <- function(vcf, omega_type = "S", folded = FALSE) {

  i <- sfsmatrix(vcf)
  nx <- notmissing_vector(vcf)
  i[i == nx] <- NA
  i[i == 0] <- NA
  omega <- omega(i, nx, omega_type = omega_type, folded = folded)
  sum_omega <- sum_omega(nx, omega_type = omega_type, folded = folded)
  #L <- colSums(i != 0, na.rm = TRUE)


  theta <- colSums(i * omega / sum_omega, na.rm = TRUE) #/ L

  return(theta)
}





#' Compute theta for subsamples of genetic sequences using the relative Site Frequency Spectrum (rSFS).
#'
#' @param vcf A vcf file.
#' @param omega_type The model to follow for omega ("S", "pi", "xi1", "H"). "S" is set as default.
#' @param folded A logical variable. If true the sfs is considered as folded (no knowledge over which is the alternative variant).
#' @param subsample The name of the subsample to compute theta for ("example"). "all" set as default.
#'
#' @returns The numerical value of theta for a subsample or a vector with the theta values for each subsample.
#' @export
rtheta <- function(vcf, omega_type = "S", folded = FALSE, subsample = "all") {

  i <- rsfsmatrix(vcf)
  j <- notmissing_matrix(vcf)
  nx <- rowSums(j, na.rm = TRUE)
  i[i == nx] <- NA
  i[i == 0] <- NA
  omega <- omega(i, nx, omega_type = omega_type, folded = folded)
  sum_omega <- sum_omega(nx, omega_type = omega_type, folded = folded)
  psi <- (factorial(nx - i) * factorial(nx - j)) / (factorial(nx) * factorial(nx - i - j))

  theta <- colSums(i * omega * psi / sum_omega, na.rm = TRUE)

  if (subsample == "all") {
    return(theta)
  }
  else {
    return(theta[subsample])
  }

}

