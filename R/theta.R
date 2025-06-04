#This is the helper omega() to build the omega weight
#
omega <- function(i, n, omega_type = "S", folded = FALSE, delta = 0) {

  if (omega_type == "xi1" & folded == FALSE) {
    i[i != 1] <- 0
    return(i)
  }

  if (omega_type == "S" & folded == FALSE) {
    return(1/i)
  }

  if (omega_type == "pi" & folded == FALSE) {
    return(n - i)
  }

  if (omega_type == "H" & folded == FALSE) {
    return(i)
  }

  if (omega_type == "xi1" & folded == TRUE) {
    i <- as.matrix(i)
    i[i != 1] <- 0
    for (x in 1:nrow(i)) {
      i[x, i[x,] == 1] <- n[x]
    }
    return(i)
  }

  if (omega_type == "S" & folded == TRUE) {
    return(n/(i * (n - i) * (1 + delta)))
  }

  if (omega_type == "pi" & folded == TRUE) {
    return(n/(1 + delta))
  }
}



#This is the helper sum_omega() to sum a vector of omega weights
#
sum_omega <- function(n, omega_type = "S", folded = FALSE, delta = 0) {

  for (x in 1:length(n)) {

    ifelse(folded == FALSE, len <- n[x]-1, len <- n[x] %/% 2)
    vector <- c(1:len)
    ifelse(folded == TRUE, delta <- delta(as.matrix(vector), n[x]), delta <- 0)
    n[x] <- sum(omega(vector, n[x], omega_type = omega_type, folded = folded, delta = delta))
  }

  return(n)

}



delta <- function(i, n) {

  for (x in 1:nrow(i)) {
    i[x, i[x,] != (n[x] - i[x,])] <- 0
    i[x, i[x,] == (n[x] - i[x,])] <- 1
  }
  return(i)

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
  del <- 0
  if (folded == TRUE) {
    del <- delta(i, nx)
  }
  omega <- omega(i, nx, omega_type = omega_type, folded = folded, delta = del)
  sum_omega <- sum_omega(nx, omega_type = omega_type, folded = folded, delta = del)

  if (folded == FALSE) {
    theta <- colSums(i * omega / sum_omega, na.rm = TRUE) #/ L
  }
  else {
    sample_names <- colnames(i)
    i <- mapply(function(x, n) {
      if (length(x) == 0 || is.na(n)) return(x)
      idx <- which(!is.na(x) & x > (n - x))
      x[idx] <- n - x[idx]
      x
    }, split(i, row(i)), nx, SIMPLIFY = FALSE)
    i <- do.call(rbind, i)
    colnames(i) <- sample_names

    phi <- nx/(i * (nx - i) * (1 + del))
    theta <- colSums((omega /phi) / sum_omega, na.rm = TRUE) #/ L
  }

  return(theta)
}





#' Compute theta for subsamples of genetic sequences using the relative Site Frequency Spectrum (rSFS).
#'
#' @param vcf A vcf file.
#' @param omega_type The model to follow for omega ("S", "pi", "xi1", "H"). "S" is set as default.
#' @param folded A logical variable. If true the rsfs is considered as folded (no knowledge over which is the alternative variant).
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
  del <- 0
  if (folded == TRUE) {
    del <- delta(i, nx)
  }
  omega <- omega(i, nx, omega_type = omega_type, folded = folded, delta = del)
  sum_omega <- sum_omega(nx, omega_type = omega_type, folded = folded, delta = del)
  psi <- 1 - (factorial(nx - i) * factorial(nx - j)) / (factorial(nx) * factorial(nx - i - j))

  if (folded == FALSE) {
    theta <- colSums((i * omega / psi) / sum_omega, na.rm = TRUE)
  }
  else {
    sample_names <- colnames(i)
    i <- mapply(function(x, n) {
      if (length(x) == 0 || is.na(n)) return(x)
      idx <- which(!is.na(x) & x > (n - x))
      x[idx] <- n - x[idx]
      x
      }, split(i, row(i)), nx, SIMPLIFY = FALSE)
    i <- do.call(rbind, i)
    colnames(i) <- sample_names

    phi <- nx/(i * (nx - i) * (1 + del))
    theta <- colSums((omega / (psi * phi)) / sum_omega, na.rm = TRUE)
  }


  if (subsample == "all") {
    return(theta)
  }
  else {
    return(theta[subsample])
  }

}

