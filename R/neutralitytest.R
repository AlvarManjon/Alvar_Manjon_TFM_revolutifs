#' Computes the difference between the theta of a subsample and the theta of another or the whole sample.
#'
#' @param vcf A vcf file.
#' @param s1 The name of the first sample to compare. A string object. It is necessary.
#' @param s2 The name of the second sample to compare. A string object. NULL by default and the whole theta of the whole sample is compared.
#' @param omega_type The model to follow for omega ("S", "pi", "xi1", "H"). "S" is set as default.
#' @param folded A logical variable. If true the sfs is considered as folded (no knowledge over which is the alternative variant).
#'
#' @returns The numerical difference between the thetas of both samples inserted.
#' @export
neutralitytest <- function(vcf, s1 , s2 = NULL, omega_type = "S", folded = FALSE) {

  theta1 <- rtheta(vcf, omega_type = omega_type, folded = folded, subsample = s1)

  if (is.null(s2)) {
    theta2 <- theta(vcf, omega_type = omega_type, folded = folded)
  }
  else {
    theta2 <- rtheta(vcf, omega_type = omega_type, folded = folded, subsample = s2)
  }

  result <- theta1 - theta2
  result <- as.numeric(result[1])

  return(result)
}
