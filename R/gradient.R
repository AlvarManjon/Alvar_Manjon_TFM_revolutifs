#' Plots the spatial distribution of theta in a plane
#'
#' @param theta a matrix with the theta values of all samples.
#' @param long a matrix or vector containing the longitude coordinates of the samples in the same order as theta
#' @param lat a matrix or vector containing the latitude coordinates of the samples in the same order as theta
#'
#' @returns the plot of theta for every sample displayed with its coordinates
#' @importFrom rlang .data
#' @export
gradient <- function(theta, long = NULL, lat = NULL) {

  theta <- as.matrix(theta)
  ifelse(is.null(long), long <- stats::runif(nrow(theta), min = -10, max = 10), long <- long)
  ifelse(is.null(lat), lat <- stats::runif(nrow(theta), min = -10, max = 10), lat <- lat)

  df <- data.frame(
    Sample = rownames(theta),
    Latitude = long,
    Longitude = lat,
    Theta = theta
  )

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Longitude, y = .data$Latitude)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Theta, color = .data$Theta), alpha = 0.8) +
    ggplot2::scale_color_gradient(low = "white", high = "blue") +
    ggplot2::scale_size(range = c(2, 8)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Sample), vjust = -1.2, size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Theta Spatial Distribution",
      x = "Longitude",
      y = "Latitude",
      color = expression(theta),
      size = expression(theta)
    ) +
    ggplot2::theme(legend.position = "right")

  return(plot)

}
