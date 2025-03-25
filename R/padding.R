###############################################################################
## Enhanced padding methods for wavelet transforms
## Added to extend the original WaveletComp package functionality
###############################################################################

#' Padding methods for wavelet transforms
#'
#' @param x Numeric vector, the time series to be padded
#' @param pad_method String, method for padding: "zero" (original method),
#'        "reflection" (MATLAB-style reflection), or "reflection_zero"
#'        (reflection followed by zero padding to next power of 2)
#' @return A list containing the padded series and the padding length
#' @keywords internal
pad_series <- function(x, pad_method = "zero") {
  # Original signal length
  series_length <- length(x)

  # Calculate padding length to next power of 2 (original method)
  pot2 <- trunc(log2(series_length) + 0.5)
  pad_length <- 2^(pot2+1) - series_length

  # Apply different padding methods
  if (pad_method == "zero") {
    # Original zero-padding method
    xpad <- c(x, rep(0, pad_length))

  } else if (pad_method == "reflection") {
    # MATLAB-style reflection padding
    # Mirror half the signal at each end: [c b a | a b c d e f g | g f e]
    reflect_length <- floor(series_length/2)
    left_reflect <- x[reflect_length:1]
    right_reflect <- x[series_length:(series_length-reflect_length+1)]
    xpad <- c(left_reflect, x, right_reflect)
    # Update pad_length for return value
    pad_length <- 2 * reflect_length

  } else if (pad_method == "reflection_zero") {
    # First apply reflection padding
    reflect_length <- floor(series_length/2)
    left_reflect <- x[reflect_length:1]
    right_reflect <- x[series_length:(series_length-reflect_length+1)]
    reflected <- c(left_reflect, x, right_reflect)

    # Then calculate additional zero padding to reach next power of 2
    reflected_length <- length(reflected)
    pot2_reflected <- trunc(log2(reflected_length) + 0.5)
    zero_pad_length <- 2^(pot2_reflected+1) - reflected_length

    # Apply zero padding
    xpad <- c(reflected, rep(0, zero_pad_length))

    # Update pad_length for return value (reflection + zero padding)
    pad_length <- 2 * reflect_length + zero_pad_length

  } else {
    stop("Unknown padding method. Use 'zero', 'reflection', or 'reflection_zero'.")
  }

  return(list(xpad = xpad, pad_length = pad_length))
}
