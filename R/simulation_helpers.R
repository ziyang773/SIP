#' Simulate Piecewise-Constant Mean Structure with Random Shifts
#'
#' Generates a length-\code{n} vector of segment means for a time series with
#' \code{n.seg} contiguous segments, each of length at least \code{min.seg.len}.
#' Segment means are drawn i.i.d. from a continuous uniform distribution over
#' \code{mean.range}. The result can be added to noise to create data with
#' frequent mean shifts.
#'
#' @param n Integer. Total length of the simulated sequence.
#' @param n.seg Integer. Number of contiguous segments.
#' @param min.seg.len Integer. Minimum length of each segment.
#' @param mean.range Numeric length-2 vector \code{(low, high)} giving the range
#'   from which segment means are drawn (uniformly).
#'
#' @return A numeric vector of length \code{n} containing the piecewise-constant
#'   mean structure; \code{NULL} if inputs are infeasible.
#'
#' @export
mean_sim <- function(n=1000, n.seg=10, min.seg.len=20, mean.range=c(-10, 10)){
  # seg.len = n/n.seg >=30
  if (n < n.seg*min.seg.len){
    warning("The length of the seq is not long enough")
    return(NULL)
  }

  # random selection and add some shift
  # n.seg should be greater to 1
  if (n.seg < 1){
    warning("n.seg should be greater than 1")
    return(NULL)
  }
  # select potential start point for each segment (exclude 1)
  seg_points <- sort(sample(1:(n-n.seg*min.seg.len), n.seg-1, replace = TRUE))
  seg_points <- seg_points + min.seg.len*(1:(n.seg-1))
  seg_points <- c(1, seg_points, n+1)

  # add mean structure
  seg_means <- sample(runif(n.seg, mean.range[1], mean.range[2])) # modified to continuous 10/04
  mean_structure <- rep(seg_means, diff(seg_points))

  return(mean_structure)
}

############################################## noise type ##################################################
noise_sim <- function(n, type="normal"){
  if (type == "normal"){
    noise <- rnorm(n, sd=1)
  }
  else if (type == "t"){
    noise <- rt(n, df=3)/sqrt(3)
  }
  else if (type == "exponential"){
    noise <- rexp(n, rate = 1) - 1
  }
  return(noise)
}

############################################## autocorrelation ##################################################
#' Simulate Noise with Target ACF Structure
#'
#' Generates a length-\code{n} series under one of \code{"iid"}, \code{"AR"},
#' \code{"MA"}, or \code{"ARMA"} (ARMA(1,1)) models with selectable noise type.
#'
#' @param n Integer. Length of the output series.
#' @param model Character. One of \code{"iid"}, \code{"AR"}, \code{"MA"}, \code{"ARMA"}.
#' @param noise.type Character. Noise distribution passed to \code{noise_sim()} (e.g., \code{"normal"}).
#' @param ar.coef Numeric. AR coefficients (vector for AR; scalar for ARMA(1,1)).
#' @param ma.coef Numeric. MA coefficients (vector for MA; scalar for ARMA(1,1)).
#'
#' @return A numeric vector of length \code{n}.
#'
#' @seealso \code{\link{noise_sim}}, \code{\link{SIP.acf}}, \code{\link{SIP.test}}
#' @importFrom stats filter na.omit
#' @export
acf_sim <- function(n=1000, model="iid", noise.type="normal", ar.coef=NULL, ma.coef=NULL){
  # noise
  noise <- noise_sim(n, type=noise.type)

  # iid
  if (model == "iid"){
    return(noise)
  }

  # AR noise
  if (model == "AR"){
    if (is.null(ar.coef)){
      print("You have to specify the coefficients for AR process")
      return(NULL)
    }
    ar_noise <- filter(noise, filter = ar.coef, method = "recursive")
    return(ar_noise)
  }

  # ARMA process, only support ARMA(1,1)
  else if (model == "ARMA"){
    if (is.null(ar.coef) || is.null(ma.coef)){
      print("You have to specify the coefficients for ARMA process")
      return(NULL)
    }
    arma_noise <- rep(0, n)
    arma_noise[1] <- noise[1]
    for(t in 2:n){
      arma_noise[t] <- ar.coef*arma_noise[t-1] + noise[t] + ma.coef*noise[t-1]
    }
    return(arma_noise)
  }

  # MA process
  else if (model == "MA"){
    if (is.null(ma.coef)){
      print("You have to specify the coefficients for MA process")
      return(NULL)
    }
    noise <- noise_sim(n+length(ma.coef), type=noise.type)
    ma_noise <- filter(noise, filter = c(1, ma.coef), sides = 1)
    ma_noise <- na.omit(ma_noise)
    return(ma_noise)
  }
}

############################################## whole sequence ##################################################
seq_simulator <- function(n=1000, model, ar.coef=NULL, ma.coef=NULL, n.seg, min.seg.len=20, mean.range=c(-10, 10), noise.type="normal"){
  mean_structure <- mean_sim(n=n, n.seg=n.seg, min.seg.len=min.seg.len, mean.range = mean.range)
  acf_structure <- acf_sim(n=n, model=model, ar.coef=ar.coef, ma.coef=ma.coef, noise.type=noise.type)
  seq <- mean_structure+acf_structure

  # estimate acf with R default acf function
  acf_true <- acf(acf_structure, lag.max = 10, plot = F)$acf

  # return a list
  return(list(sequence = seq, mean_structure = mean_structure, acf_structure = acf_structure,
              n = n, n.seg=n.seg, model=model, acf_est=acf_true))
}

# example
seq1 <- seq_simulator(10000, model = 'MA', n.seg=10, ma.coef = c(0.4))
