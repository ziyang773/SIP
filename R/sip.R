#' Shift-Immune Portmanteau (SIP) test
#'
#' SIP.test performs a Ljung–Box like autocorrelation test under frequent mean changes.
#'
#' @param seq Numeric vector. The input time series.
#' @param m Integer. Maximum lag order to test (lags 1 through \code{m}).
#' @param EVE Logical. If \code{TRUE}, uses \code{w2} from \code{W_theta_cal()} to estimate
#'   \eqn{W(\theta)/(n\,\gamma_0)}; if \code{FALSE}, uses the difference ratio
#'   \eqn{(Y_{m+2} - Y_{m+1})/\gamma_0(Y)} (\code{w1} estimator).
#' @param warnings Logical. If \code{TRUE}, warns and returns
#'   \code{p.value = 0} and \code{statistic = Inf} when the variance estimate is negative.
#'
#' @return A list with:
#' \describe{
#'   \item{statistic}{Numeric. The SIP test statistic.}
#'   \item{p.value}{Numeric. P-value from a chi-square distribution with \code{m} degrees of freedom.}
#' }
#'
#' @references Liu, Z., Hao, N., Niu, Y. S., Xiao, H., & Ding, H. (2025) Autocorrelation Test under Frequent Mean. arXiv:2510.21047
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) + c(rep(1, 50), rep(0, 50))
#' SIP.test(x, m = 4)
#'
#' @export
SIP.test <- function(seq, m=1, EVE=TRUE, warnings=TRUE){
  n <- length(seq)
  w_theta_est <- W_theta_cal(seq, pos1 = m+2)

  Y <- w_theta_est$Y
  w_theta_sigma <- w_theta_est$ratio

  # auto covariance calculation
  acf_est <- autocov_cal(Y, m)

  if (EVE == TRUE){
    if (w_theta_est$sigma<0){
      if (warnings == TRUE) warning("Variance estimation is negative, returning p.value = 0")
      return(list(p.value=0, statistic=Inf))
    } else{
      est <- acf_est[-1]/w_theta_est$sigma
    }
  } else {
    w_theta_sigma <- max(0, (Y[m+2] - Y[m+1])/acf_est[1])
    if (acf_est[1] <0) {
      if (warnings == TRUE) warning("Variance estimation is negative, returning p.value = 0")
      return(list(p.value=0, statistic=Inf))
    }
    est <- acf_est[-1]/acf_est[1]
  }

  cov_mat <- matrix(rep(0, m^2), ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      if (i==j){
        part_1 <- 2*(1+(m+2-j)*(m+1-j))
        part_2 <- 4*w_theta_sigma*(m+2-j)*(m+1-j)
        cov_mat[i, j] <- part_1 + part_2
      } else{
        k <- min(i, j)
        l <- max(i, j)
        part_1 <- 1+(m+2-j)*(m+1-i)+(m+2-i)*(m+1-j)
        part_2 <-  2*w_theta_sigma*(-4*l+2*m^2+6*m-2*l*m-2*m*k+2*k*l-2*k+4)
        cov_mat[i, j] <- part_1 + part_2
      }
    }
  }
  statistic <- n*est %*% solve(cov_mat) %*% est
  p.value <- 1-pchisq(statistic, m)
  return(list(p.value=p.value, statistic=statistic))
}


#' SIP Autocorrelation / Autocovariance (with optional plot)
#'
#' Computes (and optionally plots) sample autocorrelation or autocovariance
#' up to \code{lag.max}. Two modes are provided:
#' \itemize{
#'   \item \strong{Estimation-focused} (\code{estimation = TRUE}): uses a common
#'   \code{lag.max} for estimating all lags, typically yielding more precise
#'   estimates;
#'   \item \strong{Inference-focused} (\code{estimation = FALSE}): uses lag-specific
#'   windows for each lag, producing the same CI width across lags.
#' }
#'
#' @param x Numeric vector. Input time series (same role as \code{seq} in \code{\link{SIP.test}}).
#' @param lag.max Integer (\eqn{\ge} 0). Maximum lag to estimate/display.
#' @param type Character. Either \code{"correlation"} or \code{"covariance"}.
#'   Defaults to \code{"correlation"}.
#' @param plot Logical. If \code{TRUE}, produce the ACF-style plot with CI bands.
#'   Defaults to \code{TRUE}.
#' @param estimation Logical. If \code{TRUE}, uses the same \code{lag.max} to compute
#'   all lags (estimation-focused). If \code{FALSE}, uses lag-specific Y vector to compute each lag’s
#'   autocovariance but yields the **same** confidence interval width across lags. Defaults to \code{TRUE}.
#' @param ci.level Numeric in (0, 1). Confidence level for CI bands when \code{plot = TRUE}.
#'   Defaults to \code{0.95}.
#' @param main Optional character. Plot title.
#'
#' @return A list with components:
#' \describe{
#'   \item{acf}{Numeric vector of estimated values for lags \code{0:lag.max}.}
#'   \item{type}{The requested type (\code{"correlation"} or \code{"covariance"}).}
#'   \item{n.used}{Effective sample size used in estimation.}
#'   \item{lag}{Integer vector \code{0:lag.max}.}
#'   \item{series}{Character label for the series (if available).}
#'   \item{h0.ci}{CI half-width under the reference (white-noise) model at \code{ci.level}.}
#' }
#' If \code{plot = TRUE}, this object is returned \emph{invisibly} after plotting.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100) + c(rep(1, 50), rep(0, 50))
#' # Estimation-focused ACF (common lag.max for all lags)
#' SIP.acf(x, lag.max = 20, type = "correlation", estimation = TRUE, plot = TRUE)
#' # Inference-focused ACF (same CI width across lags)
#' SIP.acf(x, lag.max = 20, type = "correlation", estimation = FALSE, plot = TRUE)
#'
#' @seealso \code{\link{SIP.test}}, \code{\link{ACF_est}}, \code{\link{ACF_inf}}
#' @export
SIP.acf <- function(x, lag.max, type="correlation", plot=TRUE, estimation=TRUE, ci.level=0.95, main=NULL){
  if (estimation == TRUE){
    ACF_est(x, m.max=lag.max, type=type, plot=plot, ci.level=ci.level, main=main)
  } else{
    ACF_inf(x, m.max=lag.max, type=type, plot=plot, ci.level=ci.level, main=main)
  }
}


ACF_est <- function(x, m.max, type="correlation", plot=TRUE, ci.level=0.95, main=NULL){
  # m+2 <=L/2
  # get seq length and m+1 m+2 for further computation
  n=length(x)
  pos0 <- m.max+1
  pos1 <- m.max+2

  # W_theta and ratio calculation
  w_thetas <- W_theta_cal(x, pos1)
  Y <- w_thetas$Y
  W <- Y[pos1]-Y[pos0]

  # auto covariance calculation
  acf_est <- autocov_cal(Y, m.max)

  if (acf_est[1]<0){
    warning("lag.max might be too small, sample variance will be used")
    acf_est[1] <- var(x)
  }

  ratio = max(0, W/acf_est[1])

  # auto correlation
  if (type == "correlation"){
    acf_est <- acf_est/acf_est[1]
  }

  # calculate confidence interval under H0
  h0_ci <- h0ci_cal(ratio, m.max, n, ci.level=ci.level)

  # plot
  series = deparse(substitute(x))
  if (is.null(main) || is.na(main)){
    main = paste("Series", series, " ")
  }
  if (plot==TRUE) {
    if (type == "correlation") {
      plot(c(0:m.max), acf_est, type = "h", xlab = "Lag",
           ylab = "ACF", ylim = c(min(acf_est,0, -h0_ci), max(acf_est)), main = main, xaxt="n")
      axis(1, at = c(0:m.max), labels = 0:m.max)
      lines(c(1:m.max),h0_ci,col="blue", type = "o", lty=2, cex=0.1)
      lines(c(1:m.max),-h0_ci,col="blue", type = "o", lty=2, cex=0.1)
    } else if (type == "covariance") {
      plot(c(0:m.max), acf_est, type = "h", xlab = "Lag",
           ylab = "ACF (cov)", ylim = c(min(acf_est, 0), max(acf_est)), main = main, xaxt="n")
      axis(1, at = c(0:m.max), labels = 0:m.max)
    }
    abline(h=0)
    return(invisible(list(acf=acf_est, type=type, n.used=n, lag=c(0:m.max), series=series, h0.ci=h0_ci)))
  }

  return(list(acf=acf_est, type=type, n.used=n, lag=c(0:m.max), series=series, h0.ci=h0_ci))
}


ACF_inf <- function(x, m.max, type="correlation", plot=TRUE, ci.level=0.95, main=NULL){
  # m+2 <=L/2
  # get seq length and m+1 m+2 for further computation
  n=length(x)
  pos0 <- m.max+1
  pos1 <- m.max+2

  # W_theta and ratio calculation
  w_thetas <- W_theta_cal(x, pos1)
  Y <- w_thetas$Y
  W <- Y[pos1]-Y[pos0]

  # auto covariance calculation
  acf_est <- as.numeric(NULL)
  var_est <- as.numeric(NULL)
  for (i in 0:m.max){
    acf_vec <- autocov_cal(Y, i)
    if (acf_vec[1]<0){
      warning(paste("lag.max might be too small at", (i), "sample variance will be used"))
      acf_vec[1] <- var(x)
    }
    acf_est[i+1] <- acf_vec[i+1]
    var_est[i+1] <- acf_vec[1]
  }

  # auto correlation
  if (type == "correlation"){
    acf_est <- acf_est/var_est
  }
  ratio <- max(0, W/var_est[1])
  # calculate confidence interval under H0
  h0_ci <- h0ci_cal(ratio, 1, n, ci.level = ci.level)

  # plot
  series = deparse(substitute(x))
  if (is.null(main) || is.na(main)){
    main = paste("Series", series, " ")
  }
  if (plot==TRUE) {
    if (type == "correlation") {
      plot(c(0:m.max), acf_est, type = "h", xlab = "Lag",
           ylab = "ACF", ylim = c(min(acf_est,0, -h0_ci), max(acf_est)), main = main, xaxt="n")
      axis(1, at = c(0:m.max), labels = 0:m.max)
      lines(c(-1:(m.max+1)),rep(h0_ci, m.max+3),col="blue", type = "o", lty=2, cex=0.1)
      lines(c(-1:(m.max+1)),rep(-h0_ci, m.max+3),col="blue", type = "o", lty=2, cex=0.1)
    } else if (type == "covariance") {
      plot(c(0:m.max), acf_est, type = "h", xlab = "Lag",
           ylab = "ACF (cov)", ylim = c(min(acf_est, 0), max(acf_est)), main = main, xaxt="n")
      axis(1, at = c(0:m.max), labels = 0:m.max)
    }
    abline(h=0)
    return(invisible(list(acf=acf_est, type=type, n.used=n, lag=c(0:m.max), series=series, h0.ci=h0_ci)))
  }

  return(list(acf=acf_est, type=type, n.used=n, lag=c(0:m.max), series=series, h0.ci=h0_ci))
}
