########################## helper functionns to main.R ########################
# W_theta calculator which return Y_k, the estimation of W(theta)/2n, Sigma and the ratio W(theta)/2n/Sigma
# pos1 = m.max+2 where m.max is the order of dependence
W_theta_cal <- function(x, pos1){
  n <- length(x)
  xx  <-  c(x,x[1:pos1])
  Y  <-  numeric(pos1)
  for (k in 1:pos1){
    Y[k] <- mean((x-xx[(k+1):(n+k)])^2)/2           # calculate Y_k=T_k/2n
  }

  # use regression to find w_theta and sigma under H0
  lm_coeff <- lm(Y~c(1:pos1))$coefficients
  w_theta_hat <- lm_coeff[2]
  sigma_hat <-  lm_coeff[1]

  # w_theta_hat/sigma_hat
  ratio <- max(0, w_theta_hat/sigma_hat)

  return(list(Y = Y, ratio=ratio, sigma=sigma_hat, w_theta=w_theta_hat))
}

# autocovariance calculator
autocov_cal <- function(Y, m.max){
  # use T_k to compute autocovariance
  # compute variance first, a and b is the coefficient for computation
  T_m1 <- Y[m.max+1]
  T_m2 <- Y[m.max+2]
  var0 <- (m.max+2)*T_m1 - (m.max+1)*T_m2

  # update a and b to compute up to m.max order of autocovariance
  if (m.max >=1 ){
    acf_est <- as.numeric(NULL)
    for (i in 1:m.max) {
      a <- -(m.max+2-i)
      b <- m.max+1-i
      acf_est[i] <- -Y[i] - a*T_m1 - b*T_m2
    }
    var0 <- c(var0, acf_est)
  }
  return(var0)
}

# confidence interval calculator
h0ci_cal <- function(ratio, m.max, n, ci.level=0.95){
  # test on estimated auto correlation under H0
  coeff_1 <- rep(-1, m.max)
  coeff_2 <- rep(-1, m.max)
  for (h in 1:m.max) {
    # sigma 4th moment coeff
    coeff_1[h] <- 2*(1+(m.max+2-h)*(m.max+1-h))/n

    # w(theta) coeff
    coeff_2[h] <- 4*((m.max-h)**2+(m.max-h)*3+2)*ratio/n
  }
  h0_ci <- qnorm(0.5+ci.level/2)*sqrt(coeff_1+coeff_2)

  return(h0_ci)
}
