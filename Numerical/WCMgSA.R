wcm.gsa <- function(x, R = 100, min.len = NULL,
                    gappy = TRUE, Q = floor(log(length(x))^1.9), max.iter = 5,
                    p.max = 10, pen = log(length(x))^1.01) {

  n <- length(x)
  if(is.null(min.len)) min.len <- max(20, p.max + ceiling(log(n)^1.1))

  res <- t(wbs2(x, R, min.len, 1))
  res <- res[order(- res[, 4]), ]
  colnames(res) <- c('s', 'b', 'e', 'cusum')

  Q <- min(Q, dim(res)[1])
  z <- log(res[1:Q, 4])
  cs <- cumsum(z[1:Q])
  ecp.seq <- res[1:Q, 2]

  ecp.list <- list()
  current <- c()
  if(gappy){
    niter <- ii <- 0
    for(k in sort(sort(abs(diff(z)), decreasing = TRUE, index.return = TRUE)$ix[1:max.iter])){
      new.cp <- sort(ecp.seq[(ii + 1):k])
      current <- sort(c(current, new.cp))
      ecp.list <- c(ecp.list, list(current))
      ii <- k
    }
  } else{
    for(ii in 1:(Q - 1)){
      current <- sort(c(current, ecp.seq[ii]))
      ecp.list <- c(ecp.list, list(current))
    }
  }

  ll <- length(ecp.list)
  while(ll >= 1){
    if(ll == 1) current <- numeric(0) else current <- ecp.list[[ll - 1]]
    new.cp <- setdiff(ecp.list[[ll]], current)

    out <- rep(0, length(current) + 1)
    brks <- c(0, current, n)
    for(jj in 1:(length(current) + 1)){
      int <- (brks[jj] + 1):brks[jj + 1]
      ecp <- new.cp[new.cp %in% int]
      if(length(ecp) > 0){
        loc.brks <- c(brks[jj], ecp, brks[jj + 1])
        if(length(int) > p.max + log(n) + length(ecp) + 1 & max(diff(loc.brks)) > p.max + 1 + log(n)){
          ms <- ms.sc(x, s = brks[jj], e = brks[jj + 1], ecp = ecp, p.max = p.max, pen = pen)
          if(ms$cp.sc >= ms$ncp.sc) out[jj] <- -1 else out[jj] <- 1
        } else out[jj] <- -1
      }
    }
    if(any(out < 0)) flag <- FALSE else{
      rcm <- refine.chp.model(ecp.list[[ll]], x, c = .5)
      ms <- ms.sc(x, s = 0, e = n, ecp = rcm$rcp, p.max = p.max, pen = pen)
      if(ms$cp.sc >= ms$ncp.sc) flag <- FALSE else flag <- TRUE
    }
    if(flag) break
    ll <- ll - 1
  }

  if(ll == 0){
    rcp <- cp <- integer(0)
    cp.info <- matrix(0, nrow = 0, ncol = 4)
    colnames(cp.info) <- colnames(res)
  } else{
    cp <- ecp.list[[ll]]
    rcp <- refine.chp.model(cp, x, c = .5)$rcp
    cp.info <- res[1:length(cp),, drop = FALSE]
  }

  return(list(cp = cp, rcp = rcp, cp.info = cp.info))
}

## internal functions -- not to be called directly by the user

#' WBS2
#' @keywords internal
wbs2 <- function(x, R = 100, min.len = 0, top.scale = 0){

  n <- length(x)
  if (n <= 2 * min.len + 1) return(matrix(NA, 4, 0)) else {
    if(top.scale == 1){
      tmp <- rbind(systematic.cusums(x[1:round(n/2)], R, min.len)$max.val,
                   systematic.cusums(x[round(n/4) + 1:round(n/2)], R, min.len)$max.val + c(rep(round(n/4), 3), 0),
                   systematic.cusums(x[(round(n/2) + 1):n], R, min.len)$max.val + c(rep(round(n/2), 3), 0))
      cpt <- t(tmp[which.max(tmp[, 4]), , drop = FALSE])
    } else cpt <- t(systematic.cusums(x, R, min.len)$max.val)
    return(cbind(cpt, wbs2(x = x[1:cpt[2]], R, min.len, 0),
                 wbs2(x = x[(cpt[2] + 1):n], R, min.len, 0) + c(rep(cpt[2], 3), 0)))
  }

}

systematic.cusums <- function(x, R, min.len) {

  y <- c(0, cumsum(x))
  n <- length(x)
  R <- min(R, (n - 1)*n/2)

  ind <- grid.intervals(n, R)
  R <- dim(ind)[2]
  res <- matrix(0, R, 4)

  res[,1:2] <- t(ind)
  res[,3:4] <- t(apply(ind, 2, max.cusum, y, min.len))

  max.ind <- which.max(abs(res[,4]))
  max.val <- res[max.ind, c(1, 3, 2, 4), drop = FALSE]

  list(res = res, max.val = max.val, R.eff = R)

}

grid.intervals <- function(n, R){

  if(n == 2){
    ind <- matrix(c(1, 2), 2, 1)
  } else if(R >= (n-1)*n/2){
    ind <- all.intervals.flat(n)
  } else{
    k <- 1
    while (k*(k-1)/2 < R) k <- k + 1
    ind2 <- all.intervals.flat(k)
    ind2.mx <- max(ind2)
    ind <- round((ind2 - 1) * ((n - 1) / (ind2.mx - 1)) + 1)
  }

  ind
}

all.intervals.flat <- function(n){

  if (n == 2) ind <- matrix(1:2, 2, 1) else {
    R <- (n - 1)*n/2
    ind <- matrix(0, 2, R)
    ind[1,] <- rep(1:(n - 1), (n - 1):1)
    ind[2,] <- 2:(R + 1) - rep(cumsum(c(0, (n - 2):1)), (n - 1):1)
  }
  ind

}

max.cusum <- function(ind, y, min.len) {

  m <- ind[2] - ind[1] + 1
  if(m > 2 * min.len + 1){
    z <- y[(ind[1] + 1):(ind[2] + 1)] - y[ind[1]]
    ip <- sqrt(((m - 1):1) / m / (1:(m - 1))) * z[1:(m - 1)] - sqrt((1:(m - 1))/m/((m - 1):1)) * (z[m] - z[1:(m - 1)])
    ip.max <- which.max(abs(ip[(min.len + 1):(m - min.len - 1)])) + min.len
    return(c(ip.max + ind[1] - 1, abs(ip[ip.max])))
  } else{
    return(c(ind[1], 0))
  }

}

#' gSA: model selection using Schwarz criterion
#' @keywords internal
ms.sc <- function(x, s, e, ecp, p.max = 10, pen = NULL){

  n <- length(x)

  brks <- c(s, ecp, e)
  cp.lik.seq <- ncp.lik.seq <- cp.sc.seq <- ncp.sc.seq <- rep(0, p.max + 1)

  for(ii in 1:(p.max + 1)){
    pp <- ii - 1

    len <- e - max(pp, s)
    y <- x[(max(pp, s) + 1):e]
    R <- L <- c()
    if(pp > 0) for(kk in 1:pp) L <- cbind(L, x[(max(pp, s) + 1):e - kk])
    for(kk in 1:(length(brks) - 1)){
      if(brks[kk + 1] <= pp) next
      tmp <- rep(0, len)
      tmp[max(1, brks[kk] - max(pp, s) + 1):(brks[kk + 1] - max(pp, s))] <- 1
      R <- cbind(R, tmp)
    }
    X <- cbind(L, R)
    lm.fit1 <- lm(y ~ 0 + X)
    cp.lik.seq[ii] <- loglik1 <- -len/2*log(2*pi*sum(resid(lm.fit1)^2)/len) - len/2
    cp.sc.seq[ii] <- (pp + length(ecp) + 1) * pen - loglik1

    if(pp > 0) z <- y - L %*% coef(lm.fit1)[1:pp] else z <- y

    ncp.lik.seq[ii] <- loglik2 <- -len/2*log(2*pi*sum((z - mean(z))^2)/len) - len/2
    ncp.sc.seq[ii] <- (pp + 1) * pen - loglik2
  }

  p <- which.min(cp.sc.seq[-1])

  out <- list(cp.sc = cp.sc.seq[p + 1], ncp.sc = ncp.sc.seq[p + 1], p = p,
              cp.lik = cp.lik.seq, ncp.lik = ncp.lik.seq,
              cp.sc.seq = cp.sc.seq, ncp.sc.seq = ncp.sc.seq)
  return(out)
}

#' Refinement of change point location estimators
#' @keywords internal
refine.chp.model <- function(ecp, x, c = .5){

  n <- length(x)
  rcp <- ecp * 0
  nres <- c()
  for(i in 1:length(ecp)){
    if(i == 1){
      if(i < length(ecp)) int <- 1:floor(ecp[i]*(1 - c) + ecp[i + 1]*c) else int <- 1:n
    } else{
      if(i < length(ecp)) int <- floor(ecp[i]*(1 - c) + ecp[i - 1]*c + 1):floor(ecp[i]*(1 - c) + ecp[i + 1]*c) else int <- floor(ecp[i]*(1 - c) + ecp[i - 1]*c + 1):n
    }
    if(length(int) < 2){
      if(i == 1) int <- 1:ecp[i + 1] else if(i == length(ecp)) int <- (ecp[i - 1] + 1):n else int <- (ecp[i - 1] + 1):ecp[i + 1]
    }
    csx <- cumsum(x[int])
    len <- length(int)
    evalx <- (len/((1:(len - 1)) * ((len - 1):1)))^.5 * (csx[1:(len - 1)] - csx[len]/len * (1:(len - 1)))
    rcp[i] <- int[1] + which.max(abs(evalx)) - 1
    nres <- rbind(nres, c(min(int), rcp[i], max(int), max(abs(evalx))))
  }

  return(list(rcp = rcp, nres = nres))
}
