devtools::install_github("ziyang773/SIP")
library(SIP)
#############################################  Simulation  ###############################################
##################################### Table 1. type one error ##################################################
set.seed(2024)
mean_structure <- mean_sim(10000, 100, 20, c(-5, 5))
cpts = c(0, which(diff(mean_structure)!=0), 10000)
type_one_error <- as.data.frame(matrix(rep(0, 60), ncol=4))
colnames(type_one_error) <- paste("lag",c(1,2,4,8))
row.names(type_one_error) <- paste(rep(c("SIP1", "SIP2", "Box",
                                         "oracle", 'pseudo oracle'), 3), rep(c("G", "E", "T"), each = 5), sep = "-")
noise_type <- c("normal", "exponential", "t")
lags <- c(1,2,4,8)

for (k in 1:3) {
  type <- noise_type[k]
  set.seed(2024)
  for (i in 1:1000) {
    noise_structure <- acf_sim(10000, model = "iid", noise.type = type)
    seq <- mean_structure+noise_structure
    for (j in 1:4) {
      lag <- lags[j]
      type_one_error[5*k-4, j] = type_one_error[5*k-4, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE)$p.value<=0.05)

      type_one_error[5*k-3, j] = type_one_error[5*k-3, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE)$p.value<=0.05)

      type_one_error[5*k-2, j] = type_one_error[5*k-2, j] + as.numeric(Box.test(seq, lag, type = "Box-Pierce")$p.value<=0.05)

      type_one_error[5*k-1, j] = type_one_error[5*k-1, j] + as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value<=0.05)

      demeaned_seq <- as.numeric(10000)
      for (l in 2:length(cpts)){
        start <- cpts[l-1]
        end <- cpts[l]
        demeaned_seq[(start+1):(end)] <- seq[(start+1):(end)] - mean(seq[(start+1):end])
      }
      type_one_error[5*k, j] = type_one_error[5*k, j] + as.numeric(Box.test(demeaned_seq, lag, type = "Box-Pierce")$p.value<=0.05)
    }
  }
}
type_one_error <- type_one_error/1000
type_one_error


##################################### Table 2-4. power comparison ##################################################
set.seed(2024)
mean_structure <- mean_sim(10000, 100, 20, c(-5, 5))
# # MA(1)
power_MA1 <- as.data.frame(matrix(rep(0, 3*6*4), ncol=4))
colnames(power_MA1) <- c(paste("lag", c(1,2,4,8)))
coeffs <- c(-0.1, -0.05, -0.025, 0.025, 0.05, 0.1)
row.names(power_MA1) <-  paste(paste("\u03C9=", rep(coeffs, each = 3), sep = ""), rep(c("SIP1", "SIP2", "Oracle"), 6), sep = ",")
lags <- c(1,2,4,8)
set.seed(2024)
for (j in 1:6) {
  coeff <- coeffs[j]
  for (i in 1:1000) {
    noise_structure <- acf_sim(10000, model="MA", ma.coef = c(coeff))
    seq <- mean_structure+noise_structure
    for (k in 1:4) {
      lag <- lags[k]
      power_MA1[3*j-2, k] = power_MA1[3*j-2, k] + as.numeric(SIP.test(seq, lag, EVE=FALSE)$p.value<=0.05)
      power_MA1[3*j-1, k] = power_MA1[3*j-1, k] + as.numeric(SIP.test(seq, lag, EVE=TRUE)$p.value<=0.05)
      power_MA1[3*j, k] = power_MA1[3*j, k] + as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value<=0.05)
    }
  }
}
power_MA1 <- power_MA1/1000
power_MA1

# MA(4)
power_MA4 <- as.data.frame(matrix(rep(0,9*4), ncol=4))
colnames(power_MA4) <- c(paste("lag", c(1,2,4,8)))
row.names(power_MA4) <- paste(rep(c("scenario 1", "scenario 2", "scenario 3"), each=3),  rep(c("SIP1", "SIP2", "oracle"), 3), sep = ",")
lags <- c(1,2,4,8)
set.seed(2024)
for (i in 1:1000) {
  # scenario 1
  noise_structure <- acf_sim(10000, model="MA", ma.coef=c(0.5,0.4,0.3,0.2))
  seq <- mean_structure + noise_structure
  for (j in 1:4) {
    lag <- lags[j]
    power_MA4[1, j] = power_MA4[1, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE, warnings = FALSE)$p.value<=0.05)
    power_MA4[2, j] = power_MA4[2, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE, warnings = FALSE)$p.value<=0.05)
    power_MA4[3, j] = power_MA4[3, j]+as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value<=0.05)
  }

  # scenario 2
  noise_structure <- acf_sim(10000, model = "MA", ma.coef = c(0.1,0.1,0.5,-0.4))
  seq <- mean_structure + noise_structure
  for (j in 1:4) {
    lag <- lags[j]
    power_MA4[4, j] = power_MA4[4, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE, warnings = FALSE)$p.value<=0.05)
    power_MA4[5, j] = power_MA4[5, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE, warnings = FALSE)$p.value<=0.05)
    power_MA4[6, j] = power_MA4[6, j]+as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value<=0.05)
  }

  # scenario 3
  noise_structure <- acf_sim(10000, model = "MA", ma.coef = c(0, 0.1, 0, -0.8))
  seq <- mean_structure + noise_structure
  for (j in 1:4) {
    lag <- lags[j]
    power_MA4[7, j] = power_MA4[7, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE, warnings = FALSE)$p.value<=0.05)
    power_MA4[8, j] = power_MA4[8, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE, warnings = FALSE)$p.value<=0.05)
    power_MA4[9, j] = power_MA4[9, j]+as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value<=0.05)
  }
}
power_MA4 <- power_MA4/1000
power_MA4$`lag 4`<- formatC(power_MA4$`lag 4`, format = "f", digits = 3)
power_MA4$`lag 8`<- formatC(power_MA4$`lag 8`, format = "f", digits = 3)
power_MA4


# AR(1)
power_AR1 <- as.data.frame(matrix(rep(0,18*4), ncol=4))
lags <- c(1,2,4,8)
ar_coeff <- c(-0.1, -0.05,-0.025, 0.025, 0.05, 0.1)
colnames(power_AR1) <- c(paste("lag", c(1,2,4,8)))
row.names(power_AR1) <- paste(paste("\u03C6=", rep(ar_coeff, each = 3), sep = ""), rep(c("SIP1", "SIP2", "oracle"), 6), sep = ",")
set.seed(2024)
for (k in 1:6){
  ar = ar_coeff[k]
  for (i in 1:1000) {
    for (j in 1:4) {
      lag <- lags[j]
      noise_structure <- acf_sim(10000, model = "AR", ar.coef = c(ar))
      seq <- mean_structure + noise_structure
      power_AR1[3*k-2, j] = power_AR1[3*k-2, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE, warnings = FALSE)$p.value<=0.05)
      power_AR1[3*k-1, j] = power_AR1[3*k-1, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE, warnings = FALSE)$p.value<=0.05)
      power_AR1[3*k, j] = power_AR1[3*k, j] + as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value<=0.05)
    }
  }
}
power_AR1 <- power_AR1/1000
power_AR1
# save(type_one_error, power_MA1, power_MA4, power_AR1, file = "results.Rdata")

##################################### Table 6. Type one error comparison with other methods ##################################################
# install.packages('dCovTS')
source('testMethods.R')
library(dCovTS)

set.seed(2024)
mean_structure <- mean_sim(10000, 100, 20, c(-5, 5))
cpts = c(0, which(diff(mean_structure)!=0), 10000)
type_one_error <- as.data.frame(matrix(rep(0, 96), ncol=4))
colnames(type_one_error) <- paste("lag",c(1,2,4,8))
row.names(type_one_error) <- paste(rep(c("SIP1", "SIP2", "Box",
                                         "oracle", 'pseudo oracle', 'MMB', 'DWB', 'UnivTest'), 3),
                                   rep(c("G", "E", "T"), each = 8), sep = "-")
noise_type <- c("normal", "exponential", "t")
lags <- c(1,2,4,8)

for (k in 1:3) {
  type <- noise_type[k]
  for (i in 1:100) {
    noise_structure <- acf_sim(10000, model = "iid", noise.type = type)
    seq <- mean_structure+noise_structure
    for (j in 1:4) {
      lag <- lags[j]
      type_one_error[8*k-7, j] = type_one_error[8*k-7, j] + as.numeric(SIP.test(seq, lag, EVE=FALSE)$p.value<=0.05)

      type_one_error[8*k-6, j] = type_one_error[8*k-6, j] + as.numeric(SIP.test(seq, lag, EVE=TRUE)$p.value<=0.05)

      type_one_error[8*k-5, j] = type_one_error[8*k-5, j] + as.numeric(Box.test(seq, lag, type = "Box-Pierce")$p.value<=0.05)

      type_one_error[8*k-4, j] = type_one_error[8*k-4, j] + as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value<=0.05)

      demeaned_seq <- numeric(10000)
      for (l in 2:length(cpts)){
        start <- cpts[l-1]
        end <- cpts[l]
        demeaned_seq[(start+1):(end)] <- seq[(start+1):(end)] - mean(seq[(start+1):end])
      }
      type_one_error[8*k-3, j] = type_one_error[8*k-3, j] + as.numeric(Box.test(demeaned_seq, lag, type = "Box-Pierce")$p.value<=0.05)

      if (lag==1){
        mbb_res <- mbb_rt96(seq, lag = 1, b = 25, R = 500)
        type_one_error[8*k-2, j] = type_one_error[8*k-2, j] + as.numeric(mbb_res$p_value<=0.05)
      } else{
        mbb_res <- mbb_multi_lag(seq, lag = 1:lag, b = 25, R = 500)
        type_one_error[8*k-2, j] = type_one_error[8*k-2, j] + as.numeric(mbb_res$p_value<=0.05)
      }

      dwb_res <- dwb_dependence_test(seq, L = lag, B = 500, bandwidth = 25)
      type_one_error[8*k-1, j] = type_one_error[8*k-1, j] + as.numeric(dwb_res$p_value<=0.05)


      univ_res <- UnivTest(
        x = seq,
        type = "truncated",
        testType = "covariance",
        p = lag,        # bandwidth / max lag
        b = 500,       # number of bootstrap replicates
        parallel = FALSE,
        bootMethod = "Wild Bootstrap"
      )
      type_one_error[8*k, j] = type_one_error[8*k, j] + as.numeric(univ_res$p.value<=0.05)
    }
    cat("Finished Noise Type:", type, "- Iteration:", i, "\n")
  }
}
type_one_error <- type_one_error/100
type_one_error

##################################### Figure 3-4. Power curves under MA(1) and AR(1) noise ##################################################
set.seed(2024)

# # MA(1)
# Parameters
lengths <- c(2000, 5000, 10000, 20000)
lag <- 4
coeffs <- seq(-0.2, 0.2, by = 0.01)
n_sims <- 1000 # Number of simulations per coefficient

# Initialize a list to store results for each length
power_list <- list()

for (n in lengths) {

  # Initialize the matrix for the current sample size
  MA_1_power <- matrix(0, nrow = length(coeffs), ncol = 3)
  colnames(MA_1_power) <- c("SIP1", "SIP2", "Oracle")
  rownames(MA_1_power) <- paste0("ω=", coeffs)

  # Generate mean structure once for the current length
  # Assuming mean_sim and acf_sim scale with n
  mean_structure <- mean_sim(n, n/100, 20, c(-5, 5))

  for (j in 1:length(coeffs)) {
    coeff <- coeffs[j]

    for (i in 1:n_sims) {
      noise_structure <- acf_sim(n, model="MA", ma.coef = c(coeff))
      seq_data <- mean_structure + noise_structure

      # SIP1
      MA_1_power[j, 1] <- MA_1_power[j, 1] +
        as.numeric(SIP.test(seq_data, lag, EVE=FALSE)$p.value <= 0.05)

      # SIP2
      MA_1_power[j, 2] <- MA_1_power[j, 2] +
        as.numeric(SIP.test(seq_data, lag, EVE=TRUE)$p.value <= 0.05)

      # Oracle (Box-Pierce)
      MA_1_power[j, 3] <- MA_1_power[j, 3] +
        as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value <= 0.05)
    }
  }

  # Average the results to get power and store in the list
  power_list[[paste0("MA_1_power_", n)]] <- MA_1_power / n_sims
}
cols <- c("red", "green", "black")
par(oma = c(6, 2, 2, 2))
layout(matrix(1:4, nrow = 2, byrow = TRUE))

sample_sizes <- c(2000, 5000, 10000, 20000)
ylim_all <- range(unlist(power_list))

for (i in 1:4) {
  par(mar = c(5, 6, 3, 2))

  matplot(coeffs, power_list[[i]], type = "l", lty = 1, lwd = 2,
          col = cols,
          xlab = "",
          ylab = "Empirical Power",
          ylim = ylim_all,
          main = "",
          cex.axis = 1.6,
          cex.lab = 1.6)

  abline(h = 0.05, lty = 2, col = 'blue')
  text(x = min(coeffs), y = 0.05, labels = "0.05",
       col = "blue", pos = 3, cex = 1.4)

  mtext(paste("n =", sample_sizes[i]),
        side = 1, line = 3.8, cex = 1.5)

  if (i == 4) {
    legend("bottomright",
           legend = colnames(power_list[[1]]),
           col = cols,
           lty = 1,
           lwd = 2,
           bty = "n",
           cex = 1.3,
           y.intersp = 0.2,
           x.intersp = 0.2,
           seg.len = 0.6,
           inset = c(-0.35, -0.40),
           xpd = TRUE,
           adj = c(0, 0.5))
  }
}

mtext(expression("MA(1) coefficient (" * omega * ")"),
      side = 1, line = 3, outer = TRUE, cex = 1.8)
dev.off()
#####################
# # AR(1)
set.seed(2024)

# Parameters
lengths <- c(2000, 5000, 10000, 20000)
lag <- 4
coeffs <- seq(-0.2, 0.2, by = 0.01)
n_sims <- 1000

# Initialize list for results
AR_1_power_results <- list()

for (n in lengths) {

  # Initialize the matrix for the current sample size n
  AR_1_power <- matrix(0, nrow = length(coeffs), ncol = 3)
  colnames(AR_1_power) <- c("SIP1", "SIP2", "Oracle")
  rownames(AR_1_power) <- paste0("\u03C6=", coeffs)

  # Re-generate mean structure for current length
  # Note: Adjusting the 2nd param (100) if it needs to scale with n
  mean_structure <- mean_sim(n, n/100, 20, c(-5, 5))

  for (j in 1:length(coeffs)) {
    coeff <- coeffs[j]

    for (i in 1:n_sims) {
      # Simulate AR noise
      noise_structure <- acf_sim(n, model="AR", ar.coef = c(coeff))
      seq_data <- mean_structure + noise_structure

      # Test calculations
      AR_1_power[j, 1] <- AR_1_power[j, 1] +
        as.numeric(SIP.test(seq_data, lag, EVE=FALSE)$p.value <= 0.05)

      AR_1_power[j, 2] <- AR_1_power[j, 2] +
        as.numeric(SIP.test(seq_data, lag, EVE=TRUE)$p.value <= 0.05)

      AR_1_power[j, 3] <- AR_1_power[j, 3] +
        as.numeric(Box.test(noise_structure, lag, type="Box-Pierce")$p.value <= 0.05)
    }
  }

  # Scale by number of simulations and store
  AR_1_power_results[[paste0("AR_1_power_", n)]] <- AR_1_power / n_sims
}
cols <- c("red", "green", "black")
par(oma = c(6, 2, 2, 2))
layout(matrix(1:4, nrow = 2, byrow = TRUE))

sample_sizes <- c(2000, 5000, 10000, 20000)
ylim_all <- range(unlist(AR_1_power_results))

for (i in 1:4) {
  par(mar = c(5, 6, 3, 2))

  matplot(coeffs, AR_1_power_results[[i]], type = "l", lty = 1, lwd = 2,
          col = cols,
          xlab = "",
          ylab = "Empirical Power",
          ylim = ylim_all,
          main = "",
          cex.axis = 1.6,
          cex.lab = 1.6)

  abline(h = 0.05, lty = 2, col = 'blue')
  text(x = min(coeffs), y = 0.05, labels = "0.05",
       col = "blue", pos = 3, cex = 1.4)

  mtext(paste("n =", sample_sizes[i]),
        side = 1, line = 3.8, cex = 1.5)

  if (i == 4) {
    legend("bottomright",
           legend = colnames(AR_1_power_results[[1]]),
           col = cols,
           lty = 1,
           lwd = 2,
           bty = "n",
           cex = 1.3,
           y.intersp = 0.2,
           x.intersp = 0.2,
           seg.len = 0.6,
           inset = c(-0.35, -0.40),
           xpd = TRUE,
           adj = c(0, 0.5))
  }
}

mtext("AR(1) coefficient (\u03D5)",
      side = 1, line = 2, outer = TRUE, cex = 1.8)

dev.off()
################################# Table 7 ARMA(1, 1) Noise ########################################
set.seed(2024)
# Parameters
sample_sizes <- c(2000, 5000, 10000, 20000)
lags <- c(1, 2, 4, 8)
methods <- c("SIP1", "SIP2", "oracle")
n_sims <- 1000

row_names <- as.vector(t(outer(paste0("n=", sample_sizes), methods, paste, sep="_")))
ARMA_table <- matrix(0, nrow = length(row_names), ncol = length(lags))
rownames(ARMA_table) <- row_names
colnames(ARMA_table) <- paste0("lag", lags)

for (s in seq_along(sample_sizes)) {
  n <- sample_sizes[s]
  cat("Processing n =", n, "...\n")

  mean_structure <- mean_sim(n, n/100, 20, c(-5, 5))

  row_offset <- (s - 1) * 3

  for (i in 1:n_sims) {
    for (j in seq_along(lags)) {
      lag <- lags[j]

      # Simulate noise
      noise_structure <- acf_sim(n, model = "ARMA", ar.coef = 0.05, ma.coef = 0.05)

      # Combined sequence
      seq_data <- mean_structure + noise_structure

      # 1. SIP1 (EVE = FALSE) - Row offset + 1
      ARMA_table[row_offset + 1, j] <- ARMA_table[row_offset + 1, j] +
        as.numeric(SIP.test(seq_data, lag, EVE = FALSE, warnings = FALSE)$p.value <= 0.05)

      # 2. SIP2 (EVE = TRUE) - Row offset + 2
      ARMA_table[row_offset + 2, j] <- ARMA_table[row_offset + 2, j] +
        as.numeric(SIP.test(seq_data, lag, EVE = TRUE, warnings = FALSE)$p.value <= 0.05)

      # 3. Box-Pierce (Oracle) - Row offset + 3
      ARMA_table[row_offset + 3, j] <- ARMA_table[row_offset + 3, j] +
        as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value <= 0.05)
    }
  }
}
ARMA_table <- ARMA_table / 1000
ARMA_table
################################ bilinear ######################################
bilinear_sim <- function(n, a, b, burnin = 500, sd = 1) {
  # total length including burn-in
  N <- n + burnin

  # innovations
  eps <- rnorm(N, mean = 0, sd = sd)

  # initialize series
  X <- numeric(N)
  X[1] <- eps[1]  # or 0, both are fine

  # recursion
  for (t in 2:N) {
    X[t] <- (a + b * eps[t]) * X[t - 1] + eps[t]
  }

  # remove burn-in
  return(X[(burnin + 1):N])
}

set.seed(2024)
sample_sizes <- c(2000, 5000, 10000, 20000)
lags <- c(1, 2, 4, 8)
methods <- c("SIP1", "SIP2", "oracle")
n_sims <- 1000
# Initialize the transposed table: Rows = (n * Method), Cols = Lags
row_names <- as.vector(t(outer(paste0("n=", sample_sizes), methods, paste, sep="_")))
Bilinear <- matrix(0, nrow = length(row_names), ncol = length(lags))
rownames(Bilinear) <- row_names
colnames(Bilinear) <- paste0("lag", lags)

for (s in seq_along(sample_sizes)) {
  n <- sample_sizes[s]
  cat("Processing bilinear simulations for n =", n, "...\n")

  mean_structure <- mean_sim(n, n/100, 20, c(-5, 5))

  row_offset <- (s - 1) * 3

  for (i in 1:n_sims) {
    noise_structure <- bilinear_sim(n=n, a=0.1, b=0.5)
    seq_data <- mean_structure + noise_structure

    for (j in seq_along(lags)) {
      lag <- lags[j]

      # 1. SIP1 (EVE = FALSE)
      Bilinear[row_offset + 1, j] <- Bilinear[row_offset + 1, j] +
        as.numeric(SIP.test(seq_data, lag, EVE = FALSE, warnings = FALSE)$p.value <= 0.05)

      # 2. SIP2 (EVE = TRUE)
      Bilinear[row_offset + 2, j] <- Bilinear[row_offset + 2, j] +
        as.numeric(SIP.test(seq_data, lag, EVE = TRUE, warnings = FALSE)$p.value <= 0.05)

      # 3. Box-Pierce (Oracle)
      Bilinear[row_offset + 3, j] <- Bilinear[row_offset + 3, j] +
        as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value <= 0.05)
    }
  }
}
Bilinear <- Bilinear / 1000
Bilinear
################################ Table 9-10 Power Structual Change   ###############################################
# # MA(1)
set.seed(2024)

mean_values <- c(5, 10, 20)
n_segs <- c(50, 100, 200)
lags <- c(1, 2, 4, 8)
methods <- c("SIP1", "SIP2", "oracle")

grid <- expand.grid(Method = methods, Seg = n_segs, Mean = mean_values)
grid <- grid[order(grid$Mean, grid$Seg), ]

row_labels <- paste0("Mean=", grid$Mean, "_Seg=", grid$Seg, "_", grid$Method)
MA1_sc_power <- matrix(0, nrow = length(row_labels), ncol = length(lags))
rownames(MA1_sc_power) <- row_labels
colnames(MA1_sc_power) <- paste0("lag", lags)

# Counter for the current row group
row_idx <- 1

for (m_val in mean_values) {
  for (seg in n_segs) {
    cat("Processing: Mean =", m_val, "| Segments =", seg, "...\n")

    mean_structure <- mean_sim(10000, seg, 20, c(-m_val, m_val))
    for (i in 1:1000) {
      noise_structure <- acf_sim(10000, model = "MA", ma.coef = 0.1)
      seq_data <- mean_structure + noise_structure

      for (k in seq_along(lags)) {
        lag <- lags[k]

        # SIP1 (EVE = FALSE)
        MA1_sc_power[row_idx, k] <- MA1_sc_power[row_idx, k] +
          as.numeric(SIP.test(seq_data, lag, EVE = FALSE)$p.value <= 0.05)

        # SIP2 (EVE = TRUE)
        MA1_sc_power[row_idx + 1, k] <- MA1_sc_power[row_idx + 1, k] +
          as.numeric(SIP.test(seq_data, lag, EVE = TRUE)$p.value <= 0.05)

        # Oracle (Box-Pierce on noise)
        MA1_sc_power[row_idx + 2, k] <- MA1_sc_power[row_idx + 2, k] +
          as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value <= 0.05)
      }
    }
    # Move index to the next group of 3 methods
    row_idx <- row_idx + 3
  }
}

MA1_sc_power <- MA1_sc_power / 1000
MA1_sc_power
#####################
# # AR(1)
set.seed(2024)

mean_values <- c(5, 10, 20)
n_segs <- c(50, 100, 200)
lags <- c(1, 2, 4, 8)
methods <- c("SIP1", "SIP2", "oracle")

grid <- expand.grid(Method = methods, Seg = n_segs, Mean = mean_values)
grid <- grid[order(grid$Mean, grid$Seg), ]

row_labels <- paste0("Mean=", grid$Mean, "_Seg=", grid$Seg, "_", grid$Method)
AR1_sc_power <- matrix(0, nrow = length(row_labels), ncol = length(lags))
rownames(AR1_sc_power) <- row_labels
colnames(AR1_sc_power) <- paste0("lag", lags)

# Counter for the current row group
row_idx <- 1

for (m_val in mean_values) {
  for (seg in n_segs) {
    cat("Processing: Mean =", m_val, "| Segments =", seg, "...\n")

    mean_structure <- mean_sim(10000, seg, 20, c(-m_val, m_val))
    for (i in 1:n_sims) {
      noise_structure <- acf_sim(10000, model = "AR", ar.coef = 0.1)
      seq_data <- mean_structure + noise_structure

      for (k in seq_along(lags)) {
        lag <- lags[k]

        # SIP1 (EVE = FALSE)
        AR1_sc_power[row_idx, k] <- AR1_sc_power[row_idx, k] +
          as.numeric(SIP.test(seq_data, lag, EVE = FALSE)$p.value <= 0.05)

        # SIP2 (EVE = TRUE)
        AR1_sc_power[row_idx + 1, k] <- AR1_sc_power[row_idx + 1, k] +
          as.numeric(SIP.test(seq_data, lag, EVE = TRUE)$p.value <= 0.05)

        # Oracle (Box-Pierce on noise)
        AR1_sc_power[row_idx + 2, k] <- AR1_sc_power[row_idx + 2, k] +
          as.numeric(Box.test(noise_structure, lag, type = "Box-Pierce")$p.value <= 0.05)
      }
    }
    # Move index to the next group of 3 methods
    row_idx <- row_idx + 3
  }
}

AR1_sc_power <- AR1_sc_power / 1000
AR1_sc_power
##################################### Real Data ######################################################
##################################### Nanopore Data ##################################################
######################################################################################################
library(rhdf5)

file_names <- h5ls('AJO242_4fc438f9_0.fast5')$name
file_names = file_names[grepl("read", file_names)]

# filter by length
len_vec <- as.numeric(NULL)
for (i in 1:length(file_names)) {
  signal = h5read("AJO242_4fc438f9_0.fast5", file_names[i])$Raw$Signal
  len = length(signal)
  len_vec[i] <- len
}

# table 3, acf test on nanopore data
id <- which(len_vec>13760 & len_vec<150000)
test_table <- as.data.frame(matrix(rep(0,8*4), ncol=4))
colnames(test_table) <- c("Reject H0", "Accept H0", "Average p-value", "Maximum p-value")
row.names(test_table) <- paste(paste("s=", rep(c(1,2,4,8), each=2), sep=''), rep(c("SIP1", "SIP2"),4), sep = ",")
lags <- c(1,2,4,8)
for (i in 1:900) {
  signal = h5read("AJO242_4fc438f9_0.fast5", file_names[id[i]])$Raw$Signal
  signal <- signal - mean(signal)
  for (j in 1:4) {
    test1 <- SIP.test(signal, lags[j], EVE=FALSE, warnings = FALSE)
    test_table[2*j-1, 1] <- test_table[2*j-1, 1] + as.numeric(test1$p.value<0.05)
    test_table[2*j-1, 3] <- test_table[2*j-1, 3] + test1$p.value
    test_table[2*j-1, 4] <- max(test_table[2*j-1, 4], test1$p.value)

    test2 <- SIP.test(signal, lags[j], EVE=TRUE, warnings = FALSE)
    test_table[2*j, 1] <- test_table[2*j, 1] + as.numeric(test2$p.value<0.05)
    test_table[2*j, 3] <- test_table[2*j, 3] + test2$p.value
    test_table[2*j, 4] <- max(test_table[2*j, 4], test2$p.value)
  }
}
test_table$`Accept H0` <- 900 - test_table$`Reject H0`
test_table$`Average p-value` <- test_table$`Average p-value`/900
test_table

# save(type_one_error, power_MA1, power_MA4, power_AR1, test_table, file = "results.Rdata")

# # figures
## read id 33
signal33 <- h5read("AJO242_4fc438f9_0.fast5", file_names[33])$Raw$Signal
signal39 <- h5read("AJO242_4fc438f9_0.fast5", file_names[39])$Raw$Signal
# plot(signal, cex=0.5)
# figure 1
png(filename = "Nanopore Sequencing Data.png", width = 1000, height = 800)
par(mar = c(4, 5, 0.5, 1) + 0.1, mfrow=c(2,1), oma = c(0, 0, 0, 0)+0.2)
plot(x = c(10000:15000),y=signal33[10000:15000], cex=0.5, ylab="Signal Level", xlab="Index", main="",
     cex.axis = 1.5, cex.lab =1.5, cex.main=2)
plot(x = c(10000:15000),y=signal39[10000:15000], cex=0.5, ylab="Signal Level", xlab="Index", main="",
     cex.axis = 1.5, cex.lab =1.5, cex.main=2)
dev.off()

# figure 2
png(filename = "ACF Comparison on Nanopore Sequencing Data.png", width = 800, height = 600)
par(cex.axis = 1.5,cex.lab = 1.5, cex.main = 2, mar = c(4, 4.5, 0.5, 1) + 0.1, mfrow=c(2,2))
acf(signal33, main="")
SIP.acf(signal33, lag.max = 4, estimation = FALSE, main="")
acf(signal39, main="")
SIP.acf(signal39, lag.max = 4, estimation = FALSE, main="")
dev.off()

##################################### HadCET Data ##################################################
source("WCMgSA.R")
dat <- read.csv("hadcet.csv") # monthly temperatures
x <- dat$mean
x <- dat$max
x <- dat$min

matplot(dat[, 2:4], type = 'l', xaxt = 'n', xlab = '', ylab = '')
axis(1, at = seq(1, length(x), length.out = 20), label = dat$dates[seq(1, length(x), length.out = 20)])

lags <- c(1, 2, 4, 8)
vars <- c("mean", "max", "min")
temperature_table <- matrix(NA, nrow = 12, ncol = 4)
colnames(temperature_table) <- paste0("lag ", lags)
row_names <- c()
row_idx <- 1
for (v in vars) {
  x <- dat[[v]]

  row_names <- c(row_names, paste0(v, "_SIP1"))
  temperature_table[row_idx, ] <- sapply(lags, function(l) SIP.test(x, EVE = FALSE, m = l)$p.value)
  row_idx <- row_idx + 1

  row_names <- c(row_names, paste0(v, "_SIP2"))
  temperature_table[row_idx, ] <- sapply(lags, function(l) SIP.test(x, EVE = TRUE, m = l)$p.value)
  row_idx <- row_idx + 1

  row_names <- c(row_names, paste0(v, "_Box"))
  temperature_table[row_idx, ] <- sapply(lags, function(l) Box.test(x, type = 'Box-Pierce', lag = l)$p.value)
  row_idx <- row_idx + 1

  w <- wcm.gsa(x, min.len = 10, p.max = 5)

  fhat <- x * 0
  brks <- c(0, w$cp, length(x))

  for (ii in 1:(length(w$cp) + 1)) {
    idx_range <- (brks[ii] + 1):brks[ii + 1]
    fhat[idx_range] <- mean(x[idx_range], na.rm = TRUE)
  }

  res <- x - fhat

  row_names <- c(row_names, paste0(v, "_Box_res"))
  temperature_table[row_idx, ] <- sapply(lags, function(l)
    Box.test(res, type = "Box-Pierce", lag = l)$p.value)
  row_idx <- row_idx + 1
}
rownames(temperature_table) <- row_names
temperature_table

# figure
dev.new(width = 10, height = 14)   # <- key for figure height
par(
  mfrow = c(3, 3),
  mar = c(5, 5, 4, 2),
  cex.main = 2.0,
  cex.lab  = 1.8,
  cex.axis = 1.8
)

for (v in c("mean", "max", "min")) {

  # Column 1: Original Data
  plot(dat$dates, dat[[v]],
       type = "l",
       main = paste(v),
       xlab = "Date",
       ylab = "Value",
       col = "darkblue",
       lwd = 2)

  # Column 2: Default ACF
  acf(dat[[v]],
      main = paste("standard ACF on", v),
      lag.max = 8)

  # Column 3: SIP ACF (inherits par settings)
  SIP.acf(dat[[v]],
          estimation = FALSE,
          lag.max = 8,
          main = paste("shift-immune ACF on", v))
}
dev.off()

######################## pulsar timing ########################################
J1713 <- read.table("J1713_residuals.txt", header = TRUE,
                    comment.char = "#") # no sig
J1643 <- read.table("J1643_red.txt", header = TRUE,
                    comment.char = "#") # sig， with red noise
colnames(J1713) <- c("MJD", "residual_seconds")
colnames(J1643) <- c("MJD", "residual_seconds")

# Define the datasets and lags
datasets <- list("J1713" = J1713$residual_seconds,
                 "J1643" = J1643$residual_seconds)
lags <- c(1, 2, 4, 8)

# Initialize the 6x4 results matrix
pulsar_table <- matrix(NA, nrow = 6, ncol = 4)
colnames(pulsar_table) <- paste0("Lag_", lags)
row_labels <- c()

row_idx <- 1
for (name in names(datasets)) {
  x <- datasets[[name]]

  # 1. SIP.test (EVE = FALSE)
  row_labels <- c(row_labels, paste0(name, "_SIP1"))
  pulsar_table[row_idx, ] <- sapply(lags, function(l) SIP.test(x, EVE = FALSE, m = l)$p.value)
  row_idx <- row_idx + 1

  # 2. SIP.test (EVE = TRUE)
  row_labels <- c(row_labels, paste0(name, "_SIP2"))
  pulsar_table[row_idx, ] <- sapply(lags, function(l) SIP.test(x, EVE = TRUE, m = l)$p.value)
  row_idx <- row_idx + 1

  # 3. Box-Pierce Test
  row_labels <- c(row_labels, paste0(name, "_Box"))
  pulsar_table[row_idx, ] <- sapply(lags, function(l) Box.test(x, type = 'Box-Pierce', lag = l)$p.value)
  row_idx <- row_idx + 1
}

rownames(pulsar_table) <- row_labels
pulsar_table

dev.off()

dev.new(width = 20, height = 14)

# Global graphical settings (affects ALL panels, including SIP.acf)
par(
  mfrow = c(2, 3),
  mar = c(4.5, 4.5, 3.5, 1.5),  # bigger margins
  cex = 0.6,                  # scale EVERYTHING
  cex.main = 2.4,
  cex.lab  = 2.3,
  cex.axis = 2.3
)

# Define datasets
residual_data <- list(
  "J1713+0747" = J1713$residual_seconds,
  "J1643-1224" = J1643$residual_seconds
)

for (name in names(residual_data)) {
  x_res <- residual_data[[name]]

  # Column 1: Residuals
  plot(x_res,
       type = "p",
       pch = 16,
       cex = 0.7,      # slightly larger points
       main = name,
       xlab = "index",
       ylab = "residual_seconds",
       col = "darkblue")

  # Column 2: Default ACF
  acf(x_res,
      lag.max = 20,
      main = paste("standard ACF on", name))

  # Column 3: SIP ACF (inherits par())
  SIP.acf(x_res,
          estimation = FALSE,
          lag.max = 8,
          main = paste("shift-immune ACF on", name))
}

dev.off()
