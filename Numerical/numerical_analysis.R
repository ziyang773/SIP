devtools::install_github("ziyang773/SIP")
library(SIP)
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
      for (i in 2:length(cpts)){
        start <- cpts[i-1]
        end <- cpts[i]
        demeaned_seq[(start+1):(end)] <- seq[(start+1):(end)] - mean(seq[(start+1):end])
      }
      type_one_error[5*k, j] = type_one_error[5*k, j] + as.numeric(Box.test(demeaned_seq, lag, type = "Box-Pierce")$p.value<=0.05)
    }
  }
}
type_one_error <- type_one_error/1000
type_one_error


##################################### Table 2. power comparison ##################################################
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

##################################### Nanopore Data ##################################################
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
test_table <- as.data.frame(matrix(rep(0,6*4), ncol=4))
colnames(test_table) <- c("Reject H0", "Accept H0", "Average p-value", "Maximum p-value")
row.names(test_table) <- paste(paste("s=", rep(c(1,2,4), each=2), sep=''), rep(c("SIP1", "SIP2"),3), sep = ",")
lags <- c(1,2,4)
for (i in 1:900) {
  signal = h5read("AJO242_4fc438f9_0.fast5", file_names[id[i]])$Raw$Signal
  signal <- signal - mean(signal)
  for (j in 1:3) {
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

save(type_one_error, power_MA1, power_MA4, power_AR1, test_table, file = "results.Rdata")

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
