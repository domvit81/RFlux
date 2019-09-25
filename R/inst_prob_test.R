
inst_prob_test <- function(x){

flucts <- x - mean(x, na.rm=TRUE)
sigma_f <- max(0.01, Qn(na.omit(flucts)))
n_spike1 <- length(which(abs(flucts)>5*sigma_f))
n_spike2 <- length(which(abs(flucts)>10*sigma_f))

d0 <- diff(x)
ind <- which(abs(d0) < 10^-3)
x_r <- replace(x, ind+1, NA)
d1 <- diff(x_r)

K_VM97 <- as.numeric(3+kurtosis(detrend(na.omit(x)), na.rm=TRUE))
S_VM97 <- as.numeric(skewness(detrend(na.omit(x)), na.rm=TRUE))
KID0 <- as.numeric(3+kurtosis(d0, na.rm=TRUE))
KID1 <- as.numeric(3+kurtosis(diff(na.omit(x_r)), na.rm=TRUE))

ifelse((length(na.omit(d1)) > 1000),{sigma_d <- max(0.01, Qn(na.omit(d1))); n_spike3 <- length(which(abs(d0)>5*sigma_d)); n_spike4 <- length(which(abs(d0)>10*sigma_d))}, {n_spike3 <- NA; n_spike4 <- NA})

return(list("Skew"= S_VM97, "Kurt"=K_VM97, "KID0"=KID0, "KID1"=KID1, "HF4"=n_spike1, "HF1"=n_spike2, "HD4"=n_spike3, "HD1"=n_spike4))
}
