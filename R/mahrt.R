
mahrt <- function(x) {
	x <- na.omit(x,na.rm=TRUE)
	ifelse(nrow(x)>20000, WL <- 5995, WL <- 2995)
	COV <- cov(x[,1],x[,2],  use="complete.obs")
	COVs <- rollapply(x, FUN=function(x) cov(x[,1], x[,2],  use="complete.obs"), width=WL, by=WL,by.column=FALSE)
	COVw <- rollapply(x, FUN=function(x) cov(x[,1], x[,2],  use="complete.obs"), width=WL/6, by=WL/6,by.column=FALSE)
	sigmaB <- sqrt(sum((COVs - COV)^2)/(length(COVs)-1))
	ifelse(length(COVw)>=6, sigmaW1 <- sqrt(1/5*sum((COVw[1:6]-COVs[1])^2)), sigmaW1 <- NA)
	ifelse(length(COVw)>=12, sigmaW2 <- sqrt(1/5*sum((COVw[7:12]-COVs[2])^2)), sigmaW2 <- NA)
	ifelse(length(COVw)>=18, sigmaW3 <- sqrt(1/5*sum((COVw[13:18]-COVs[3])^2)), sigmaW3 <- NA)
	ifelse(length(COVw)>=24, sigmaW4 <- sqrt(1/5*sum((COVw[19:24]-COVs[4])^2)), sigmaW4 <- NA)
	ifelse(length(COVw)>=30, sigmaW5 <- sqrt(1/5*sum((COVw[25:30]-COVs[5])^2)), sigmaW5 <- NA)
	ifelse(length(COVw)>30 & length(COVw)<=36, sigmaW6 <- sqrt(1/length(COVw)*sum((COVw[31:length(COVw)]-COVs[6])^2)), sigmaW6 <- NA)
	sigmaWi <- c(sigmaW1, sigmaW2, sigmaW3, sigmaW4, sigmaW5, sigmaW6)
	sigmaW <- mean(sigmaWi, na.rm=TRUE)
	stat <- sigmaB/(sigmaW/sqrt(length(na.omit(sigmaWi))))
	return(list("M98"=stat))
	}
