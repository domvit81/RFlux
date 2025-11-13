
inst_prob_test <- function(x, detrend=FALSE, hz=c(10,20), plot=FALSE, var_name=c("U","V","W","T_SONIC","CO2","H2O"), cex.test=1){

ifelse(length(x)>20000, nl <- 36000, nl <- 18000)

ifelse(detrend==FALSE, flucts <- x - mean(x, na.rm=TRUE), flucts <- detrend(as.vector(na.omit(na.approx(x))))) 
sigma_f <- max(0.01, Qn(na.omit(flucts)))
HF5 <- length(which(abs(flucts)>5*sigma_f))
HF10 <- length(which(abs(flucts)>10*sigma_f))

d0 <- diff(x)
ind <- which(abs(d0) < 10^-3)
x_r <- replace(x, ind+1, NA)
d1 <- diff(x_r)
ifelse((length(na.omit(d1)) > 1000),{sigma_d <- max(0.01, Qn(na.omit(d1))); HD5 <- length(which(abs(d0)>5*sigma_d)); HD10 <- length(which(abs(d0)>10*sigma_d))}, {HD5 <- NA; HD10 <- NA})

K_VM97 <- as.numeric(3+kurtosis(detrend(as.vector(na.omit(x))), na.rm=TRUE))
S_VM97 <- as.numeric(skewness(detrend(as.vector(na.omit(x))), na.rm=TRUE))
KID0 <- as.numeric(3+kurtosis(d0, na.rm=TRUE))
KID <- as.numeric(3+kurtosis(diff(na.omit(x_r)), na.rm=TRUE))

ifelse(is.numeric(x) & min(x, na.rm=TRUE)!=max(x, na.rm=TRUE), AL1 <- acf(x, na.action=na.pass, plot=FALSE)$acf[2], AL1 <- NA)

ifelse(is.numeric(x) & min(x, na.rm=TRUE)!=max(x, na.rm=TRUE), DDI <- max(hist(x, breaks="FD", plot=FALSE)$counts), DDI <- nl) 

ifelse(is.numeric(flucts) & min(flucts, na.rm=TRUE)!=max(flucts, na.rm=TRUE), DIP <- dip.test(flucts)$p.val, DIP <- NA)


if (plot==TRUE){
 	mod_thr <- c(NA, NA, NA, 30, hz*60*30*0.04/2, hz*60*30*0.01/2, hz*60*30*0.04/2, hz*60*30*0.01/2, 0.75, hz*60*2.5, 0.05)
 	sev_thr <- c(NA, NA, NA, 50, hz*60*30*0.04, hz*60*30*0.01, hz*60*30*0.04, hz*60*30*0.01, 0.5, hz*60*5, 0.01)
 	if (var_name=="U" | var_name=="V" | var_name=="W") units <- expression(W~(m*s^-1))
 	if (var_name=="T_SONIC") units <- "T SONIC (K)"
 	if (var_name=="CO2") units <- expression(CO[2]~(mu*mol/mol))
 	if (var_name=="H2O") units <- expression(H[2]*O~(mmol/mol))
 	if (var_name=="U" | var_name=="V" | var_name=="W") dunits <- expression(Delta*W~(m*s^-1))
 	if (var_name=="T_SONIC") dunits <- expression(Delta*T~SONIC~(K))
 	if (var_name=="CO2") dunits <- expression(Delta*CO[2]~(mu*mol/mol))
 	if (var_name=="H2O") dunits <- expression(Delta*H[2]*O~(mmol/mol))
  	nf <- layout(matrix(rbind(c(1,1,1,1,2), c(3,3,3,3,4), c(5,5,5,5,6), c(7,7,7,7,8)), nrow=4, ncol=5))
 	#layout.show(nf)
 	par(mar=c(4,1,1,1), oma=c(0,6,3,1), las=1, cex.axis=1.25, cex.lab=1.25)
  	plot(c(1:length(x))/(hz*60), x, type="l", main="", xlab="", ylab="")
 	mtext(side=1, "Minutes", line=2.5, las=0, cex=.75)	
 	mtext(side=2, units, line=3.5, las=0, cex=.75)
 	box(lwd=2)
 	par(mar=rep(0,4))
 	plot(1:10, 1:10, type="n", axes=FALSE, xlab="", ylab="")
 	mtext(side=3, line=1.5, "SevEr", col="red")
 	mtext(side=3, line=0, "ModEr", col="orange")
 	mtext(side=3, line=-1.5, "NoEr", col="gray")
  	legend("left", 
 	legend=paste0(
 	c("HF5", "HF10"), ": ", 
 	as.character(round(c(HF5, HF10),0))),
 	text.col=c(
 	ifelse(HF5>sev_thr[5],"red",ifelse(HF5>mod_thr[5] & HF5<sev_thr[5], "orange","gray48")),
 	ifelse(HF10>sev_thr[6],"red",ifelse(HF10>mod_thr[6] & HF10<sev_thr[6], "orange","gray48"))),
 	bty="n", cex=cex.test)
 
 	par(mar=c(4,1,1,1))
 	plot(c(1:length(diff(x)))/(hz*60), diff(x), type="l", main="", xlab="", ylab="")
 	mtext(side=1, "Minutes", line=2.5, las=0, cex=.75)	
	mtext(side=2, dunits, line=3.5, las=0, cex=.75)
 	box(lwd=2)
 	par(mar=rep(0,4))
 	plot(1:10, 1:10, type="n", axes=FALSE, xlab="", ylab="")
 	legend("left",
 	legend=paste0(c("KID", "HD5", "HD10"), ": ",
 	as.character(round(c(KID, HD5, HD10),0))),
 	text.col=c(
 	ifelse(KID>sev_thr[4],"red",ifelse(KID>mod_thr[4] & KID<sev_thr[4], "orange","gray48")),
 	ifelse(HD5>sev_thr[7],"red",ifelse(HD5>mod_thr[7] & HD5<sev_thr[7], "orange","gray48")),
 	ifelse(HD10>sev_thr[8],"red",ifelse(HD10>mod_thr[8] & HD10<sev_thr[8], "orange","gray48"))),
 	bty="n", cex=cex.test)
 	
 	par(mar=c(4,1,1,1))
 	Acf(x, lag.max=100, na.action=na.contiguous, plot=TRUE, ci=0, main="", xlab="", ylab="")
	mtext(side=1, "Lags", line=2.5, las=0, cex=.75)	
	mtext(side=2, "ACF", line=3.5, las=0, cex=.75)
  	box(lwd=2)
 	par(mar=rep(0,4))
 	plot(1:10, 1:10, type="n", axes=FALSE, xlab="", ylab="")
 	legend("left",
 	legend=paste0("AL1: ",
 	as.character(round(AL1,3))),
 	text.col=ifelse(AL1<sev_thr[9],"red",ifelse(AL1<mod_thr[9], "orange","gray48")), bty="n", cex=cex.test)
 	
 	par(mar=c(4,1,1,1))
 	hist(x, breaks="FD", main="", xlab="", ylab="")
	mtext(side=1,units, line=2.5, las=0, cex=.75)	
	mtext(side=2, "Frequency", line=3.5, las=0, cex=.75)
	box(lwd=2)
 	par(mar=rep(0,4))
 	plot(1:10, 1:10, type="n", axes=FALSE, xlab="", ylab="")
 	legend("left",
 	legend=paste0(c("DDI", "DIP"), ": ",
 	as.character(round(c(DDI, DIP),3))),
 	text.col=c(
 	ifelse(DDI>sev_thr[10],"red",ifelse(DDI>mod_thr[10] & DDI<sev_thr[10], "orange","gray48")),
 	ifelse(DIP<sev_thr[11],"red",ifelse(DIP<mod_thr[11] & DIP>sev_thr[11], "orange","gray48"))),
 	bty="n", cex=cex.test)
 	}

return(list("Skew"= S_VM97, "Kurt"=K_VM97, "KID0"=KID0, "KID"=KID, "HF5"=HF5, "HF10"=HF10, "HD5"=HD5, "HD10"=HD10, "AL1"=AL1, "DDI"=DDI, "DIP"=DIP))
}
