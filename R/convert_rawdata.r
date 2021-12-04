convert_rawdata <- function(
	file_path_in,
	file_path_out,
	info_U=list(NULL,NULL,1,0),   
	info_V=list(NULL,NULL,1,0),
	info_W=list(NULL,NULL,1,0),
	info_T_SONIC=list(NULL,NULL,1,0),
	info_CO2=list(NULL,NULL,1,0),
	info_H2O=list(NULL,NULL,1,0),
	info_T_CELL=list(NULL,NULL,1,0),
	info_T_CELL_IN=list(NULL,NULL,1,0),
	info_T_CELL_OUT=list(NULL,NULL,1,0),
	info_PRESS_CELL=list(NULL,NULL,1,0),
	info_SA_DIAG=NULL,
	info_GA_DIAG=NULL,
	na.strings, 
	nrow.header,
	sep,
	dec,
	timestamp_loc=c(NULL, NULL),
	timestamp_format,
	siteID="Cc-Xxx",
	info_LN="LLL",
	info_FN="FNN")
	{
	
	molw_co2 <- 44.01 ## g/mol
	molw_h2o <- 18.02 ## g/mol
	molw_ch4 <- 16.04 ## g/mol
	molw_n2o <- 44.01 ## g/mol
	
	timestamp_orig <- str_sub(file_path_in, start=timestamp_loc[1], end=timestamp_loc[2])
	data_orig <- fread(file_path_in, na.strings=na.strings, data.table=FALSE, skip=nrow.header, header=FALSE, sep=sep, dec=dec)
	
	#Wind Components - outputs in m s-1
	if(!is.null(info_U[[1]])){
	if(info_U[[2]]=="mm s-1") U <- (data_orig[,info_U[[1]]]*info_U[3]+info_U[[4]])/10^3
	if(info_U[[2]]=="cm s-1") U <- (data_orig[,info_U[[1]]]*info_U[3]+info_U[[4]])/10^2;
	if(info_U[[2]]=="m s-1")  U <- data_orig[,info_U[[1]]]*info_U[[3]]+info_U[[4]];
	}
	if(is.null(info_U[[1]])) U <- rep(NA, nrow(data_orig))

	if(!is.null(info_V[[1]])){
	if(info_V[[2]]=="mm s-1") V <- (data_orig[,info_V[[1]]]*info_V[[3]]+info_V[[4]])/10^3
	if(info_V[[2]]=="cm s-1") V <- (data_orig[,info_V[[1]]]*info_V[[3]]+info_V[[4]])/10^2;
	if(info_V[[2]]=="m s-1")  V <- data_orig[,info_V[[1]]]*info_V[[3]]+info_V[[4]];
	}
	if(is.null(info_V[[1]])) V <- rep(NA, nrow(data_orig))

	if(!is.null(info_W[[1]])){
	if(info_W[[2]]=="mm s-1") W <- (data_orig[,info_W[[1]]]*info_W[[3]]+info_W[[4]])/10^3
	if(info_W[[2]]=="cm s-1") W <- (data_orig[,info_W[[1]]]*info_W[[3]]+info_W[[4]])/10^2;
	if(info_W[[2]]=="m s-1")  W <- data_orig[,info_W[[1]]]*info_W[[3]]+info_W[[4]];
	}
	if(is.null(info_W[[1]])) W <- rep(NA, nrow(data_orig))

	#T_SONIC - output in degree Kelvin
	if(!is.null(info_T_SONIC[[1]])){
	if(info_T_SONIC[[2]]=="K")  T_SONIC <- data_orig[,info_T_SONIC[[1]]]*info_T_SONIC[[3]]+info_T_SONIC[[4]]
	if(info_T_SONIC[[2]]=="C") T_SONIC <- (data_orig[,info_T_SONIC[[1]]]*info_T_SONIC[[3]]+info_T_SONIC[[4]])+273.15
	}
	if(is.null(info_T_SONIC[[1]])) T_SONIC <- rep(NA, nrow(data_orig))

	##CO2 molar density in mmol m-3, CO2 molar fraction in umol/mol (ppm), CO2 mixing ratio in umol/mol (ppm)
	if(!is.null(info_CO2[[1]])){
	if(info_CO2[[2]]=="mmol m-3" | info_CO2[[2]]=="ppm") CO2 <- data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]];
	if(info_CO2[[2]]=="ppt") CO2 <- (data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]])*10^3;
	if(info_CO2[[2]]=="umol m-3" | info_CO2[[2]]=="ppb") CO2 <- (data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]])/10^3;
	if(info_CO2[[2]]=="g m-3") CO2 <- 24.45*(data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]])*10^3/molw_co2;
	if(info_CO2[[2]]=="mg m-3") CO2 <- 24.45*(data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]])/molw_co2;
	if(info_CO2[[2]]=="ug m-3") CO2 <- 24.45*(data_orig[,info_CO2[[1]]]*info_CO2[[3]]+info_CO2[[4]])*10^-3/molw_co2
	}
	if(is.null(info_CO2[[1]])) CO2 <- rep(NA, nrow(data_orig))

	##H2O molar density in mmol m-3, H2O molar fraction in mmol/mol (ppt), H2O mixing ratio in mmol/mol (ppt)
	if(!is.null(info_H2O[[1]])){
	if(info_H2O[[2]]=="mmol m-3" | info_H2O[[2]]=="ppt") H2O <- data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]];
	if(info_H2O[[2]]=="umol m-3" | info_H2O[[2]]=="ppm") H2O <- (data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]])/10^3;
	if(info_H2O[[2]]=="ppb") H2O <- (data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]])/10^6;
	if(info_H2O[[2]]=="g m-3") H2O <- 24.45*(data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]])*10^3/molw_h2o;
	if(info_H2O[[2]]=="mg m-3") H2O <- 24.45*(data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]])/molw_h2o;
	if(info_H2O[[2]]=="ug m-3") H2O <- 24.45*(data_orig[,info_H2O[[1]]]*info_H2O[[3]]+info_H2O[[4]])*10^-3/molw_h2o
	}
	if(is.null(info_H2O[[1]])) H2O <- rep(NA, nrow(data_orig))


	#SA_DIAG - dimensionless
	if(!is.null(info_SA_DIAG)) SA_DIAG <- data_orig[,info_SA_DIAG[[1]]]
	if(is.null(info_SA_DIAG)) SA_DIAG <- rep(NA, nrow(data_orig))

	#GA_DIAG - dimensionless
	if(!is.null(info_GA_DIAG)) GA_DIAG <- data_orig[,info_GA_DIAG[[1]]]
	if(is.null(info_GA_DIAG)) GA_DIAG <- rep(NA, nrow(data_orig))


	#T_CELL - output in degree Celsius
	if(!is.null(info_T_CELL[[1]])){
	if(info_T_CELL[[2]]=="C") T_CELL <- data_orig[,info_T_CELL[[1]]]*info_T_CELL[[3]]+info_T_CELL[[4]];
	if(info_T_CELL[[2]]=="K")  T_CELL <- (data_orig[,info_T_CELL[[1]]]*info_T_CELL[[3]]+info_T_CELL[[4]])-273.15
	}
	if(is.null(info_T_CELL[[1]])) T_CELL <- rep(NA, nrow(data_orig))

	if(!is.null(info_T_CELL_IN[[1]])){
	if(info_T_CELL_IN[[2]]=="C")  T_CELL_IN <- data_orig[,info_T_CELL_IN[[1]]]*info_T_CELL_IN[[3]]+info_T_CELL_IN[[4]];
	if(info_T_CELL_IN[[2]]=="K") T_CELL_IN <- (data_orig[,info_T_CELL_IN[[1]]]*info_T_CELL_IN[[3]]+info_T_CELL_IN[[4]])-273.15
	}
	if(is.null(info_T_CELL_IN[[1]])) T_CELL_IN <- rep(NA, nrow(data_orig))
	
	if(!is.null(info_T_CELL_OUT[[1]])){
	if(info_T_CELL_OUT[[2]]=="C")  T_CELL_OUT <- data_orig[,info_T_CELL_OUT[[1]]]*info_T_CELL_OUT[[3]]+info_T_CELL_OUT[[4]];
	if(info_T_CELL_OUT[[2]]=="K") T_CELL_OUT <- (data_orig[,info_T_CELL_OUT[[1]]]*info_T_CELL_OUT[[3]]+info_T_CELL_OUT[[4]])-273.15
	}
	if(is.null(info_T_CELL_OUT[[1]])) T_CELL_OUT <- rep(NA, nrow(data_orig))
	
	#PRESS_CELL - output in kPa
	if(!is.null(info_PRESS_CELL[[1]])){
	if(info_PRESS_CELL[[2]]=="KPa") PRESS_CELL <- data_orig[,info_PRESS_CELL[[1]]]*info_PRESS_CELL[[3]]+info_PRESS_CELL[[4]];
	if(info_PRESS_CELL[[2]]=="Pa")  PRESS_CELL <- (data_orig[,info_PRESS_CELL[[1]]]*info_PRESS_CELL[[3]]+info_PRESS_CELL[[4]])/10^3;
	if(info_PRESS_CELL[[2]]=="hPa")  PRESS_CELL <- (data_orig[,info_PRESS_CELL[[1]]]*info_PRESS_CELL[[3]]+info_PRESS_CELL[[4]])/10
	}
	if(is.null(info_PRESS_CELL[[1]])) PRESS_CELL <- rep(NA, nrow(data_orig))

	data_2exp <- data.frame(U, V, W, T_SONIC, CO2, H2O, SA_DIAG, GA_DIAG, T_CELL, T_CELL_IN, T_CELL_OUT, PRESS_CELL)

	fwrite(data_2exp, paste0(file_path_out, "/", siteID, "_EC_", format(strptime(timestamp_orig, format=timestamp_format, tz="GMT"), format="%Y%m%d%H%M", tz="GMT"),"_", info_LN, "_", info_FN,".csv" ), col.names=TRUE, row.names=FALSE, quote=FALSE, na="-9999",dec=".", sep=",")
	
}

