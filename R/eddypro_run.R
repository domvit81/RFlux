
eddypro_run <- function(siteID, path_eddypro_bin, path_eddypro_projfiles, showLOG=TRUE){

workdir <- getwd()
setwd(path_eddypro_projfiles)
setwd("..")
eddypro.path <- getwd()
dir.create(paste0(getwd(),"/tmp"))

if (Sys.info()[['sysname']]=="Darwin") command_OpSys <- " -s mac "
if (Sys.info()[['sysname']]=="Windows") command_OpSys <- " -s win "
if (Sys.info()[['sysname']]=="Linux") command_OpSys <- " -s linux "

file.copy(from=path_eddypro_bin, to=eddypro.path, recursive=TRUE, overwrite=TRUE)

setwd(paste0(eddypro.path, "/bin"))
getwd()

rp.command <- paste0("./eddypro_rp", command_OpSys, path_eddypro_projfiles, "/", siteID, ".eddypro")
rp <- system(rp.command, intern=!showLOG)

fcc.command <- paste0("./eddypro_fcc", command_OpSys, path_eddypro_projfiles, "/", siteID, ".eddypro")
fcc <- system(fcc.command, intern=!showLOG)

unlink(paste0(eddypro.path, "/bin"), recursive=TRUE)
unlink(paste0(eddypro.path, "/tmp"), recursive=TRUE)
setwd(workdir)
cat("\n The results of EddyPro run are stored in ", path_eddypro_projfiles, "\n", sep="")
}
