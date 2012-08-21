# R --no-restore --no-save --args ... <  cghseg.examples.args.R
# with ... = 
# 		profMode parallelMode simul M N writeMode
# 		profMode parallelMode load M N compareMode

args <- commandArgs(trailingOnly = TRUE)
profMode = as.logical(args[1])
nbp = as.integer(args[2])
typeMode = args[3]
iM = as.integer(args[4])
iN = as.integer(args[5])
if (typeMode == "simul"){
	writeMode = as.logical(args[6])	
}
if (typeMode == "load"){
	compareMode = as.logical(args[6])	
}
optimMode = as.logical(args[7])


##############################################################################
#install.packages("/home/vmiele/Work/Developpement/cghseg/cghseg_0.0.1.tar.gz")
#unloadNamespace("cghseg")
library(cghseg)

## OPTIMIZATION environment selection
#is_optimization_mode <- function()
#{
#	exists(".OPTIMIZATION", envir = globalenv()) && 
#			get(".OPTIMIZATION", envir = globalenv())
#}
#set_optimization_mode <- function(on = FALSE)
#{
#	old_value <- is_optimization_mode()
#	.OPTIMIZATION <<- on
#	invisible(old_value)
#}
#set_optimization_mode(optimMode)

## Parallel environment selection
#is_parallel_mode <- function()
#{
#	exists(".PARALLEL", envir = globalenv()) && 
#			get(".PARALLEL", envir = globalenv())
#}
#set_parallel_mode <- function(on = FALSE)
#{
#	old_value <- is_parallel_mode()
#	.PARALLEL <<- on
#	invisible(old_value)
#}
#set_parallel_mode(parMode)
#if (is_parallel_mode()){
#	library(parallel)
#}

##############################################################################
## M : nombre d'individus
## n : nombre de positions
## pour les microarrays tiling n ~ 1e5-1e6
## pour une etude M ~ 100-500
## k.mean : nombre de segments en moyenne par profil
## SNR : Signal to Noise Ratio (+grand +facile)
## lambda : SNR bis (auxiliaire)

#simul        = simulprofiles(M=100,n=1000,k.mean=10,SNR=5,lambda=10)
if (typeMode == "simul"){
	simul        = simulprofiles(M=iM,n=iN,k.mean=10,SNR=5,lambda=10)
	if (isTRUE(writeMode)){
		filename = paste("tests/simul_M",iM,"_N",iN,".RData", sep="")
		save(simul, file=filename)
	}
}
if (typeMode == "load"){
	filename = paste("tests/simul_M",iM,"_N",iN,".RData", sep="")
	load(file=filename)
}

## simul$Y : matrice des (n x M) signaux

##############################################################################
## CGHdata : classe d'objet pour les data
## CGHoptions : classe d'options
## CGHresults : classe pour les resultats

CGHd         = new("CGHdata",Y=simul$Y)
CGHo         = new("CGHoptions")

# Multicore
nbprocs(CGHo) = nbp

## pour d�terminer combien de segments au max par profil
## plus Kmax augmente pour chaque profil plus ca coute cher (l'algo est de complexit� (Kmax x n) pour chaque profil)
if (isTRUE(profMode)){
	Rprof("profiling/getuniKmax.out")
}
uniKmax      = getuniKmax(CGHd,CGHo)
if (isTRUE(profMode)){
	Rprof()
	summaryRprof("profiling/getuniKmax.out")
}
## pour d�terminer combien de segments au max au total
multiKmax    = getmultiKmax(CGHd,CGHo,uniKmax)


##############################################################################
## fonction generale pour la segmentation multivariee
#ptm = proc.time()
if (isTRUE(profMode)){
	Rprof("profiling/multiseg.out")
}
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)
if (isTRUE(profMode)){
	Rprof()
	summaryRprof("profiling/multiseg.out")
}
#proc.time() - ptm

if (typeMode == "simul"){
	if (isTRUE(writeMode)){
		filename = paste("tests/res_M",iM,"_N",iN,".RData", sep="")
		save(CGHr, file=filename)
	}
}
if (typeMode == "load"){
	if (isTRUE(compareMode)){
		thisCGHr = CGHr
		filename = paste("tests/res_M",iM,"_N",iN,".RData", sep="")
		load(file=filename)
		cat("Integration test ",iN,"_",iM," is ",as.numeric(thisCGHr@loglik)," == ",as.numeric(CGHr@loglik)," = ",as.numeric(thisCGHr@loglik) ==  as.numeric(CGHr@loglik),"\n", sep="")
	}
}

## options :
## select(CGHo)   = "none" ou "mBIC" : selection du nombre de segment. Soit pas de selection (donc Kmax segments) ou critere BIC modifi�
## calling(CGHo)  = T/F : utilisation du mod�le de segmentation/clustering (T) pour r�duire le nombre de niveaux des segments
## wavenorm(CGHo) = "none","spline","position" : m�thodes de normalisation de l'effet vague (wavenormalization)


## toutes les fonctions multi... utilisent du uni... qu'on pourrait paralleliser
## l'algo de recherche du minimum du BIC par l'algo des cordes : golden.search.R
