#install.packages("/home/vmiele/Work/Developpement/cghseg/cghseg_0.0.1.tar.gz")
#unloadNamespace("cghseg")
library(cghseg)

is_parallel_mode <- function()
{
	exists(".PARALLEL", envir = globalenv()) && 
			get(".PARALLEL", envir = globalenv())
}
set_parallel_mode <- function(on = FALSE)
{
	old_value <- is_parallel_mode()
	.PARALLEL <<- on
	invisible(old_value)
}
set_parallel_mode(FALSE)

## M : nombre d'individus
## n : nombre de positions
## pour les microarrays tiling n ~ 1e5-1e6
## pour une etude M ~ 100-500
## k.mean : nombre de segments en moyenne par profil
## SNR : Signal to Noise Ratio (+grand +facile)
## lambda : SNR bis (auxiliaire)

simul        = simulprofiles(M=50,n=100,k.mean=10,SNR=5,lambda=10)
## simul$Y : matrice des (n x M) signaux

## CGHdata : classe d'objet pour les data
## CGHoptions : classe d'options
## CGHresults : classe pour les resultats

CGHd         = new("CGHdata",Y=simul$Y)
CGHo         = new("CGHoptions")

## pour d�terminer combien de segments au max par profil
## plus Kmax augmente pour chaque profil plus ca coute cher (l'algo est de complexit� (Kmax x n) pour chaque profil)
uniKmax      = getuniKmax(CGHd,CGHo)
## pour d�terminer combien de segments au max au total
multiKmax    = getmultiKmax(CGHd,CGHo,uniKmax)


## fonction generale pour la segmentation multivariee
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)

## options :
## select(CGHo)   = "none" ou "mBIC" : selection du nombre de segment. Soit pas de selection (donc Kmax segments) ou critere BIC modifi�
## calling(CGHo)  = T/F : utilisation du mod�le de segmentation/clustering (T) pour r�duire le nombre de niveaux des segments
## wavenorm(CGHo) = "none","spline","position" : m�thodes de normalisation de l'effet vague (wavenormalization)


## toutes les fonctions multi... utilisent du uni... qu'on pourrait paralleliser
## l'algo de recherche du minimum du BIC par l'algo des cordes : golden.search.R
