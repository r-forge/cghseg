library(cghseg)

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

## pour déterminer combien de segments au max par profil
## plus Kmax augmente pour chaque profil plus ca coute cher (l'algo est de complexité (Kmax x n) pour chaque profil)
uniKmax      = getuniKmax(CGHd,CGHo)
## pour déterminer combien de segments au max au total
multiKmax    = getmultiKmax(CGHd,CGHo,uniKmax)


## fonction generale pour la segmentation multivariee
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)

## options :
## select(CGHo)   = "none" ou "mBIC" : selection du nombre de segment. Soit pas de selection (donc Kmax segments) ou critere BIC modifié
## calling(CGHo)  = T/F : utilisation du modèle de segmentation/clustering (T) pour réduire le nombre de niveaux des segments
## wavenorm(CGHo) = "none","spline","position" : méthodes de normalisation de l'effet vague (wavenormalization)


## toutes les fonctions multi... utilisent du uni... qu'on pourrait paralleliser
## l'algo de recherche du minimum du BIC par l'algo des cordes : golden.search.R

