library(cghseg)

## M : nombre d'individus
## n : nombre de positions
## pour les microarrays tiling n ~ 1e5-1e6
## pour une etude M ~ 100-500
## k.mean : nombre de segments en moyenne par profil
## SNR : Signal to Noise Ratio (+grand +facile)
## lambda : SNR bis (auxiliaire)

simul        = simulprofiles(M=50,n=1000,k.mean=10,SNR=1,lambda=10)
## simul$Y : matrice des (n x M) signaux

## CGHdata : classe d'objet pour les data
## CGHoptions : classe d'options
## CGHresults : classe pour les resultats

CGHd         = new("CGHdata",Y=simul$Y)
CGHo         = new("CGHoptions")
nbprocs(CGHo) = 4

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

## others
calling(CGHo)  = T
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)

calling(CGHo)  = F
wavenorm(CGHo) = "spline"
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)

calling(CGHo)  = T
wavenorm(CGHo) = "spline"
CGHr         = multiseg(CGHd,CGHo,uniKmax,multiKmax)



## Individual segmentations for patients 
            
Res = lapply(names(.Object@Y), FUN = function(m){
  n     = length(which(!is.na(.Object@Y[[m]])))
  Kmax  = uniKmax[[m]]
  out   = unisegmixt(.Object@Y[[m]],CGHo,Kmax,phi)
  J.est = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
  invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
}) 
names(Res) = names(.Object@Y)


Res = lapply(.Object@Y, FUN = function(y,K){
  cat(names(y),"-",K,"\n")
},uniKmax)




names(Res) = names(.Object@Y)


library(cghseg)
n = 100
x = rnorm(n,0,1)
K    = 4
rupt = matrix(Inf,K,2)
rupt[,1] =c(1,5,50,80)
rupt[,2] =c(4,49,79,100)
phi = c(0,2,1,1,0.5,0.5)
nk = rupt[,2]-rupt[,1]+1
xk = sapply(1:K,FUN=function(k){mean(x[rupt[k,1]:rupt[k,2]])})
x2k = sapply(1:K,FUN=function(k){mean(x[rupt[k,1]:rupt[k,2]]^2)})
EMinit(x,rupt,P=2,vh=TRUE)
compactEMinit(xk,x2k,nk,P=2,vh=TRUE)
compactEMalgo(xk,x2k,phi,nk,P=2,vh=TRUE)

