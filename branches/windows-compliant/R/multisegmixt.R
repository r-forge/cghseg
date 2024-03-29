setMethod(f = "multisegmixt",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax,phi){

            P            = length(phi)/3
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            multiKselect = multiKmax

######   Individual segmentations for patients 
			if (CGHo@nbprocs>1){		
				#cat("multisegmixt //                         \r")	
				if (Sys.info()["sysname"] == "Windows"){ 				
					unisegmixt.proxy <- function(m){
						Y.ref	 	= get("Y.ref", envir = .GlobalEnv)
						CGHo.ref	= get("CGHo.ref", envir = .GlobalEnv)
						uniKmax.ref	= get("uniKmax.ref", envir = .GlobalEnv)
						phi.ref		= get("phi.ref", envir = .GlobalEnv)
						n     = length(which(!is.na(Y.ref[[m]])))
						Kmax  = uniKmax.ref[[m]]
						out   = unisegmixt(Y.ref[[m]],CGHo.ref,Kmax,phi.ref)
						J.est = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
						invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
					}
					environment(unisegmixt.proxy) <- .GlobalEnv
					clusterExport(CGHo@cluster, "unisegmixt")	# to be know in unisegmixt.proxy
					assign("phi.ref", phi, envir = .GlobalEnv)
					clusterExport(CGHo@cluster, "phi.ref")
					Y.ref	 	= get("Y.ref", envir = .GlobalEnv)			
					Res = parLapply(CGHo@cluster, names(Y.ref), fun = unisegmixt.proxy)
					names(Res) = names(.Object@Y)
				}
				else{				
					Res = mclapply(names(.Object@Y), FUN = function(m){
								n     = length(which(!is.na(.Object@Y[[m]])))
								Kmax  = uniKmax[[m]]
								out   = unisegmixt(.Object@Y[[m]],CGHo,Kmax,phi)
								J.est = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
								invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
							}, mc.cores = CGHo@nbprocs)
					names(Res) = names(.Object@Y)  
				}
			}
			else{	
				#cat("multisegmixt                            \r")	
				Res = lapply(names(.Object@Y), FUN = function(m){
							n     = length(which(!is.na(.Object@Y[[m]])))
							Kmax  = uniKmax[[m]]
							out   = unisegmixt(.Object@Y[[m]],CGHo,Kmax,phi)
							J.est = n*exp(-((2/n)*out$loglik+log(2*pi)+1))
							invisible(list(t.est = out$t.est, loglik = out$loglik,J.est=J.est))
						})
				names(Res) = names(.Object@Y)  
			}
  
######   Segment Repartition segments across patients 
			#cat("multisegmixt finishing                  \r")
      
            J.est              = lapply(Res,FUN = function(x){x$J.est})
            nbdata             = lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) 
            nbdata             = sum(unlist(nbdata))    
            out.ibp            = segibp(data.frame(J.est = unlist(J.est)),unlist(uniKmax),multiKmax)
            multiloglik        = -(nbdata/2)*(log(2*pi*out.ibp[[1]]/nbdata)+1)
            seg.rep            = out.ibp[[2]]
            row.names(seg.rep) = names(.Object@Y)
    
######   Outputs  

            multiKselect    = multiKmax  
            mu              = multisegout(.Object,seg.rep,Res,multiKselect)
            select(CGHo)    = select.tmp
			#cat("multisegmixt finished                   \r")             
            invisible(list(mu=mu,loglik = multiloglik[length(multiloglik)],nbiter=0)) 
          }
          )
