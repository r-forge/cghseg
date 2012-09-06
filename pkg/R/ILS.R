setMethod(f = "ILS",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax,multiKmax){

            tol          = 1e-2
            select.tmp   = CGHo["select"]
            select(CGHo) = "none"
            command      = parse(text = "invisible(list(mu = mu, theta = B,loglik = loglik,nbiter = iter))")
    
            nbdata   = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}) )
            M        = length(names(.Object@Y))
            n.com    = length(.Object@Y[[1]])
            eps      = Inf
            iter     = 0
			
			if (CGHo@nbprocs>1){
				## Initial data sends, will be reused but not resend
				## Data are emulated to belong to .GlobalEnv
				## since worker function will also belong to .GlobalEnv
				assign("Y.ref", .Object@Y, envir = .GlobalEnv)
				clusterExport(CGHo@cluster, "Y.ref")
				assign("uniKmax.ref", uniKmax, envir = .GlobalEnv)
				clusterExport(CGHo@cluster, "uniKmax.ref")
				assign("CGHo.ref", CGHo, envir = .GlobalEnv)
				clusterExport(CGHo@cluster, "CGHo.ref")
			}
            
            mu       = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
            B        = list(waveffect = rep(0,n.com), GCeffect = rep(0,n.com))
	    mu.tmp   = mu 
	            
            while ( (eps > tol) & (iter < CGHo@itermax)){
              iter                = iter+1
              B                   = getbias(.Object,CGHo,mu,B)		
              removebias(.Object) = B$waveffect+B$GCeffect
              mu                  = multisegmean(.Object,CGHo,uniKmax,multiKmax)$mu
              revertbias(.Object) = B$waveffect+B$GCeffect
              eps = max(sapply(names(.Object@Y),FUN=function(m,x,y){xk = rep(x[[m]]$mean,x[[m]]$end-x[[m]]$begin+1); yk =rep(y[[m]]$mean,y[[m]]$end-y[[m]]$begin+1) ; return(max(abs((xk-yk)/xk)))},mu.tmp,mu))
	      mu.tmp = mu
            } # end while
            
            loglik       = ILS.loglik(.Object,mu,B)
            select(CGHo) = select.tmp
            eval(command)
            
          })


setMethod(f = "ILS.loglik",signature = "CGHdata",
          definition = function(.Object,mu,bias){
            n   = Reduce("sum",lapply(.Object@Y,FUN = function(y){length(y[!is.na(y)])}))
            RSS = lapply(names(.Object@Y),FUN = function(m){
              nk      = mu[[m]]$end -  mu[[m]]$begin + 1
              rss     = sum( (.Object@Y[[m]] - rep(mu[[m]]$mean,nk) - bias$waveffect - bias$GCeffect)^2, na.rm = TRUE)
            })
            RSS    = 0.5* Reduce("sum",RSS)
            loglik = -  (n/2)*(log(2*pi*RSS/n)+1)
            invisible(loglik)
          })


