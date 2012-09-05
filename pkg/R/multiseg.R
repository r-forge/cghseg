setMethod(f = "multiseg",signature = "CGHdata",
          definition = function(.Object,CGHo,uniKmax=NULL,multiKmax=NULL){

            uniKmax    = getuniKmax(.Object,CGHo,uniKmax)
            multiKmax  = getmultiKmax(.Object,CGHo,uniKmax,multiKmax)  
            CGHr       = new("CGHresults",CGHd=.Object,CGHo=CGHo)
            from(CGHr) = "multiseg"
            Res        = list()
            
            if ((CGHo["GCnorm"]=="linear") & (is.null(.Object@GCcontent))){
              cat("[check multiseg] GCnorm can not be performed \n")
              cat("[check multiseg] check that GCcontent is a record in the data\n")
              stop()
            }
            
            if (CGHo["select"] != "none"){
              if ( (CGHo["calling"]==FALSE) & (CGHo["wavenorm"]=="none")  & (CGHo@GCnorm=="none")){
                cat("[multiseg] multisegmean running \n")
                Res = multisegmean(.Object,CGHo,uniKmax,multiKmax)
              } else {
				  if (CGHo@nbprocs>1){		
					  CGHo@cluster <- makeCluster(getOption("cl.cores", CGHo@nbprocs))
				  }
				cat("Golden search                           \n")
                Kh        = golden.search(.Object,CGHo,uniKmax,multiKmax)
				cat("Best segmentation                       \n")
                multiKmax = Kh
                eval(fun2run(CGHo))
				if (CGHo@nbprocs>1){
					stopCluster(CGHo@cluster)
				}
              }              
            } else {
              eval(fun2run(CGHo))
            }
            
            if ( (CGHo["wavenorm"]!="none") | (CGHo["GCnorm"]!="none")  ){
              theta(CGHr) =  Res$theta
            }

            if (!is.null(.Object["genomic.position"])){
              x      = .Object["genomic.position"]
              Res$mu = lapply(Res$mu,FUN = function(y){
                       y$begin = x[y$begin]
                       y$end   = x[y$end]
                       return(y)
                     })
            }
            
            mu(CGHr)     = Res$mu 
            loglik(CGHr) = list(multiloglik = Res$loglik)
            nbiter(CGHr) = Res$nbiter
            return(CGHr)
                        
          } # end function
          )
