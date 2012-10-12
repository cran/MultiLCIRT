test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10){
	
# test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10)
    J = dim(S)[2]
	out0 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi0,tol=tol)  # unidimensional model
	out1 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi1,tol=tol)  # bidimensional model
	dev = 2*(out1$lk-out0$lk)
	df = out1$np-out0$np
	pv = 1-pchisq(dev,df)
	cat("\n")
	cat(c("Log-likelihood of the constrained model =  ",out0$lk,"\n"))
	cat(c("Log-likelihood of the unconstrained model =",out1$lk,"\n"))
	cat(c("Deviace =                                  ",dev,"\n"))
	cat(c("Degree of freedom =                        ",df,"\n"))
	cat(c("P-value =                                  ",pv,"\n"))
	out = list(out0=out0,out1=out1,dev=dev,df=df,pv=pv)	
}