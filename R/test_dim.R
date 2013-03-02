test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10,disp=FALSE){
	
# disp = to display log-likelihood evolution step by step

# test_dim<-function(S,yv,k,link=1,disc=0,difl=0,multi0=1:J,multi1,tol=10^-10,disp=FALSE)
    J = dim(S)[2]
	out0 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi0,tol=tol,disp=disp)  # unidimensional model
	out1 = est_multi_poly(S,yv,k,link=link,disc=disc,difl=difl,multi=multi1,tol=tol,disp=disp)  # bidimensional model
	dev = 2*(out1$lk-out0$lk)
	df = out1$np-out0$np
	pv = 1-pchisq(dev,df)
    table = c(round(out0$lk,4),round(out0$aic,4),round(out0$bic,4),round(out0$np,4),
              round(out1$lk,4),round(out1$aic,4),round(out1$bic,4),round(out1$np,4),
              round(dev,4),df,round(pv,4))
    table = matrix(table,11,1)
    colnames(table) = ""
    rownames(table) = c("Log-likelihood of the constrained model","AIC of the constrained model",
                        "BIC of the constrained model","N.parameters of the constrained model",
                        "Log-likelihood of the unconstrained model","AIC of the unconstrained model",
                        "BIC of the unconstrained model","N.parameters of the unconstrained model",
                        "Deviance","Degrees of freedom","p-value")
    table = as.table(table)
	print(table)	
	out = list(out0=out0,out1=out1,dev=dev,df=df,pv=pv,table=table)	
}