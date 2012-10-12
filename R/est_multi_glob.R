est_multi_glob<-function(yv,X,model,ind=rep(1,length(yv)),w=rep(1,length(yv)),be=NULL){

#        [be,lk,P] = est_multi_glob(yv,X,model,ind=rep(1,length(yv)),w=rep(1,length(yv)),be=NULL)
#
# Fit multinomal logit and proportional odds models
#
# INPUT:
# y:     vector of responses (n x 1)
# X:     matrix of covariates (l-1 x nc x ni)
# model: "m" = multinomial with first reference cateogory
#        "g" = version based on global logits 
#        "l" = version based on local logits 
# ind:   list of indices of X for each observation in y
# w:     vector of weights
# be:    initial parameter vector if available
#
# OUTPUT:
# be:  parameters estimates
# lk:  log-likelihood
# P:   matrix of probabilities

# Preliminaries
	l = nrow(X)+1
	ncov = ncol(X)
	n = length(yv)
	Y = matrix(0,n,l)
	for(y in 0:(l-1)) Y[yv==y,y+1] = 1
	if(model=="g" || model=="l"){
		out = matr_glob(l,model)
		Co = out$Co; Ma = out$Ma 
	}
	G = rbind(matrix(0,1,l-1),diag(l-1))
	nd = max(ind);
# Initial estimate
	if(is.null(be)){
		dist = colMeans(Y)
		if(model=="m") eta = log(dist[2:end]/dist[1])
		if(model=="g" || model=="l") eta = Co%*%log(Ma%*%dist)
		Xm = 0
		for(h in 1:nd) Xm = Xm+X[,,h]*sum(w[ind==h])
		Xm = matrix(Xm/sum(w),l-1,ncov)
		be = ginv(t(Xm)%*%Xm)%*%t(Xm)%*%eta
  	}
# Compute log-likelihood
	P = matrix(0,nd,l)
	for(h in 1:nd){
		Xh = matrix(X[,,h],l-1,ncov)
		if(model=="m"){
			P[h,] = exp(G%*%(Xh%*%be))
			P[h,] = P[h,]/sum(P[h,])
		}
		if(model=="g" | model=="l") P[h,] = inv_glob(Xh%*%be,model)$p
	}
	P = pmax(P,10^-10)
	pm = rep(0,n)
	for(h in 1:nd) for(i in which(ind==h)) pm[i] = P[h,yv[i]+1]
	lk = sum(w*log(pm))
# Iterate until convergence
	lko = lk; it = 0
	while((abs(lk-lko)>10^-5 | it==0) & it<10^4){
		it = it+1
		s = 0; FI = 0
		for(h in 1:nd){
			indh = which(ind==h)
			p = P[h,]
			Xh = matrix(X[,,h],l-1,ncov)
			if(model=="m") D = t(X[,,h])%*%t(G)
			if(model=="g" || model=="l"){
				ve = 1/as.vector(Ma%*%p)
				D = solve(Co%*%diag(ve)%*%Ma%*%diag(p)%*%G)
				D = t(Xh)%*%(t(D)%*%t(G))
			}
      		for(i in indh) s = s+w[i]*D%*%(Y[i,]-p)
      		FI = FI+sum(w[indh])*(D%*%(diag(p)-p%*%t(p))%*%t(D))
		}
		if(rcond(FI)<10^-15){
			cat("matrix close to singularity in est_multi_glog\n")
			dbe = ginv(FI)%*%s
			}else dbe = solve(FI)%*%s
		mdbe = max(abs(dbe))
		if(mdbe>0.25) dbe = dbe/mdbe*0.25
		be = be+dbe
# Compute log-likelihood
    		P = matrix(0,nd,l)
    		for(h in 1:nd){
    			if(model=="m") P[h,] = exp(G%*%(X[,,h]%*%be)); P[h,] = P[h,]/sum(P[h,])
      		if(model=="g" || model=="l") P[h,] = inv_glob(X[,,h]%*%be,model)$p
      	}
	    P = pmax(P,10^-15)
	    pm = rep(0,n)
	    for(h in 1:nd) for(i in which(ind==h)) pm[i] = P[h,yv[i]+1]
	    lko = lk; lk = sum(w*log(pm))
	}
# output
	out = list(be=be,lk=lk,P=P)	
}