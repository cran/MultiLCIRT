est_multi_poly <- function(S,yv=rep(1,ns),k,X=NULL,start=0,link=0,disc=0,difl=0,
                           multi=1:J,piv,Th,bec,gac,fort=FALSE,tol=10^-10){

#        [piv,Th,Bec,gac,fv,Phi,Pp,lk,np,aic,bic] = est_multi_poly(S,yv,k,start,link,disc,difl,multi,piv,Th,bec,gac)
#
# Fit Latent Class model and some restricted versions with k classes for ordinal (NA for missing data)
# 
# S    : matrix of available configurations
# X    : matrix of corresponding covariates affecting the ability
# yv   : frequencies of the available configurations
# k    : number of latent classes
# X    : matrix of covariates for the multinomial logit on the class weights
# start: type of startgine values (0 = deterministic, 1 = random)
# link : type of link (0 = LC, 1 = GRM, 2 = PCM)
# disc : discriminating indices (0 = constrained, 1 = free)
# difl : difficulty levels (0 = free, 1 = additive decomposition)
# lk   : maximum log-likelihood
# piv  : weights of the latent classes
# Phi  : conditional distributions given the latent classes
# np   : number of free parameters
# bic  : Bayesian information criterion
# Th,be,ga : parameters for the model 4 (Th=Psi if start==3)
# fv   : list of items constrained
# fort : T for using fortran code for covariates, F using R code only
# tol  : relative tolerance level for convergence

# With k=1
	cat("*-------------------------------------------------------------------------------*\n")
	if(k==1){
	  cat("fit only for LC model with no other imput\n")
	  X = NULL
	}
# Preliminaries
# check problems with standard errors
    cov = !is.null(X)
    if(cov) X = as.matrix(X)
    miss = any(is.na(S))
	ns = nrow(S); J = ncol(S)
    if(miss){
    	cat("Missing data in the dataset, units and items without responses are removed\n")
    	ind = which(apply(is.na(S),1,all))
    	if(length(ind)>0){
        	S = S[-ind,]; yv = yv[-ind]
        	if(!is.null(X)) X = as.matrix(X[-ind,])
        	ind = which(apply(is.na(S),2,all))
        	if(length(ind)>0){
        	    S = S[,-ind]
	            miss = any(is.na(S))
	        }
	    }
    }
    if(miss){R=1*(!is.na(S)); S[is.na(S)]=0}
	l = max(S)+1
	ns = nrow(S); J = ncol(S)
	n = sum(yv)
# checks about the covariates
    if(cov){
    	ncov = ncol(X)
    	out = aggr_data(X)
    	Xdis = out$data_dis; Xlabel = out$label; Xndis = max(out$label)
    	XXdis = array(0,c(k,(k-1)*(ncov+1),Xndis))
    	if(k==2) II = as.matrix(c(0,1)) else {II = diag(k); II = II[,-1]}
    	for(i in 1:Xndis) XXdis[,,i] = II%x%t(c(1,Xdis[i,]))
    }
# about models
	if(link==1) ltype = "g" else if(link==2) ltype = "l"
	if(link == 1 || link == 2){
		if(is.vector(multi)) rm = 1
		else rm = nrow(multi)
		De = matrix(0,J,rm)
		if(rm==1){
			De = 1
			fv = multi[1]
		}else{
			for(r in 1:rm){
				ind = multi[r,]
				ind = ind[ind>0]
				De[ind,r] = 1      
			}
			fv = multi[,1]     # list of constrained items
		}
		fve = (fv-1)*(l-1)+1
		indga = 1:J; indga = indga[-fv]
		indth = 1:(k*rm)
		if(difl==0){
			indbe = k*rm+(1:(J*(l-1)-rm))
			indbec = 1:(J*(l-1)); indbec = indbec[-fve]
		}else{
			indbe = k*rm+(1:(J-rm+l-2));
			indbec = 1:J; indbec = indbec[-fv]
		}
# find non redundant response configurations for each item
		conf1 = NULL; conf2 = NULL
		for(j in 1:J) for(h in 0:l){
			if(any(S[,j]==h)){
				conf1 = c(conf1,j)
				conf2 = c(conf2,h)
			}
		}
		nconf = length(conf1)
# abililities for each item
		if(rm==1) abils=rep(1,J) else{
			abils = rep(0,J)
			for(h in 1:rm){
				ind = multi[h,]; ind = ind[ind>0]
				abils[ind] = h
			}
		}
# design matrix
		if(difl==0) ZZ = array(0,c(l-1,k*rm+(l-1)*J-rm,J*k)) else{
			if(difl==1) ZZ = array(0,c(l-1,k*rm+J-rm+l-2,J*k))
		}
		cont = 0; refitem = matrix(0,J*k,1)       # reference item of that design matrix
		for(c in 1:k){
			u1 = matrix(0,1,k); u1[c] = 1
			for(j in 1:J){
				u2 = matrix(0,1,rm); u2[abils[j]] = 1
				v = matrix(0,1,J); v[j] = 1
				cont = cont+1
				if(difl==0){
					Te = v%x%diag(l-1)
					Te = matrix(Te[,-fve],l-1,dim(Te)[2]-length(fve))
				}else if(difl==1){
					Te = cbind(v%x%rep(1,l-1),diag(l-1))
					Te = Te[,-c(fv,J+1)]
				}				
				ZZ[,,cont] = cbind(rep(1,l-1)%*%(u1%x%u2),-Te)
				refitem[cont] = j
			}
		}
	    confe2 = rep(conf2,k)   # response configuration
		confe1 = conf1          # corresponding design matrix in ZZ
		for(c in 2:k) confe1 = c(confe1,max(confe1)+conf1)
   		ZZ0 = ZZ
	}
# When there is just 1 latent class
#   if k == 1,
#     piv = 1;
#     P = zeros(2,J);
#     for j in 1:J,
#       for jb = 0:1,
#         ind = which(S[,j]==jb);
#         P(jb+1,j) =  sum(yv(ind))/n;
#       end
#     end
#     Psi = P;
#     psi = ones(ns,1);
#     for j in 1:J,
#       psi = psi.*P(S[,j]+1,j);
#     end
#     lk = yv"*log(psi);
#     np = J;
#     aic = -2*lk+2*np;
#     bic = -2*lk+np*log(n);
#     bec = NULL; gac = NULL;
#     Pp = ones(ns,1);
#     Th = NULL;
#     return
#   end
	out = matr_glob(l); Co = out$Co; Ma = out$Ma
# Starting values
	if(start == 0){
		if(cov){
		    de = rep(0,(k-1)*(ncov+1)); piv = NULL
	    }else{
			be = NULL; piv = rep(1,k)/k
		} # latent class probabilities
		if(k==1) grid = 0 else grid = seq(-k,k,2*k/(k-1))
		Phi = array(0,c(l,J,k)) # probability every response
		for(j in 1:J){
			dist = rep(0,l)
			for(y in 0:(l-1)) dist[y+1] = sum(yv[S[,j]==y])/n
			eta = Co%*%log(Ma%*%dist)
			for(c in 1:k) Phi[,j,c] = inv_glob(eta+grid[c])$p
		}
	}
	if(start == 1){
		if(cov){
			de = rnorm((k-1)*(ncov+1))/5
			piv = NULL
		}else{
			piv = runif(k)
			piv = piv/sum(piv)
		}
		Phi = array(runif(l*J*k),c(l,J,k))
		for(c in 1:k) for(j in 1:J) Phi[,j,c] = Phi[,j,c]/sum(Phi[,j,c]);
	}  
	if(link==0) ga = NULL else if(link==1 || link==2) ga = rep(1,J-rm)
# Compute log-likelihood
	Psi = matrix(1,ns,k) # probability observed response
	if(miss){
		for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))
	}else{
		for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
	}
	if(cov){
	    Piv = matrix(0,ns,k)
	    out = prob_multilogit(XXdis,de,Xlabel,fort=fort)
	    Pdis = out$Pdis; Piv = out$P
	}else Piv = rep(1,ns)%o%piv
	if(k==1) Pj=Psi else Pj = Psi*Piv
	pm = rowSums(Pj)
	lk = sum(yv*log(pm))
	cat(c("Model with multdimensional structure\n"))
	print(multi)
	cat(c("Model of type =                ",link,"\n"))
	cat(c("Discrimination index =         ",disc,"\n"))
	cat(c("Constraints on the difficulty =",difl,"\n"))
	if(disc==0 || length(ga)==0){
    	cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    	cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n");
    	cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
	}else if(disc==1){
		cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
		cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n");
		cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");

	}
	cat(sprintf("%11g",c(0,k,link,lk)),"\n",sep=" | ")
 	it = 0; lko = lk-10^10; dis = 0; par = 0; dga = NULL
# Iterate until convergence
	while(((abs(lk-lko)/abs(lko)>tol) && it<10^4) || it<2){
		it = it+1
		paro = par; gao = ga; pivo = piv; deo = de; lko = lk
# ---- E-step ----
		V = ((yv/pm)%o%rep(1,k))*Piv*Psi; sV = colSums(V)
# ---- M-step ----
		if(link==0){  # LC model
			if(miss){
				for(j in 1:J) for(y in 1:l){
					ind = (S[,j]==(y-1))
					for(c in 1:k) Phi[y,j,c] = sum(V[ind,c]*R[ind,j])/sum(V[,c]*R[,j])
				}
			}else{
				for(j in 1:J) for(y in 1:l){
					ind = (S[,j]==y-1)
					for(c in 1:k) Phi[y,j,c] = sum(V[ind,c])/sV[c]
				}
			}
		}else{          # other models
			w = matrix(0,nconf,k)
			if(miss){
				for(conf in 1:nconf){
					j = conf1[conf]; h = conf2[conf]
					ind = (S[,j]==h)
					w[conf,] = colSums(V[ind,]*R[ind,j])
				}				
			}else{
				for(conf in 1:nconf){
					j = conf1[conf]; h = conf2[conf]
					w[conf,] = colSums(V[S[,j]==h,])
				}
			}
			if(disc==1){
				if(it>1 & rm<J){
					ZZ1 = array(0,c(l-1,J,J*k))
					count = 0
					for(c in 1:k) for(j in 1:J){
						count = count+1;
						ZZ1[,j,count] = ZZ0[,,count]%*%par
					}
					dimz = dim(ZZ1)
					dimz[2] = dimz[2]-length(fv)
					ZZ1 = array(ZZ1[,-fv,],dimz)
					ga = est_multi_glob(confe2,ZZ1,ltype,confe1,as.vector(w),ga)$be
				}
				gac = rep(1,J); gac[indga] = ga
				ZZ = ZZ0
				for(j in 1:J){
					ind = (refitem==j)
					ZZ[,,ind] = ZZ[,,ind]*gac[j]
				}
			}
#			print(length(confe2))
#			print(confe2)
#			print(dim(ZZ))
#			print(ZZ)
#			print(length(confe1))
#			print(confe1)
#			print(dim(w))
#			print(w)
#			return
			if(it==1) out = est_multi_glob(confe2,ZZ,ltype,confe1,as.vector(w))   # update par
			else out = est_multi_glob(confe2,ZZ,ltype,confe1,as.vector(w),par)
			par = out$be; lkc = out$lk; P = out$P
			Phi = array(t(P),c(l,J,k))
    		}
# Update piv
        if(cov){
            out = est_multilogit(V,XXdis,Xlabel,de,Pdis,fort=fort)
			de = out$be; Pdis = out$Pdis; Piv = out$P
        }else piv = sV/n
# Compute log-likelihood
		Psi = matrix(1,ns,k);
		if(miss){
			for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*(Phi[S[,j]+1,j,c]*R[,j]+(1-R[,j]))	
		}else{
			for(j in 1:J) for(c in 1:k)	Psi[,c] = Psi[,c]*Phi[S[,j]+1,j,c]
		}
		if(k==1) Pj=Psi else Pj = Psi*Piv
        pm = rowSums(Pj)
	    lk = sum(yv*log(pm))
		dis = max(c(abs(par-paro),abs(ga-gao),abs(piv-pivo)))
		if(it/10==floor(it/10)){
			if(disc==0 || length(ga)==0) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ") else{
				if(disc==1) cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
			}
		}
	}
	if(it/10>floor(it/10)){
		if(disc==0 || length(ga)==0){
			cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(par),max(par))),"\n",sep=" | ")
		    	cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
		}else if(disc==1){
			cat(sprintf("%11g",c(it,k,link,lk,lk-lko,dis,min(ga),max(ga),min(par),max(par))),"\n",sep=" | ")
			cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
		}
	}
# Compute number of parameters  
	if(link == 0){
	  np = k*J*(l-1)
      if(cov) np = np+(k-1)*(ncov+1) else np = np+k-1
    }else if(link==1 || link==2){
      np = k*rm+disc*(J-rm)
      if(cov) np = np+(k-1)*(ncov+1) else np = np+k-1
		if(difl==0) np = np+(l-1)*J-rm
		else if(difl==1) np = np+J-rm+l-2
	}
# extract parameter estimates  
  aic = -2*lk+2*np;
  bic = -2*lk+np*log(n);
  if(link==0){
    Th = NULL; Bec = NULL; gac = NULL; fv = NULL
  }
  else if(link==1 || link==2){
    th = par[indth]; be = par[indbe]
    if(difl==0){
      bec = rep(0,J*(l-1)); bec[indbec] = be
      Bec = t(matrix(bec,l-1,J))
    }
    else{
      bec1 = rep(0,J); bec1[indbec] = be[1:(J-rm)]
      bec2 = rep(0,l-1); bec2[2:(l-1)] = be[J-rm+(1:(l-2))]
      Bec = list(bec1=bec1,bec2 = bec2)
    }
    gac = rep(1,J); gac[indga] = ga
    Th = matrix(th,rm,k)
  }
  Pp = ((1./pm)%o%rep(1,k))*Piv*Psi
  if(cov) De = matrix(de,ncov+1,k-1) else De = NULL
  out = list(piv=piv,Th=Th,Bec=Bec,gac=gac,fv=fv,Phi=Phi,De=De,Piv=Piv,Pp=Pp,lk=lk,np=np,aic=aic,bic=bic)
}