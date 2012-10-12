prob_multilogit <- function(Xdis,be,label,fort=F){

     k = dim(Xdis)[1]
     ndis = max(label)
     n = length(label)
     ncov = length(be)
     Pdis = matrix(0,ndis,k); P = matrix(0,n,k)
     if(fort==F){
         for(i in 1:ndis){
            pdis = exp(Xdis[,,i]%*%be); pdis = pdis/sum(pdis)
            Pdis[i,] = pdis
     		mul = sum(label==i)
     		P[label==i,] = rep(1,mul)%o%pdis
          }
     }else{
	     out = .Fortran("prob_multilogit",Xdis=Xdis,be=be,label=as.integer(label),Pdis=Pdis,P=P,k=as.integer(k),ndis=as.integer(ndis),
	     ns=as.integer(n),ncov=as.integer(ncov))
	     Pdis = out$Pdis; P = out$P
     }
     out = list(Pdis=Pdis,P=P)
     return(out)
}