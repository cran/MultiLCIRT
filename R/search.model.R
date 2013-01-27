search.model <- function(S, yv = rep(1,ns), kv, X = NULL, link = 0, disc = 0,
	difl = 0, multi = 1:J, fort = FALSE, tol = 10^-10)
{

# function that search for the global maximum of the log-likelihood
# vector of kv to try for
ns = dim(S)[1]
J = dim(S)[2]
out = vector("list",max(kv))
bicv = rep(0,max(kv)); lkv = bicv
for(k in kv){
  print("***************************************************************************")
  print(k)
  out[[k]] = est_multi_poly(S=S,yv=yv,k=k,X=X,start=0,link=link,disc=disc,difl=difl,multi=multi,fort=fort,tol=tol)
  lktrace = out[[k]]$lk
  lkv[k] = out[[k]]$lk
  bicv[k] = out[[k]]$bic
  print(sort(lktrace))
  print(lkv)
  print(bicv)
  if(k>1) for(h in 1:(2*(k-1))){
    print("***************************************************************************")
  	print(c(k,h))
    outh = est_multi_poly(S=S,yv=yv,k=k,X=X,start=1,link=link,disc=disc,difl=difl,multi=multi,fort=fort,tol=tol)
    lktrace = c(lktrace,outh$lk)
    if(outh$lk>out[[k]]$lk) out[[k]] = outh  	
    lkv[k] = out[[k]]$lk
    bicv[k] = out[[k]]$bic
    print(sort(lktrace))
    print(lkv)
    print(bicv)
  }
  out[[k]]$lktrace = lktrace
}
out = list(out.single=out,bicv=bicv,lkv=lkv)

}