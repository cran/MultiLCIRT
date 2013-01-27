      subroutine matr_YY(jj,k,ns,S,V,YY)

      integer j,jj,c,k,i,ns,S(ns,k),count,ind
      double precision V(ns,k),YY(jj*k,2)
            
      YY=0
      count = 0
	  do c=1,k
	    do j=1,jj
	      count = count+1
	      do i=1,ns
	        ind = S(i,j)+1
            YY(count,ind) = YY(count,ind)+V(i,c)
		  end do
		end do
      end do

      end 



