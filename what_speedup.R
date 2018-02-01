#Note patten_in_lp(K-1) is corresponding to K
patten_in_lp<- function(k){
  tmp = diag(k)
  z = matrix( rep(0, k*k), nrow = k    ) 
  mat = NULL
  for(i in 1:k){ 
    mat=cbind( mat,  tmp[,i], z[,i]       )  
  }
  return(mat) 
}

#---------------------------------------------------------------------------

getwhat<-function(n, d, K, beta, pos, Hs, lambda_proj){    
  nuisance_pos=setdiff(1:d, pos)
  Hab = Hs[  pos , nuisance_pos ,drop=FALSE]
  Hbb = Hs[ nuisance_pos  , nuisance_pos ]
  what = matrix(0, (d-K), K);
  tmp = matrix(rep(0,(d-K)*(d-K)),nrow = d-K)
  A1  = cbind(Hbb,tmp)
  A2  = cbind(-Hbb,tmp)
  A3  = cbind(diag(d-K),-diag(d-K))
  A4  = cbind(-diag(d-K),-diag(d-K))
  A_onegp   = rbind(A1,A2,A3,A4)
  A = kronecker(diag(K), A_onegp)
  if(K>1){
    sym = NULL
    for(q in 1:(K-1)){
      if(q>1){  previous_zeros = matrix(rep(0, 2*(K-q)*2*(q-1)*(d-K)), nrow = 2*(K-q))  }else{previous_zeros=NULL}
      curw=rbind(  Hab[(q+1):K, ,drop=FALSE] , -Hab[(q+1):K, ,drop=FALSE]     )
      curw=cbind(curw, matrix(  0, dim(curw)[1], dim(curw)[2]  )    )
      otherw = rbind(  patten_in_lp(K-q) , -patten_in_lp(K-q)     )
      otherw = kronecker(otherw, -Hab[q, ,drop=FALSE])
      sym = rbind(sym,  cbind(previous_zeros,curw, otherw ))
    }
    A = rbind(A, sym)
  }
  obj = NULL;   for(q in 1:K){   obj=c(obj, c(rep(0,d-K),rep(1,d-K)))   }
  dir = rep("<=",    ( (d-K)*4 )*K+ K*(K-1)   )   
    
  rhs = NULL;   for(q in 1:K){    rhs=c(rhs, Hab[q,] + lambda_proj , -Hab[q,] + lambda_proj , rep(0,d-K),rep(0,d-K))              }
  rhs = c(rhs, rep(   0,   K*(K-1))     )
  lpsol = lp(direction = "min",objective.in = obj,const.mat = A,const.dir = dir, const.rhs = rhs)
  #if(lpsol$status==2){ stop("Infeasible for the lp")}
  if(lpsol$status==2){
    lambda_proj_surrogate=lambda_proj
    while(lpsol$status==2){
      lambda_proj_surrogate=lambda_proj_surrogate*2
      rhs = NULL;   for(q in 1:K){    rhs=c(rhs, Hab[q,] + lambda_proj_surrogate , -Hab[q,] + lambda_proj_surrogate , rep(0,d-K),rep(0,d-K))              }
      rhs = c(rhs, rep(   0,   K*(K-1))     )
      lpsol = lp(direction = "min",objective.in = obj,const.mat = A,const.dir = dir, const.rhs = rhs)
    }
  }
  #if(sum(abs( lpsol$solution) )<1e-15){stop("w is 0!!!")}
  for(q in 1:K){  what[,q] = lpsol$solution[  ( 1+2*(q-1)*( d-K ) ):  ( ( d-K ) +2*(q-1)*( d-K ))  ]  }
  return(what)
}





 