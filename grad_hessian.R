

grad_hessian<-function(x,time,status, n, d, otime, beta){    
    S0  = 0; S1  = rep(0,d) ; S2  = matrix(rep(0,d*d),nrow = d)
    Hs  = matrix(rep(0,d*d),nrow = d)                                # Hessian
    derivative=matrix(rep(0,d*1),ncol = d) 
    previous_position = n
    invseq = seq(from=n, to = 1, by=-1)
    for(i in invseq){
      if (status[otime[i]]==1){
        ind = otime[ i:previous_position ] 
        for (j in ind){
          tmp = as.numeric(exp(x[j,]%*%beta))/n
          S0  = S0 + tmp
          S1  = S1 + tmp%*%t(x[j,])
          S2  = S2 + tmp*x[j,]%*%t(x[j,])
        }
        derivative= derivative + (- (x[otime[i],] - S1/S0) / n )
        #if( S0<=0 ){stop("S0 is not positive")}
        #if( ! is.positive.semi.definite( t(S1)%*%S1 ) ){stop("S1S1' is not psd")}
        #if( ! is.positive.semi.definite( S2 ) ){stop("S2 is not psd")}
        #if( ! is.positive.semi.definite( S2/S0- (t(S1)%*%S1)/(S0^2) ) ){stop("term is not psd")}
        Hs= Hs + ( ( S2/S0- (t(S1)%*%S1)/(S0^2) ) / n )
        previous_position=i-1
      }
    }
    #la = derivative[,pos]
    #lb = derivative[, setdiff(1:d, pos)   ]
    r = list()
    r[[1]]=derivative; r[[2]]=Hs;
    return(r)
}  
    
      
      
        
        
      