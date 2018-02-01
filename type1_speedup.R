if(FALSE){
  rm(list=ls())

  library("MASS")
  library("mvtnorm")
  library("lpSolve")
  library("grpreg")
  library("matrixcalc")
  library(magic)
  source("/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/rjointtesting/speedup/grad_hessian.R")
  source("/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/rjointtesting/speedup/what_speedup.R")
  source("/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/rjointtesting/speedup/groupstat_speedup.R")

  
  patten_in_lp<- function(k){
    tmp = diag(k)
    z = matrix( rep(0, k*k), nrow = k    ) 
    mat = NULL
    for(i in 1:k){ 
      mat=cbind( mat,  tmp[,i], z[,i]       )  
    }
    return(mat) 
  }

  n = 50
  ngp = 50
  s=2
  K = 4
  rou = 0.25
  task="type-1"
  coeff=1.5*sqrt(log(ngp)/n)
  alternative=1
  Sigma = (rou)^(toeplitz (1:(ngp*K))-1)
  x=MASS::mvrnorm(n, matrix(0, (ngp*K), 1), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  beta_gold = matrix(0,K*ngp, 1);
  activefeatures = 2:(s+1)      
  pos = matrix(0,K*ngp, 1)
  pos[  (K*(ngp-s)+1):(K*ngp), ]=1
  beta_gold[pos==1,]=coeff
  survival_time =rexp(  n ,rate = exp(x%*%beta_gold)  )
  censor_time = rexp(  n ,rate = (  exp(x%*%beta_gold)*runif(n, min=0, max=1) ) )   
  time = pmin(censor_time, survival_time) 
  status = as.integer((censor_time>survival_time) == TRUE)#n
  groups <- ceiling(   (1:(K*ngp))   / K) 
  otime = order(time)  
  d    = dim(x)[2]        
  lambdas = seq(from = 0.1,to=500,by = 3) 
  lambdas = lambdas*log(d)/n                 
  cvfit <- cv.grpsurv(x, Surv(time, status), groups, returnY=FALSE,trace=FALSE, penalty="grLasso", lambda = lambdas, maxit = 10000,nfolds = 20,alpha=1)
  beta    = coef(cvfit)   # Initial Estimator
  beta_undernull=beta
  coi = 1
  pos = seq(from = (coi-1)*K+1,to=(coi-1)*K+K,by = 1)
  beta_undernull[pos] = 0
  
  
  #--------------------------------------------------------------------------
  tmp=grad_hessian(x,time,status, n, d, otime, beta, pos)
  la=tmp[[1]];lb=tmp[[2]]; Hs=tmp[[3]];
  if(length(tmp)>3){
    tmp=grad_hessian(x,time,status, n, d, otime, beta_undernull, pos)
    Hs_undernull=tmp[[3]]
  }else{Hs_undernull=Hs}
  tmp=getwhat(n, d, s, K, beta, pos, Hs, Hs_undernull)
}





# Note when we build betatilde we use score and hessin with alphahat
# But, when normalize the statistic, we should use hessin with 0 instead of alphahat!!!
groupstat <- function(x,time,status, groups, task, alternative, coeff)
{
  
  Pval = c()
  d    = dim(x)[2]                    # Note d equals the product of K and number of groups
  n    = dim(x)[1]
  K = d/length(unique(groups))          # Groupsize
  numgp = d/K                           # Number of groups
  if(   (numgp%%1!=0)  )  {stop("Number of groups are not integers!")}
  if(task=="fdr"){     indx = seq(from = 1,to=numgp,by = 1)      }
  if(task=="type-1"){  indx = c(alternative)                     }
  lambdas = seq(from = 0.1,to=1.5,by = 0.03) 
  lambdas = lambdas*log(d)/n                 ##Set of tuning parameters
  
  cvfit <- cv.grpsurv(x, Surv(time, status), groups, returnY=FALSE,trace=FALSE, penalty="grLasso", lambda = lambdas, maxit = 10000,nfolds = 20,alpha=1)
  lambdamin = cvfit$lambda.min/(log(d)/n  )  # cross-validate lambda for grplasso

  Hs=list()
  Hs_undernull=list()
  la=list()
  lb=list()
  S0=list()
  S0_undernull=list()
  S1=list()
  S1_undernull=list()
  S2=list()
  S2_undernull=list()
  statistics=list()
  pval = list()
  
  for (coi in indx)                                     # coi stands for the group we are looking for.
  {
    beta    = coef(cvfit)   # Initial Estimator
    betas   = beta
    pos = seq(from = (coi-1)*K+1,to=(coi-1)*K+K,by = 1)  # positions corresponding to group coi
    beta_undernull=beta
    beta_undernull[pos] = 0
    stime = sort(time)                                   # Sorted survival/censored times
    otime = order(time)                                  # Order of time
    Vs  = matrix(rep(0,d*d),nrow = d)
    Hs[[coi]]  = Vs                                 # Hessian
    Hs_undernull[[coi]]=Vs
    
    ind = 0
    
    la[[coi]]  = matrix(rep(0,(K)*n),nrow = n)      # Gradient w.r.t parameter of interest
    lb[[coi]]  = matrix(rep(0,(d-K)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
    i   = 1
    while( i<=n)
    {
      if (status[otime[i]]==1)
      {
        ind = which(time >= stime[i])
        S0[[coi]]  = 0
        S0_undernull[[coi]] = 0
        S1[[coi]]  = rep(0,d)
        S1_undernull[[coi]]=rep(0,d)
        S2[[coi]]  = matrix(rep(0,d*d),nrow = d)
        S2_undernull[[coi]] =  matrix(rep(0,d*d),nrow = d)
        
        if (length(ind)>0)
        {
          for (j in 1:length(ind))
          {
            tmp = exp(x[ind[j],]%*%betas)
            tmp_undernull=exp(x[ind[j],]%*%beta_undernull)
            S0[[coi]]  = S0[[coi]] + tmp
            S0_undernull[[coi]]  = S0_undernull[[coi]] + tmp_undernull
            
            S1[[coi]]  = S1[[coi]] + tmp%*%t(x[ind[j],])
            S1_undernull[[coi]]  = S1_undernull[[coi]] + tmp_undernull%*%t(x[ind[j],])
            
            tmp = apply(tmp,1,as.numeric)   
            tmp_undernull = apply(tmp_undernull,1,as.numeric) 
            
            S2[[coi]]  = S2[[coi]] + tmp*x[ind[j],]%*%t(x[ind[j],])
            S2_undernull[[coi]]  = S2_undernull[[coi]] + tmp_undernull*x[ind[j],]%*%t(x[ind[j],])
          }
        }
        S0[[coi]] = apply(S0[[coi]],1,as.numeric)    #We finish computing S0,S1,S2 here for group coi.
        S0_undernull[[coi]] = apply(S0_undernull[[coi]],1,as.numeric)
        
        la[[coi]][i, ]  = -(x[otime[i],pos] - S1[[coi]][pos]/S0[[coi]])
        if (coi == 1)
        {
          lb[[coi]][i,] = -(x[otime[i],c(  (coi*K+1):d)] - S1[[coi]][c((coi*K+1):d)]/S0[[coi]])
        } else if (coi == numgp){
          lb[[coi]][i,] = -(x[otime[i],c(1:(  (coi-1)*K  )  )] - S1[[coi]][c(1:(  (coi-1)*K  ))]/S0[[coi]])
        } else {
          lb[[coi]][i,] = -(x[otime[i],c(1:(  (coi-1)*K  ) , (  (coi)*K+1  ) :d)] - S1[[coi]][c(1:(  (coi-1)*K  ) , (  (coi)*K+1  ) :d)]/S0[[coi]])
        }                                             #We finish computing la,lb here for group coi.
        V   = S0[[coi]]*S2[[coi]] - t(S1[[coi]])%*%(S1[[coi]])
        
        V_undernull   = S0_undernull[[coi]]*S2_undernull[[coi]] - t(S1_undernull[[coi]])%*%(S1_undernull[[coi]])
        
        Hs[[coi]]  = Hs[[coi]] + V/(n*S0[[coi]]^2)       
        Hs_undernull[[coi]]  = Hs_undernull[[coi]] + V_undernull/(n*S0_undernull[[coi]]^2)    
        
      }
      i = i + 1
    }
    
    #if( ! (is.positive.semi.definite( Hs[[coi]], tol=1e-8))  ){stop(paste0("The ",coi," -th Hs is not positive definite"))}
    #if(max(eigen(Hs[[coi]])$values)<=1e-8){Hs[[coi]]=diag(d)*(1e-8)}
    #if(  ! (is.positive.semi.definite( Hs_undernull[[coi]], tol=1e-8))  ) {stop(paste0("The ",coi," -th Hs(undernull) is not positive definite"))}
    #if(max(eigen(Hs_undernull[[coi]])$values)<=1e-8){Hs_undernull[[coi]]=diag(d)*(1e-8)}
    
    #Hs = Hs/n
    
    # Compute Hab, Hbb for computing what by Dantzig
    if (coi == 1)
    {
      Hab = Hs[[coi]][  (coi*K+1):d , ((coi-1)*K+1):(coi*K),drop=FALSE]
      Hbb = Hs[[coi]][  (coi*K+1):d  ,  (coi*K+1):d  ]
      Hab_undernull = Hs_undernull[[coi]][  (coi*K+1):d , ((coi-1)*K+1):(coi*K),drop=FALSE]
      Hbb_undernull = Hs_undernull[[coi]][  (coi*K+1):d  ,  (coi*K+1):d  ]
    } else if (coi == numgp){
      Hab = Hs[[coi]][  1:((coi-1)*K)  ,  ((coi-1)*K+1):(coi*K) ,drop=FALSE ]
      Hbb = Hs[[coi]][1:((coi-1)*K),1:((coi-1)*K)]
      Hab_undernull = Hs_undernull[[coi]][  1:((coi-1)*K)  ,  ((coi-1)*K+1):(coi*K) ,drop=FALSE ]
      Hbb_undernull = Hs_undernull[[coi]][  (coi*K+1):d  ,  (coi*K+1):d  ]
    } else{
      Hab = Hs[[coi]][c(  1:((coi-1)*K),(coi*K+1):d  ),  ((coi-1)*K+1):(coi*K),drop=FALSE]
      Hbb = Hs[[coi]][c(  1:((coi-1)*K),(coi*K+1):d  ),c(  1:((coi-1)*K),(coi*K+1):d  )]
      Hab_undernull = Hs_undernull[[coi]][c(  1:((coi-1)*K),(coi*K+1):d  ),  ((coi-1)*K+1):(coi*K),drop=FALSE]
      Hbb_undernull = Hs_undernull[[coi]][c(  1:((coi-1)*K),(coi*K+1):d  ),c(  1:((coi-1)*K),(coi*K+1):d  )]
    }
    what = matrix(0, (d-K), K);
    what_undernull = matrix(0, (d-K), K);
    for(i in 1:K){
      tmp = matrix(rep(0,(d-K)*(d-K)),nrow = d-K)
      A1  = cbind(Hbb,tmp)
      A2  = cbind(-Hbb,tmp)
      A3  = cbind(diag(d-K),-diag(d-K))
      A4  = cbind(-diag(d-K),-diag(d-K))
      
      A   = rbind(A1,A2,A3,A4)
      obj = c(rep(0,d-K),rep(1,d-K))
      dir = rep("<=",(d-K)*2 + (d-K)*2)
      lambda = 0#s*sqrt(log(d)/n)
      rhs = c(Hab[i,] + lambda, -Hab[i,] + lambda, rep(0,d-K),rep(0,d-K))
      tmp = lp(direction = "min",objective.in = obj,const.mat = A,const.dir = dir, const.rhs = rhs)
      what[,i] = tmp$solution[1:(d-K)]
      #--------------------------------------------------------------------------------------------------------
      tmp_undernull = matrix(rep(0,(d-K)*(d-K)),nrow = d-K)
      A1_undernull  = cbind(Hbb_undernull,tmp_undernull)
      A2_undernull  = cbind(-Hbb_undernull,tmp_undernull)
      A3_undernull  = cbind(diag(d-K),-diag(d-K))
      A4_undernull  = cbind(-diag(d-K),-diag(d-K))
      
      A_undernull   = rbind(A1_undernull,A2_undernull,A3_undernull,A4_undernull)
      obj_undernull = c(rep(0,d-K),rep(1,d-K))
      dir_undernull = rep("<=",(d-K)*2 + (d-K)*2)
      lambda = s*sqrt(log(d)/n)
      rhs_undernull = c(Hab_undernull[,i] + lambda, -Hab_undernull[,i] + lambda, rep(0,d-K),rep(0,d-K))
      tmp_undernull = lp(direction = "min",objective.in = obj_undernull,const.mat = A_undernull,const.dir = dir_undernull, const.rhs = rhs_undernull)
      what_undernull[,i] = tmp_undernull$solution[1:(d-K)]
    }
    
    # Decorrelated Wald
    # pos = seq(from = (coi-1)*K+1,to=coi*K,by = 1) 
    wholeset = 1:d;
    nuisance_pos = wholeset[wholeset!=pos]
    var   = Hs[[coi]][pos,pos] - t(what)%*%Hs[[coi]][nuisance_pos,pos]
    var_undernull   = Hs_undernull[[coi]][pos,pos] - t(what_undernull)%*%Hs_undernull[[coi]][nuisance_pos,pos]
    
    svd_var = svd(var)
    if(max(svd_var$d)<1e-8){var=(svd_var$u)%*%(diag(d)*1e-8)%*%t(svd_var$v)}
    if(   min(svd_var$d)<0 )  {stop("var is not psd")}
    
    svd_var_undernull=svd(var_undernull)
    if(max(svd_var_undernull$d)<1e-8){var=(svd_var_undernull$u)%*%(diag(d)*1e-8)%*%t(svd_var_undernull$v)}
    if(   min(svd_var$d)<0 )  {stop("var_null is not psd")}
    
    
    S = beta[pos] - solve(  var  , (  colMeans(la[[coi]])  - t(what)%*%(colMeans(lb[[coi]]))  )  )
    #statistics[[coi]]   = (n)*S^2*(max(var,1e-8))
    #print lasso est, debias term and var
    #------------------------------------------------------------------------------------------------------------------------
    # cat('lasso est: ')
    # cat(beta[pos])
    # cat(',   debias term: ')
    # cat(- solve(  var  , (  colMeans(la[[coi]])  - t(what)%*%(colMeans(lb[[coi]]))  )  ))
    # cat(',  var_undernull')
    # print(as.numeric(var_undernull))
    #------------------------------------------------------------------------------------------------------------------------
    
    if(indx!=1){statistics[[coi]]   = (n)*(t(S-coeff))%*%(var_undernull)%*%((S-coeff))} else{statistics[[coi]]   = (n)*(t(S))%*%(var_undernull)%*%(S)}
    #print(statistics[[coi]])
    #print(paste0("K is ", K))
    #print(paste0("statistics is ",statistics[[coi]]))
    #pval[[coi]]  = 1-pchisq(statistics[[coi]],df=K)
    #Pval = c(Pval,pval)
  }
  if(task=="type-1"){return(statistics[[coi]])}
  if(task=="fdr"){return(statistics)}
}





comparision_simulation <- function(n, alpha, s, ngp, method='D', K, task,  alternative, coeff){
  rou = 0.25
  #toeplitz matrix
  Sigma = (rou)^(toeplitz (1:(ngp*K))-1)
  #Sigma = diag(ngp*K)
  sample=MASS::mvrnorm(n, matrix(0, (ngp*K), 1), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  #set true beta
  beta = matrix(0,K*ngp, 1);
  activefeatures = 2:(s+1)                  #fix the beta!!!
  pos = matrix(0,K*ngp, 1)
  pos[  2:(s*K+1)  ]=1
  
  
  if(method =='D'){
    beta[pos==1,]=coeff
  }else{
    beta[pos==1,]=coeff*runif(s*K, min=1/3, max=4/3) 
  }
  #generate delta
  survival_time =rexp(  n ,rate = exp(sample%*%beta)  )
  censor_time = rexp(  n ,rate = (  exp(sample%*%beta)*runif(n, min=0, max=1) ) )   #(runif(n, min=0, max=1) )
  time = pmin(censor_time, survival_time) 
  status = as.integer((censor_time>survival_time) == TRUE)#n
  groups <- ceiling(   (1:(K*ngp))   / K) 
  #------------------------------------------------------------------------------------------------------------------------------------
  rgp = groupstat(sample,time,status, groups, task,alternative, coeff)
  return(rgp)
  
  if(task=="fdr"){
    return(fdr(rgp, ngp, K, alpha, activefeatures))
  }
}