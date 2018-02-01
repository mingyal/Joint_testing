library("MASS")
library("survival")
library("mvtnorm")
library("lpSolve")
library("grpreg")
library("matrixcalc")
source("./grad_hessian.R")
source("./what_speedup.R")


var_check<-function(Hs, pos, nuisance_pos, what){
  var   = Hs[pos,pos] - t(what)%*%Hs[nuisance_pos,pos]
  eig_var = eigen(var)
  if(max(eig_var$values)<1e-8){var=(eig_var$vectors)%*%(diag(K)*1e-8)%*%t((eig_var$vectors))}
  #if(   min(eig_var$values)<0 ) {stop("var is not pd")}
  return(var)
}

groupstat <- function(x,time,status, groups, task, alternative, coeff, lambdasto, lambdasby, lambda_proj){

  d    = dim(x)[2]; n = dim(x)[1]; K = d/length(unique(groups)); numgp = d/K     
  if(   (numgp%%1!=0)  )  {stop("Number of groups are not integers!")}
  if(task=="fdr"){  indx = seq(from = 1,to=numgp,by = 1)    }
  if(task=="type-1"){  indx = c(alternative)}
  lambdas = seq(from = 0.1,to=lambdasto,by = lambdasby); lambdas = lambdas*log(d)/n          
  cvfit <- cv.grpsurv(x, Surv(time, status), groups, returnY=FALSE,trace=FALSE, penalty="grLasso", lambda = lambdas, maxit = 10000,nfolds = 20,alpha=1)
  #lambdamin = cvfit$lambda.min  
  statistics=c()
  beta    = coef(cvfit)   # Initial Estimator
  #--------------------------------------------------------------------------
  # cat('non-zero gps are: ')
  # nonzerogp=which(beta[seq(from=1, to=d-K+1, by=K)]!=0)
  # print(as.vector(matrix(nonzerogp, nrow = 1)))
  #--------------------------------------------------------------------------
  otime = order(time) 
  tmp=grad_hessian(x,time,status, n, d, otime, beta)
  derivative = tmp[[1]];  Hs = tmp[[2]]
  for (coi in indx){
     
      cat(coi)
      pos = seq(from = (coi-1)*K+1,to=(coi-1)*K+K,by = 1)  # positions corresponding to group coi
      coi_not_zero = (sum( abs( beta[pos] ) )>1e-10)
      nuisance_pos = setdiff(1:d, pos)
      beta_undernull=beta
      beta_undernull[pos] = 0
      #--------------------------------------------------------------------------
      la=derivative[,pos];  lb=derivative[, setdiff(1:d, pos)   ];
      what=getwhat(n, d,  K, beta, pos, Hs, lambda_proj)
      var   = var_check(Hs, pos, nuisance_pos, what)
      if(  coi_not_zero  ){
        tmp=grad_hessian(x,time,status, n, d, otime, beta_undernull);
        Hs_undernull=tmp[[2]];
        what_undernull=getwhat(n, d,  K, beta, pos, Hs_undernull, lambda_proj)
        var_undernull=var_check(Hs_undernull, pos, nuisance_pos, what_undernull)
      }else{   var_undernull=var    }
      #--------------------------------------------------------------------------
      S = beta[pos] - solve(  var  , (  la  - t(what)%*%lb   )  )  
      #print(solve(  var  , (  la  - t(what)%*%lb   )  )  )   # print the debiased term
      #print(n*(t(S)) %*% var_undernull %*% S)
      statistics   = c(statistics, n*(t(S)) %*% var_undernull %*% S)
      #pval[[coi]]  = 1-pchisq(statistics[[coi]],df=K)
  }
  return(statistics)
}