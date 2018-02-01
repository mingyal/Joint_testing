rm(list=ls())
Sys.setenv(TZ='EST')
library("MASS")
library("survival")
library("mvtnorm")
library("lpSolve")
library("grpreg")
library("matrixcalc")
library('foreach')
library('doParallel')
mac_path <- "/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/code/speedup/"
win_path <- "/Users/mingyal/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/code/speedup/"
setwd(mac_path)
source("grad_hessian.R")
source("what_speedup.R")
source("groupstat_speedup.R")
source("fdr_speedup.R")




Main<-function( n ,ngp,s,K ,rou ,task,coeff,alternative, lambdasto, lambdasby, lambda_proj){
  
  #n = 5;ngp = 3;s=2;K = 2;rou = 0.25;task="type-1";coeff=1;alternative=1
#--------------------------------------------------------------------------
  #Model 1
  #Sigma = (rou)^(toeplitz (1:(ngp*K))-1)   
  #Sigma = diag(ngp*K)
#--------------------------------------------------------------------------
  #Model 3
  # v = runif(ngp*K,1,3); D = diag(v);
  # Omega = diag(ngp*K);
  # tmp1 = diag(ngp*K-1)*0.6; tmp1 = cbind(matrix(0, nrow = ngp*K, ncol = 1), rbind(tmp1,matrix(0, nrow = 1, ncol = ngp*K-1) )  ); sub1 = tmp1+t(tmp1);
  # tmp2 = diag(ngp*K-2)*0.3; tmp2 = cbind(matrix(0, nrow = ngp*K, ncol = 2), rbind(tmp2,matrix(0, nrow = 2, ncol = ngp*K-2) )  ); sub2 = tmp2+t(tmp2);
  # Omega = Omega + sub1+sub2;
  # Omega = sqrt(D)%*%Omega%*%sqrt(D);
  # Sigma = solve(Omega);
#--------------------------------------------------------------------------
  #Model 2
  
  groupmat = matrix(withingp, nrow = K, ncol = K); diag(groupmat)<-1;
  Sigma = kronecker(diag(ngp), groupmat);
  Sigma[Sigma == 0] <- abs(rnorm(ngp*(ngp-1)*K^2, 0, 1*between))
  Sigma = ( Sigma+t(Sigma) )/2
  
#--------------------------------------------------------------------------
  
  x=MASS::mvrnorm(n, matrix(0, (ngp*K), 1), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  beta_gold = matrix(0,K*ngp, 1);
  pos = matrix(0,K*ngp, 1)
  pos[  (K*(ngp-s)+1):(K*ngp), ]=1
  beta_gold[pos==1,]=coeff
  survival_time =rexp(  n ,rate = exp(x%*%beta_gold)  )
  censor_time = rexp(  n ,rate = (  exp(x%*%beta_gold)* runif(n, min=0, max=0.5) ) )      #
  time = pmin(censor_time, survival_time) 
  status = as.integer((censor_time>survival_time) == TRUE)#status=1(sur<=cen)
  groups <- ceiling(   (1:(K*ngp))   / K) 
  return( groupstat(x,time,status, groups, task, alternative, coeff, lambdasto, lambdasby, lambda_proj) )
}
#-------------------------------------------------------------------------
withingp=0.25; 
between = 0
rou = 0.75;  
#--------------------------------------------------------------------------
K = 1
n = 100;    ngp = 100;   s= 5; d = K*ngp  ; ncore=3
#--------------------------------------------------------------------------
task="type-1"; 
#task="fdr";
#--------------------------------------------------------------------------
coeff=0.5;       alternative=1;     Notconservative=1; 
#--------------------------------------------------------------------------
lambdasto=500;     lambdasby=3;    lambda_proj_mulple = 5;
#---------2----------------------------------------------------------------
num_of_batch= 5
num_per_batch = 200
#--------------------------------------------------------------------------
#coeff=1.5*sqrt(log(ngp)/n)
lambda_proj=lambda_proj_mulple*sqrt(log(K*ngp)/n)
if(task=="type-1"){
  type1=c()
  for(t in 1:num_of_batch){
    start_time<-Sys.time()
    cat('(task:'); cat(task);cat(',#cores='); cat(ncore);cat(',n='); cat(n); cat(',ngp='); cat(ngp); cat(',s='); cat(s); cat(',K='); cat(K); cat(',num_per_batch='); cat(num_per_batch);cat(',lmdato='); cat(lambdasto);cat(',lmdaby='); cat(lambdasby);cat(',coef='); cat(coeff);cat(',lmdamtple='); cat(lambda_proj_mulple);cat(', start_time: ');cat(as.character(start_time)); cat(')');  
#    if(FALSE){
      cl<-makeCluster(ncore)
      registerDoParallel(cl)
      batch<- foreach( times=1:num_per_batch, .combine = c,.packages=c("MASS", "survival", "mvtnorm", "lpSolve", "grpreg", "matrixcalc")  ) %dopar%
                      Main(n ,ngp,s, K ,rou ,task,coeff,alternative, lambdasto, lambdasby,lambda_proj)
      type1=c(type1, batch)
      stopCluster(cl)
      end_time<-Sys.time()
      print(end_time-start_time)
#    }
#    type1=c(type1, Main(n ,ngp,s, K ,rou ,task,coeff,alternative, lambdasto, lambdasby,lambda_proj)  )
    print( sum(type1>qchisq(0.95, K)) )
    
  }
  hist(type1,  breaks=50, prob=TRUE, xlab = "x", ylab="Prob",main = "Empircial histogram", col='blue')
  curve( dchisq(x, df=K), col='red', add=TRUE)
  #write.csv(type1, '/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/report/data/new/fdr/100_x_x_x.csv', row.names=FALSE)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------
fdrsummary<-function(ALLstat, ngp, s, K, alpha, Notconservative){
  chi2absnorm=qnorm(1-  (  (  pchisq(as.matrix(  ALLstat  ), K, ncp = 0, lower.tail = FALSE, log.p = FALSE)  )  /2  ),0,1)
  wholeset = 1:(ngp); H0=1:(ngp-s); H1 = setdiff(wholeset, H0); 
  falsediscover_power<-apply(chi2absnorm, 1, fdr, ngp,  alpha, H0, H1,Notconservative )
  return(rowMeans(falsediscover_power))
}


if(task=="fdr"){
  ALLstat = NULL
  for(t in 1:num_of_batch){
    start_time<-Sys.time()
    #Model2
    cat('(task:'); cat(task);cat(',#cores='); cat(ncore);cat(',n='); cat(n); cat(',ngp='); cat(ngp); cat(',s='); cat(s); cat(',K='); cat(K); cat(',num_per_batch='); cat(num_per_batch);cat(',lmdato='); cat(lambdasto);cat(',lmdaby='); cat(lambdasby);cat(',coef='); cat(coeff);cat(',lmdamtple='); cat(lambda_proj_mulple);cat(', withingp='); cat(withingp); cat(', between='); cat(between);cat(', start_time: ');cat(as.character(start_time)); cat(')');  
    #Model1
    #cat('(task:'); cat(task);cat(',#cores='); cat(ncore);cat(',n='); cat(n); cat(',ngp='); cat(ngp); cat(',s='); cat(s); cat(',K='); cat(K); cat(',num_per_batch='); cat(num_per_batch);cat(',lmdato='); cat(lambdasto);cat(',lmdaby='); cat(lambdasby);cat(',coef='); cat(coeff);cat(',lmdamtple='); cat(lambda_proj_mulple);cat(', rou='); cat(rou);cat(')');  cat('  ');cat(t); cat(': ');
    
    cl<-makeCluster(ncore)
    registerDoParallel(cl)
    batch<- foreach( times=1:num_per_batch, .combine = rbind,.packages=c("MASS", "survival", "mvtnorm", "lpSolve", "grpreg", "matrixcalc")  ) %dopar%
      Main(n ,ngp,s, K ,rou ,task,coeff,alternative, lambdasto, lambdasby,lambda_proj)
    ALLstat=rbind(ALLstat, batch)
    stopCluster(cl)
    #ALLstat=rbind(ALLstat, Main(n ,ngp,s, K ,rou ,task,coeff,alternative)  )
    end_time<-Sys.time()
    print(end_time-start_time)
    print( fdrsummary(ALLstat, ngp, s, K, 0.05, Notconservative) )
  }
  #write.csv(fdr, '/Users/lmy/Dropbox/Mingyang-HLi/projects/joint_testing_and_FDR_control/code/speedup/result/fdr/100_100_5_3.csv', row.names=FALSE)
  row.names(ALLstat)<-1:(num_of_batch*num_per_batch)
  #print(ALLstat)
}

#model1(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,rou = 0.25)
#100/100/5/1  0.027/0.901(0.53, 1.41) 0.03083333 0.64800000(401:500)(checked)
#100/100/5/2  0.034/0.889(0.98, 1.73) 0.03066667 0.65400000(201:300)(checked)              |   0.05 0.32  (1:100)(checked)
#100/100/5/3  0.047/0.818(1.51, 1.93) 0.05483333 0.50800000(1:100)  (checked)              |   0.068 0.127(1:100)(checked)
#100/100/5/4  0.059/0.777(1.94, 2.13) 0.07366667 0.49600000(101:200)(checked)              |   0.058 0.069(1:100)(checked)
#100/200/5/1  0.028/0.883(0.53, 1.52) 0.0335 0.5740        (1:100)  (checked)         
#100/200/5/2  0.043/0.886(1.02, 1.84) 0.05433333 0.54400000(101:200)(checked)              |   0.048 0.222(1:100)(checked)
#100/200/5/3  0.046/0.820(1.58, 1.94) 0.05816667 0.44600000(1:100)  (checked)              |   0.057 0.069(1:100)(checked)
#100/200/5/4  0.067/0.769(1.97, 2.07) 0.06733333 0.36200000(51:150) (checked)              |   0.053 0.042(1:100)(checked)

#model2(censor = 0,0.5, lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.056 0.528       0.034     
#100/100/5/3  0.044 0.495       0.037          0.09583333 0.06333333
#100/100/5/4  0.075 0.488       0.049          0.1145833 0.0340000
#100/200/5/1  0.064 0.395(1000) 0.033
#100/200/5/2  0.049 0.427(1000) 0.046
#100/200/5/3  0.058 0.421       0.041
#100/200/5/4  0.094 0.405       0.050 
#----------------------------------------------------------------------------------------------------------------------------------------------------
#model2(censor = 0,0.5, lmdato=800,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.044 0.541       0.027
#100/200/5/4                    0.049

#model2(censor = 0,0.5, lmdato=800,lmdaby=5,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.063 0.418       0.053 
#100/200/5/4                    0.064

#model2(censor = 0,0.5, lmdato=800,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1                    0.027


#model2(lmdato=500,lmdaby=3,coef=4*sqrt(log(ngp)/n),lmdamtple=5,within = 0.5, between = 1e-4)
#100/100/5/1  0.02333333 0.78000000                                   0.016
#100/100/5/2  0.01833333 0.68800000   |
#100/100/5/2  0.03533333 0.52400000
#100/100/5/4  0.07064103 0.44400000   |    0.04666667 0.03900000
#100/200/5/1  0.03366667 0.71200000
#100/200/5/2  0.07633333 0.52400000   |    0.08333333 0.15200000
#100/200/5/3  0.055 0.424           
#100/200/5/4                                                          0.056

#model2(censor = 0,0.5, lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/200/6/1  0.045 0.338       0.022
#100/200/6/2  
#100/200/6/3  
#(X)#100/200/6/4 0.169 0.305                   0.070
#100/200/7/1  0.041 0.267       0.033
#100/200/7/2  0.087 0.229
#100/200/7/3  
#100/200/7/4  0.208 0.171
#100/200/10/1                   0.029
#100/200/10/2  
#100/200/10/3  
#(X)100/200/10/4  0.253 0.094                 0.048

#(X)model2(censor = 0,0.5, lmdato=500,lmdaby=1,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.057 0.633       0.030
#100/200/5/4                    0.005

#(X)model2(censor = 0,0.5, lmdato=300,lmdaby=1.5,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.077 0.597        0.036        
#(BAD)100/200/5/4                     0.01

#model2(censor = 0,0.5, lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.75, between = 0)
#100/100/5/1  0.04191667 0.56800000

#model2(censor = 0,0.5, lmdato=500,lmdaby=3,coef=1,lmdamtple=5,within = 0.25, between = 0)
#100/100/5/1  0.008916667 0.974000000
#100/100/5/4  0.07116667 0.56000000
#100/200/5/1  0.00475 0.92900




#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 0)
#100/200/5/4  0.086 0.31













#model2(lmdato=500,lmdaby=3,coef=4*sqrt(log(ngp)/n),lmdamtple=5,within = 0.25, between = 1e-4)
#100/100/5/1  0.03666667 0.60800000
#100/100/5/4  0.07966667 0.42000000
#100/200/5/1  0.02666667 0.55600000
#100/200/5/2  0.06166667 0.43200000











#model2(lmdato=500,lmdaby=3,coef=3*sqrt(log(ngp)/n),lmdamtple=5,within = 0.25, between = 1e-3)

#100/100/5/1  0.03133333 0.58800000(ALL)
#100/100/5/4  0.07233333 0.42800000(ALL)
#100/200/5/1  0.05166667 0.48800000(ALL)















#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 1e-3)
#100/100/5/1  0.032/                  0.05993333 0.42680000(ALL)
#100/100/5/2                          0.04957143 0.42228571(ALL)
#100/100/5/3                          0.07366667 0.39600000(middle)
#100/100/5/4  0.052/                  0.09608333 0.38650000(ALL)
#100/200/5/1  0.029/                  0.06083333 0.30200000(101:200)
#100/200/5/2                          0.06995833 0.32500000(ALL)
#100/200/5/3    
#100/200/5/4  0.057
















#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 1e-2)
#100/100/5/1  0.06333333 0.42880000(ALL)
#100/100/5/3  0.06744444 0.49066667(ALL)
#100/200/5/1  0.07783333 0.30160000(ALL)

#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.25, between = 1e-3)
#100/100/5/1  0.059 0.423(ALL)
#100/200/5/1  0.071 0.276(ALL)


#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.1, between = 1e-4)
#100/100/5/1  0.047 0.411(ALL)
#100/100/5/2  0.071 0.399(ALL)
#100/100/5/3  0.064 0.368(ALL) 
#100/100/5/4  0.099 0.371(ALL)
#100/200/5/1  0.096 0.283(ALL)
#100/200/5/2  0.090 0.308(ALL)
#100/200/5/2  0.128 0.295(ALL)

#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.1, between = 0)
#100/100/5/1  0.05016667 0.40680000(ALL)
#100/100/5/3            
#100/200/5/1  

#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.01, between = 0)
#100/100/5/1  0.026/            0.060 0.432
#100/100/5/3  0.039/            
#100/200/5/1  0.027/                 

#model2(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5,within = 0.01, between = 0.001)
#100/100/5/1  0.032/0.782       0.070/0.413
#100/100/5/3  0.050/            
#100/200/5/1  /                 0.068/0.270

#model3(lmdato=500,lmdaby=3,coef=0.5,lmdamtple=5)
#100/100/5/1  0.047/0.3125
#100/100/5/3  0.044/
#100/100/5/3  0.051/
#100/100/5/4  0.063/
#100/200/5/1  0.049/
#100/200/5/2  0.059/
#100/200/5/3  0.044/
#100/200/5/4  0.066/























#(num_per_batch=200,lmdato=500,lmdaby=3,coef=0.3218949,lmdamtple=1)
#100/100/5/1  0.035 (0.52, 1.73)   0.05416667 0.32800000        0.07523333 0.33360000
#100/100/5/2  0.041 (0.96, 1.89)   0.04083333 0.39000000  
#100/100/5/3  0.049 (1.46, 2.03)   0.1251667 0.3820000
#100/100/5/4  0.062 (1.92, 2.23)
#type1: (num_per_batch=200,lmdato=500,lmdaby=3,coef=0.3218949,lmdamtple=5)
#100/100/5/1  0.036 (0.49, 1.72)   0.08333333 0.34200000        0.0703 0.3404
#100/100/5/2  0.047 (0.95, 1.94)                                0.0596 0.3928
#100/100/5/3  0.049 (1.48, 2.03)   0.0925 0.3960                0.1063 0.3924
#100/100/5/4  0.047 (2.11, 1.95)   0.1332 0.3504
#100/200/5/1  0.033 (0.49, 1.74)
#100/200/5/2  0.039 (1.02, 1.77)
#100/200/5/3  0.047 (1.49, 1.96)
#100/200/5/4  0.060 (1.94, 2.17)






fdrtest<-function(n, ngp, s, K, ncp){
  ALLstat= cbind( matrix(rchisq(n*(ngp-s), K, 0), ncol = ngp-s),  matrix(rchisq(n*s, K, ncp), ncol = s)  )
  fdrsummary(ALLstat, ngp, s, K, 0.05, 0)
}  

if(FALSE){  
  #N=chi2absnorm[1,]; numgp=ngp;  alpha=0.05; otconservative=1
  #fdr(N, numgp,  alpha, H0, H1, Notconservative)
  error = c(); power = c(); tcur=c();l=c()
  n=500; ngp=100; s=5; K=1; ncp=1:25;
  for(i in ncp){
    print(i)
    l=c(l, i)
    tmp=fdrtest(n, ngp, s, K, i)
    error = c(error, tmp[1])
    power = c(power, tmp[2])
    tcur = c(tcur, tmp[3])
  }
  df = data.frame(l, error, power, tcur)
}

fdrtest1<-function(ngp){
  statistics = abs(rnorm(ngp,0,1))
  t=sqrt(2*log(ngp))
  print(length(statistics[statistics>t])/ngp)
}  