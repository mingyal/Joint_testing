

fdr <- function(N, numgp,  alpha, H0, H1, Notconservative){
  tset = rev(sort(N)); 
  tset = tset[  tset<=sqrt(  2*log(numgp)-2*log(log(numgp))*Notconservative  )   ]
  tcur = sqrt((2)*log(numgp))
  for(candidate in tset){
    #print(candidate)
    #print((2*(numgp)*(  1-pnorm(candidate, 0, 1)  )/(  max(  sum(N>=candidate) , 1  )  )) )
    if(    (2*(numgp)*(  1-pnorm(candidate, 0, 1)  )/(  max(  sum(N>=candidate) , 1  )  )) <=alpha   ){
      
      tcur = candidate
    }
  }
   # if(tcur>sqrt(  2*log(numgp)-2*log(log(numgp))*Notconservative    )){
   #   tcur = sqrt(2*log(numgp))
   # }
  #print(tcur)
  FDR =  sum(N[H0]>tcur)  /  max(  sum(N>tcur), 1  )  
  Power = sum(N[H1]>tcur)  /  max(length(H1), 1)
  return(c(FDR, Power, tcur))
}