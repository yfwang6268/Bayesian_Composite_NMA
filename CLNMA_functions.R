############################################################################################
###########Functions for Composite Likelihood Network Meta-Analysis (CLNMA)#################
###########Author: Rui Duan  Last Updated: 08/23/2020 ######################################

library(mvtnorm)



####CL functions#####

RElik = function (t,a, dataout,narm){
  #t----tau^2, for unequal t, there are (narm choose 2) different t for each outcome.
  #a----mu
  #narm----how many arms in total.
  #note : the reference drug has to be indicated by the largest number (narm)
  #narm = max(c(dataout$t1,dataout$t2))
  l=0
  penalty = 0
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      mu = a[i] - a[j]
      indext = i*(narm-1) - i*(i-1)/2+j-i
      l.temp = -0.5*sum(log(t[indext]+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t[indext]+data$sd^2))
      l = l+l.temp
      penalty = penalty + log(sum(t[indext]+data$sd^2)^{-1})
    }
  }
  for(k in 1:(narm-1)){
    data = subset(dataout,dataout$t1==i&dataout$t2==narm)
    mu =a[k]
    l.temp = -0.5*sum(log(t[k]+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t[k]+data$sd^2))
    l = l+l.temp
    penalty = penalty + log(sum(t[k]+data$sd^2)^{-1})
  }
reml = l-0.5*(penalty)
  return(-reml)
}
simu = function(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2){
  dataout = gendata(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2)
  out1 = tryCatch(CLNMA(dataout),error=function(e) rep(NA,20))
  out2 = tryCatch(CLNMA.equal.tau(dataout),error=function(e) rep(NA,12))
  return(c(out1,out2))
}


#Format dataset
select = function(i,dataout){
  #Function to select all compares to i
  data1 = subset(dataout,dataout$t2==i)
  data2 = subset(dataout,dataout$t1==i)
  data2$outcome = -data2$outcome
  data = rbind(data1,data2)
  return(data)
}


#estimate effect sizes given tau
get.mu = function (t,dataout,narm){
  Y = W = rep(0,(narm-1))
  for (i in 1:(narm-1)){
    datatemp = select(i,dataout)
    y = datatemp$outcome
    w = 1/(datatemp$sd^2+t[i])
    Y[i] = sum(y*w)
    W[i] = sum(w)
  }
  
  A = diag(W)
  
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      indext = i*(narm-1) - i*(i-1)/2+j-i
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      A[i,j] = A[j,i] =-sum(1/(data$sd^2+t[indext]))
    }
  }
  mu = solve(-A, Y)
  return(mu)
}

#Iterative algorithm for estimating all parameters
Iterate.REML = function (t0,tol=10^(-6),maxiter=500,dataout,narm){
  delta = 10
  m=0
  while(delta>tol & m<maxiter){
    A = get.mu(t0,dataout,narm) 
    re = function(t){RElik(t,a=A,dataout,narm)}
    t1 = optim(t0, re, method = 'L-BFGS-B', lower = rep(0,length(t0)))$par
    delta = max(abs(t1-t0))
    t0 = t1
    #print(t0)
    m=m+1
  }
  {if(m>=maxiter){t1 =NA}}
  return(t1)
}


#Calculate  Hessian
Hessian = function (a,t,dataout,narm){
  #narm = max(c(dataout$t1,dataout$t2))
  H_mu = matrix(0,nrow = narm-1, ncol = narm-1)#submtrix for mu
  #H_tau =matrix(0,nrow = narm*(narm-1)/2, ncol = narm*(narm-1)/2)#submtrix for tau
  H_tm =matrix(0,nrow = narm-1, ncol = narm*(narm-1)/2)
  T1 = T2 = temp_mu =temp_tau= temp_mutau= rep(0,narm*(narm-1)/2)
  
  for (i in 1:(narm-1)){
      indext = i
      data = subset(dataout,dataout$t1==i&dataout$t2==narm)
      y = data$outcome
      vv = data$sd^2
      mu = a[i]
      q1 = sum(1/(vv+t[indext]))
      q2 = sum(1/(vv+t[indext])^2)
      q3 = sum((y-mu)/(vv+t[indext])^2)
      q4 = sum((y-mu)^2/(vv+t[indext])^3)
      temp_mu[i] = q1
      temp_mutau[i] = q3
      temp_tau[i] = q2/2-q4
      T1[i] = i
      T2[i] = narm
  }
   
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      indext = i*(narm-1) - i*(i-1)/2+j-i
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      y = data$outcome
      vv = data$sd^2
      mu = a[i]-a[j]
      q1 = sum(1/(vv+t[indext]))
      q2 = sum(1/(vv+t[indext])^2)
      q3 = sum((y-mu)/(vv+t[indext])^2)
      q4 = sum((y-mu)^2/(vv+t[indext])^3)
      temp_mu[indext] = q1
      temp_mutau[indext] = q3
      temp_tau[indext] = q2/2-q4
      H_mu[i,j] = H_mu[j,i] = q1
      T1[indext] = i
      T2[indext] = j
    }
  }
  
  
  for (i in 1:(narm-1)){
    index1 = as.numeric(T1 ==i)
    index2 = as.numeric(T2 ==i)
    indmu = index1+index2
    H_mu[i,i] = -sum(indmu*temp_mu)
    indtm = -index1+index2
    H_tm[i,] = indtm*temp_mutau
  }
  H_tau = diag(temp_tau)
  H = rbind(cbind(H_mu,H_tm),cbind(t(H_tm),H_tau))
  return(H)
}

#Meat in sandwich estimator
Scoresquare = function(a1,t1,a2,t2,dataout,narm){
  #narm = max(c(dataout$t1,dataout$t2))

  I1 = rep(0,narm-1)
  I2 = rep(0,narm*(narm-1)/2)
  SS = matrix(0,nrow =2*(length(I1)+length(I2)),ncol = 2*(length(I1)+length(I2)))
  for (i in 1:dim(dataout)[1]){
    if(dataout$t2[i]==narm){
    indext = dataout$t1[i]
    m1 = a1[indext]#mean outcome1
    m2 = a2[indext]#mean outcome2
    I1[indext] =1
    I2[indext] =1
    }
    else{
      ind1 = dataout$t1[i]
      ind2 = dataout$t2[i]
      m1 = a1[ind1] - a1[ind2]
      m2 = a2[ind1] - a2[ind2]
      I1[ind1] = 1
      I1[ind2] = -1 
      indext = ind1*(narm-1) - ind1*(ind2-1)/2+ind2-ind1
      I2[indext] = 1
    }
    S_mu1 = I1*((dataout$outcome1[i]-m1)/(t1[indext]+dataout$sd1[i]^2))
    S_tau1 =I2*((dataout$outcome1[i]-m1)^2/2/(t1[indext]+dataout$sd1[i]^2)^2-1/2/(t1[indext]+dataout$sd1[i]^2))
    S_mu2 = I1*((dataout$outcome2[i]-m2)/(t2[indext]+dataout$sd2[i]^2))
    S_tau2 =I2*((dataout$outcome2[i]-m2)^2/2/(t2[indext]+dataout$sd2[i]^2)^2-1/2/(t2[indext]+dataout$sd2[i]^2))
    S = c(S_mu1,S_tau1,S_mu2,S_tau2)
    temp = t(t(S))%*%S
    SS = temp+SS
  }
return(SS)
}


#CLNMA with unequal variances
CLNMA = function(dataout){
  narm = max(c(dataout$t1,dataout$t2))
  dataout1 = dataout[,c("ID","t1","t2","outcome1","sd1")]
  colnames(dataout1)[4:5] = c("outcome","sd")
  dataout2 = dataout[,c("ID","t1","t2","outcome2","sd2")]
  colnames(dataout2)[4:5] = c("outcome","sd")
  
  #get estimations
  tau1=Iterate.REML(rep(0.5,narm*(narm-1)/2),tol=10^(-5),maxiter=500,dataout1,narm)
  mu1=get.mu(tau1,dataout1,narm)
  
  tau2=Iterate.REML(rep(0.5,narm*(narm-1)/2),tol=10^(-5),maxiter=500,dataout2,narm)
  mu2=get.mu(tau2,dataout2,narm)
  
  #getHessian
  H1 = Hessian(mu1,tau1,dataout1,narm)
  H2 = Hessian(mu2,tau2,dataout2,narm)
  Zero = matrix(0,nrow = dim(H1)[1], ncol = dim(H1)[1])
  H = rbind(cbind(H1,Zero),cbind(Zero,H2))
  #get SS
  SS = Scoresquare(mu1,tau1,mu2,tau2,dataout,narm)
  #get variance
  V = solve(H)%*%SS%*%solve(H)
  
  return(c(mu1,mu2,tau1,tau2,diag(V)))
}

#Generate data in simulation 
generate.data = function(mu1,mu2,betweenv,rho_w,ss1,ss2){
  #library(MASS)
s1 =abs(rnorm(2,ss1,0.1))
s2 =abs(rnorm(2,ss2,0.1))
s = c(s1,s2)^2
withinv = diag(s)
for(i in 1:3){
  for(j in (i+1):4){
    withinv[i,j] = withinv[j,i]=rho_w[i,j]*sqrt(s[i]*s[j])
  }
}
y = MASS::mvrnorm(1,c(mu1,mu2),withinv+betweenv)
return(c(y,s1,s2))
}

gendata = function(nab,nac,nbc,nabc,mu1,mu2,betweenv,rho_w,ss1,ss2){
  n = nab+nac+nbc+nabc
  data = replicate(n,generate.data(mu1,mu2,betweenv,rho_w,ss1,ss2))
  
  data = data.frame(t(data))
  colnames(data) = c("Y_BA_1","Y_CA_1","Y_BA_2","Y_CA_2",
                     "S_BA_1","S_CA_1","S_BA_2","S_CA_2")
  
  data$Y_BC_1 = data$Y_BA_1-data$Y_CA_1
  data$Y_BC_2 = data$Y_BA_2-data$Y_CA_2
  data$S_BC_1 = sqrt(data$S_BA_1^2+data$S_CA_1^2-2*rho_w[1,2]*data$S_CA_1*data$S_BA_1)
  data$S_BC_2 = sqrt(data$S_BA_2^2+data$S_CA_2^2-2*rho_w[3,4]*data$S_CA_2*data$S_BA_2)
  #BAarm
  dataBA = data.frame(ID = 1:nab, t1 = rep(1,nab),t2 = rep(3,nab))
  dataBA$outcome1 = data[1:nab,"Y_BA_1"]
  dataBA$sd1 = data[1:nab,"S_BA_1"]
  dataBA$outcome2 = data[1:nab,"Y_BA_2"]
  dataBA$sd2 = data[1:nab,"S_BA_2"]
  
  #CAarm
  dataCA = data.frame(ID = nab+1:nac, t1 = rep(2,nac),t2 = rep(3,nac))
  dataCA$outcome1 = data[nab+1:nac,"Y_CA_1"]
  dataCA$sd1 = data[nab+1:nac,"S_CA_1"]
  dataCA$outcome2 = data[nab+1:nac,"Y_CA_2"]
  dataCA$sd2 = data[nab+1:nac,"S_CA_2"]
  
  #BCarm
  dataBC = data.frame(ID = nab+nac+1:nbc, t1 = rep(1,nbc),t2 = rep(2,nbc))
  dataBC$outcome1 = data[nab+nac+1:nbc,"Y_BC_1"]
  dataBC$sd1 = data[nab+nac+1:nbc,"S_BC_1"]
  dataBC$outcome2 = data[nab+nac+1:nbc,"Y_BC_2"]
  dataBC$sd2 = data[nab+nac+1:nbc,"S_BC_2"]
  
  #ABCarm
  dataABC1 = data.frame(ID = nab+nac+nbc+1:nabc, t1 = rep(1,nabc),t2 = rep(3,nabc))
  dataABC1$outcome1 = data[nab+nac+nbc+1:nabc,"Y_BA_1"]
  dataABC1$sd1 = data[nab+nac+nbc+1:nabc,"S_BA_1"]
  dataABC1$outcome2 = data[nab+nac+nbc+1:nabc,"Y_BA_2"]
  dataABC1$sd2 = data[nab+nac+nbc+1:nabc,"S_BA_2"]
  
  dataABC2 = data.frame(ID = nab+nac+nbc+1:nabc, t1 = rep(2,nabc),t2 = rep(3,nabc))
  dataABC2$outcome1 = data[nab+nac+nbc+1:nabc,"Y_CA_1"]
  dataABC2$sd1 = data[nab+nac+nbc+1:nabc,"S_CA_1"]
  dataABC2$outcome2 = data[nab+nac+nbc+1:nabc,"Y_CA_2"]
  dataABC2$sd2 = data[nab+nac+nbc+1:nabc,"S_CA_2"]
  
  
  dataABC3 = data.frame(ID = nab+nac+nbc+1:nabc, t1 = rep(1,nabc),t2 = rep(2,nabc))
  dataABC3$outcome1 = data[nab+nac+nbc+1:nabc,"Y_BC_1"]
  dataABC3$sd1 = data[nab+nac+nbc+1:nabc,"S_BC_1"]
  dataABC3$outcome2 = data[nab+nac+nbc+1:nabc,"Y_BC_2"]
  dataABC3$sd2 = data[nab+nac+nbc+1:nabc,"S_BC_2"]
  
  dataout = rbind(dataBA,dataCA,dataBC,dataABC1,dataABC2,dataABC3)
  return(dataout)
}

######equal between study variance########

RElik.equal.tau = function (t,a, dataout,narm){
  #t----tau^2, for unequal t, there are (narm choose 2) different t for each outcome.
  #a----mu
  #narm----how many arms in total.
  #note : the reference drug has to be indicated by the largest number (narm)
  l=0
  #penalty = 0
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      mu = a[i] - a[j]
      #indext = i*(narm-1) - i*(i-1)/2+j-i
      l.temp = -0.5*sum(log(t+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t+data$sd^2))
      l = l+l.temp
      #penalty = penalty 
    }
  }
  for(k in 1:(narm-1)){
    data = subset(dataout,dataout$t1==i&dataout$t2==narm)
    mu =a[k]
    l.temp = -0.5*sum(log(t+data$sd^2))-0.5*sum((data$outcome-mu)^2/(t+data$sd^2))
    l = l+l.temp
    #penalty = penalty 
  }
  reml = l-0.5*log(sum(t+dataout$sd^2)^{-1})
  return(-reml)
}
get.mu.equal.tau = function (t,dataout,narm){
#narm = max(c(dataout$t1,dataout$t2))
  Y = W = rep(0,narm-1)
  for (i in 1:(narm-1)){
    datatemp = select(i,dataout)
    y = datatemp$outcome
    w = 1/(datatemp$sd^2+t)
    Y[i] = sum(y*w)
    W[i] = sum(w)
  }
  
  A = diag(W)
  
  for (j in 1:(narm-2)){
    for (k in (j+1):(narm-1)){
      #indext = i*(narm-1) - i*(i-1)/2+j-
        data = subset(dataout,dataout$t1==j&dataout$t2==k)
        A[j,k] = A[k,j] =-sum(1/(data$sd^2+t))
    }
  }
  mu = solve(-A, Y)
  return(mu)
}

Iterate.REML.equal.tau = function (t0,tol=10^(-5),maxiter=50,dataout,narm){
  delta = 10
  m=0
  while(delta>tol & m<maxiter){
    A = get.mu.equal.tau(t0,dataout,narm) 
    re = function(t){RElik.equal.tau(t,a=A,dataout,narm)}
    t1 = nlminb(t0, re, lower = 0, upper= Inf)$par
    delta = max(abs(t1-t0))
    t0 = t1
    #print(t0)
    m=m+1
  }
  {if(m>=maxiter){t1 =NA}}
  return(t1)
}

Hessian.equal.tau = function (a,t,dataout,narm){
  H_mu = matrix(0,nrow = narm-1, ncol = narm-1)#submtrix for mu
  H_tm =matrix(0,nrow = narm-1, ncol = 1)
  T1 = T2 = temp_mu =temp_tau= temp_mutau= rep(0,narm*(narm-1)/2)
  for (i in 1:(narm-1)){
    data = subset(dataout,dataout$t1==i&dataout$t2==narm)
    y = data$outcome
    vv = data$sd^2
    mu = a[i]
    q1 = sum(1/(vv+t))
    q2 = sum(1/(vv+t)^2)
    q3 = sum((y-mu)/(vv+t)^2)
    q4 = sum((y-mu)^2/(vv+t)^3)
    temp_mu[i] = q1
    temp_mutau[i] = q3
    temp_tau[i] = q2/2-q4
    T1[i] = i
    T2[i] = narm
  }
  for (i in 1:(narm-2)){
    for (j in (i+1):(narm-1)){
      indext = i*(narm-1) - i*(i-1)/2+j-i
      data = subset(dataout,dataout$t1==i&dataout$t2==j)
      y = data$outcome
      vv = data$sd^2
      mu = a[i]-a[j]
      q1 = sum(1/(vv+t))
      q2 = sum(1/(vv+t)^2)
      q3 = sum((y-mu)/(vv+t)^2)
      q4 = sum((y-mu)^2/(vv+t)^3)
      temp_mu[indext] = q1
      temp_mutau[indext] = q3
      temp_tau[indext] = q2/2-q4
      H_mu[i,j] = H_mu[j,i] = q1
      T1[indext] = i
      T2[indext] = j
    }}
    for (i in 1:(narm-1)){
      index1 = as.numeric(T1 ==i)
      index2 = as.numeric(T2 ==i)
      ind = index1+index2
      H_mu[i,i] = -sum(ind*temp_mu)
      H_tm[i,] = -sum(ind*temp_mutau)
    }
      H_tau = sum(temp_tau)
      H = rbind(cbind(H_mu,H_tm),c(t(H_tm),H_tau))
  return(H)
}


#for two outcomes
SS.equal.tau = function(a1,a2,t1,t2,dataout,narm){
  #narm = max(c(dataout$t1,dataout$t2))
  SS = matrix(0,nrow = narm*2,ncol=narm*2)
  for (i in 1:dim(dataout)[1]){
    I = rep(0,narm-1)
    {if (dataout$t2[i] == narm){
      m1 = a1[as.numeric(dataout$t1[i])]
      m2 = a2[as.numeric(dataout$t1[i])]
      I[as.numeric(dataout$t1[i])] =1
    }
      else {
        m1 = a1[as.numeric(dataout$t1[i])] - a1[as.numeric(dataout$t2[i])]
        m2 = a2[as.numeric(dataout$t1[i])] - a2[as.numeric(dataout$t2[i])]
        I[as.numeric(dataout$t1[i])] = 1
        I[as.numeric(dataout$t2[i])] = -1
      } } 
    
    S1 = ((dataout$outcome1[i]-m1)/(t1+dataout$sd1[i]^2))*I
    st1 = (dataout$outcome1[i]-m1)^2/2/(t1+dataout$sd1[i]^2)^2-1/2/(t1+dataout$sd1[i]^2)
    
    S2 = ((dataout$outcome2[i]-m2)/(t2+dataout$sd2[i]^2))*I
    st2 = (dataout$outcome2[i]-m2)^2/2/(t2+dataout$sd2[i]^2)^2-1/2/(t2+dataout$sd2[i]^2)
    
    S = c(S1,st1,S2,st2)
    #print(S)
    SS.temp = t(t(S))%*%S
    #print(SS.temp)
    SS = SS.temp+SS
  }
  return(SS)}


CLNMA.equal.tau = function(dataout){
  narm = max(c(dataout$t1,dataout$t2))
  dataout1 = dataout[,c("ID","t1","t2","outcome1","sd1")]
  colnames(dataout1)[4:5] = c("outcome","sd")
  dataout2 = dataout[,c("ID","t1","t2","outcome2","sd2")]
  colnames(dataout2)[4:5] = c("outcome","sd")
  
  #get estimations
  tau1=Iterate.REML.equal.tau(0.5,tol=10^(-5),maxiter=500,dataout1,narm)
  mu1=get.mu.equal.tau(tau1,dataout1,narm)
  
  tau2=Iterate.REML.equal.tau(0.5,tol=10^(-5),maxiter=500,dataout2,narm)
  mu2=get.mu.equal.tau(tau2,dataout2,narm)
  
  #getHessian
  H1 = Hessian.equal.tau(mu1,tau1,dataout1,narm)
  H2 = Hessian.equal.tau(mu2,tau2,dataout2,narm)
  Zero = matrix(0,nrow = dim(H1)[1], ncol = dim(H1)[1])
  H = rbind(cbind(H1,Zero),cbind(Zero,H2))
  #get SS
  SS = SS.equal.tau(mu1,mu2,tau1,tau2,dataout,narm)
  #get variance
  V = solve(H)%*%SS%*%solve(H)
  
  return(c(mu1,mu2,tau1,tau2,diag(V)))
}

