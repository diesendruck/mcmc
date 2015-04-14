#setwd("D:\\Courses\\Spring15\\MonteCarloMthd\\Homework\\hw4");
library(ElemStatLearn);
library(gtools)
library(coda)
library(truncnorm)
data(galaxy)
#hist(galaxy$velocity,breaks=80,freq=F)
#lines(density(galaxy$velocity))

x = (galaxy$velocity-min(galaxy$velocity))/(max(galaxy$velocity)-min(galaxy$velocity));
n = length(x);
eps=1e-6;
alpha=1;tau=1;lambda=5;a=5;b=0.05;c=1;

rTruncPois=function(lambda){
  r=rpois(1,lambda);
  while(r==0){
    r=rpois(1,lambda);
  }
  return(r);
}

jointpost=function(x,J,mu,w,sigma2){
  #print("calculating posterior with param:");
  #print(list(J=J,mu=mu,w=w,sigma2=sigma2));
  n = length(x);
  
  # Posterior (Dirichlet and Mu's to follow)
  po=-lambda+J*log(lambda)-log(factorial(J))+log(sigma2)*(1-a)-b/sigma2-log(sigma2)*(n/2);
  # Dirichlet and Mu's
  for(j in 1:J){
    po=po+log(w[j]+eps)*(alpha-1)-mu[j]^2/2/tau;
  }
  for(i in 1:n){
    po=po+log(sum(w*exp(-(x[i]-mu)^2/2/sigma2)));
  }
  
  return(po);
}

mixture = function(x,J,mu,w,sigma2){
  return(sum(w*exp(-(x-mu)^2/2/sigma2)));
}

split = function(x,J,mu,w,sigma2){
  # Original
  ori = list(J=J,mu=mu,w=w,sigma2=sigma2);
  jointpost_last=jointpost(x,J,mu,w,sigma2);
  
  # J star
  js=sample.int(J,1);
  sortedmu=sort(mu,index.return=T);
  
  # Position of mu in sorted array
  ijs=match(mu[js],sortedmu$x);
  if(ijs==1){
    # Mu left
    mul=-99;
  }else{
    mul=sortedmu$x[ijs-1];
  }
  if(ijs==J){
    # Mu right
    mur=99;
  }else{
    mur=sortedmu$x[ijs+1];
  }
  
  # Sample new mu that's no farther than the nearest neighboring mu.
  u=rtruncnorm(1,a=0,b=min(mur-mu[js],mu[js]-mul),mean=0,sd=c);
  print("truncated normal")
  print(u)
  
  # To generate weights for new mu's (u).
  v=rbeta(1,1,1);
  
  # Densities of new u and new v.
  qu=ptruncnorm(u,a=0,b=min(mur-mu[js],mu[js]-mul),mean=0,sd=c);
  qv=pbeta(v,1,1);
  
  # Adds new mu at end of mu array.
  mu=c(mu,mu[js]+u);
  # Replaces the deceased mu[js] with the other new mu.
  mu[js]=mu[js]-u;
  
  # Weights updated similarly.
  w=c(w,(1-v)*w[js]); # Associated with mu[js]+u
  w[js]=w[js]*v; # Associated with mu[js]-u
  
  # Proposal
  rho=exp(jointpost(x,J+1,mu,w,sigma2)-jointpost_last);
  if(J==1){rho=rho*2;} # Pr(split|J=1) = 1; Pr(split|J>1) = 1/2
  rho=rho*(J/(J+1)); # TODO: Should this be rho*(J/J)?
  rho=rho/qu/qv;
  rho=rho*2*w[js];
  J=J+1;
  acc = min(1,rho);
  
  if(runif(1,0,1)<acc){
    return(list(J=J,mu=mu,w=w,sigma2=sigma2));
  }else{return(ori);}
}

merge = function(x,J,mu,w,sigma2){
  ori = list(J=J,mu=mu,w=w,sigma2=sigma2);
  # Storing the current state, in case we don't accept.
  jointpost_last=jointpost(x,J,mu,w,sigma2);
  sortedmu = sort(mu,index.return=T);
  j1 = sample.int((J-1),1);
  if(j1==1){
    mul=-99;
  }else{
    mul = sortedmu$x[j1-1];
  }
  if(j1==J-1){
    mur=99;
  }else{
    mur=sortedmu$x[j1+2];
  }
  j2 = sortedmu$ix[j1+1];
  j1 = sortedmu$ix[j1];
  
  # Half the distance between merge candidates.
  u = abs((mu[j1]-mu[j2])/2);
  mu = c(mu,mean(mu[j1],mu[j2]));
  mu = mu[-c(j1,j2)];
  w = c(w,(w[j1]+w[j2]));
  v = w[j1]/(w[j1]+w[j2])
  w = w[-c(j1,j2)];
  # Densities of new u and new v.
  qu=ptruncnorm(u,a=0,b=min(mur-mu[J-1],mu[J-1]-mul),mean=0,sd=c);
  qv=pbeta(v,1,1);
  rho=exp(jointpost(x,J-1,mu,w,sigma2)-jointpost_last);
  if(J==2){rho=rho/2;}
  rho=rho*(J/(J-1)); # TODO: Is this (J-1)/(J-1)?
  rho=rho*qu*qv;
  rho=rho/2/w[J-1];
  J=J-1;
  acc = min(1,rho);
  if(runif(1,0,1)<acc){
    return(list(J=J,mu=mu,w=w,sigma2=sigma2));
  }else{return(ori);}
}

sample.s <- function(J,w,mu,sig2)
{ ## sample s[i] from p(s[i] | ...)
  n <- length(x)
  sd <- sqrt(sig2)
  s <- rep(0,n) # initialize
  for(i in 1:n){
    pr <- w*dnorm(x[i],m=mu,sd=sd)
    s[i] <- sample(1:J,1,replace=T,prob=pr)
  }
  return(s)
}
sample.mu <- function(J,s,w,sig2)
{ ## sample mu[j]
  mu <- rep(0,J) # initialize
  for (j in 1:J){
    Aj <- which(s==j)
    nj <- length(Aj)
    if (nj==0){
      m <- 0; V <- tau
    } else {
      xbar <- mean(x[Aj])
      V <- 1/(1/tau + nj/sig2)
      m <- V*xbar*nj/sig2
    }
    mu[j] <- rnorm(1,m=m,sd=sqrt(V))
  }
  return(mu)
}
sample.w <- function(J,s,sig2)
{## sample w
  a <- rep(0,J) # initialize
  for(j in 1:J){
    Aj <- which(s==j)
    nj <- length(Aj)
    a[j] <- alpha+nj
  }
  w <- rdirichlet(1, a)
  w <- c(w)
  return(w)
}
sample.sig <- function(J,s,w,mu)
{## sample sig2
  S2 <- sum((x-mu[s])^2)
  a1 <- a+n/2
  b1 <- b+S2/2
  sig2 <- 1.0/rgamma(1,shape=a1,rate=b1)
  return(sig2)
}

rjmcmc = function(iter){
  # Initialization
  J=rTruncPois(lambda);
  alphavec=rep(alpha,J);
  w=rdirichlet(1,alphavec);
  mu=rnorm(J,0,sqrt(tau));
  sigma2=1/(rgamma(1,shape=a,rate=b));
  
  # Iteration
  Jlist = rep(0,iter);
  siglist = rep(0,iter);
  f1 = rep(0,iter);
  f2 = rep(0,iter);
  f3 = rep(0,iter);
  for (i in 1:iter){
    print("Iteration:")
    print(i);
    # Generating s,mu,w,sigma2
    s=sample.s(J,w,mu,sigma2);
    mu = sample.mu(J,s,w,sigma2);
    w = sample.w(J,s,sigma2);
    sigma2 = sample.sig(J,s,w,mu);
    # RJ Move
    tmp=runif(1,0,1);
    if(J==1 || tmp<0.5){
      th = split(x,J,mu,w,sigma2);
      J=th$J;mu=th$mu;w=th$w;sigma2=th$sigma2;
    }
    else{
      th = merge(x,J,mu,w,sigma2);
      J=th$J;mu=th$mu;w=th$w;sigma2=th$sigma2;
    }
    Jlist[i]=th$J;
    siglist[i]=th$sigma2;
    f1[i]=mixture(0.1,J,mu,w,sigma2);
    f2[i]=mixture(0.5,J,mu,w,sigma2);
    f3[i]=mixture(0.9,J,mu,w,sigma2);
  }
  plot(c(1:iter),Jlist,'l',main="Trajectory of J");
  summary(mcmc(data.frame(Jlist,siglist)));
  xgrid=seq(0,1,0.01);
  f=xgrid;
  for(i in 1:length(xgrid)){
    f[i]=mixture(xgrid[i],J,mu,w,sigma2);
  }
  plot(c(1:length(f)),f,'l',main='Estimation for f');
  return(mcmc(data.frame(Jlist,siglist,f1,f2,f3)));
}

mcmc1 = rjmcmc(500);
#mcmc2 = rjmcmc(200);
#rj.list = mcmc.list(list(mcmc1,mcmc2));
#gelman.diag(rj.list);
#gelman.plot(rj.list);
traceplot(mcmc1);
geweke.diag(mcmc1[200:500,]);
