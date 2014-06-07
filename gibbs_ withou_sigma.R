library(HyperbolicDist)
library(quantreg)
library(MASS)
simulate_gibbs<-function(n=100){
x<<-rep(1,n)
x<<-cbind(x,rnorm(n),rnorm(n))
#### data:x
epsilon1<-rnorm(n)
epsilon2<-rt(n,3)
epsilon3<-(x[,2]+1)*rnorm(n)
###three epsilon,simulate from a normal,t,heterosceadastic normal
y1<<-x%*%c(1,1,1)+epsilon1
y2<<-x%*%c(1,1,1)+epsilon2
y3<<-x%*%c(1,1,1)+epsilon3
}


### test

t<-c(0.1,0.5,0.9)

mean_fit1<-0
mean_fit2<-0
mean_fit3<-0
RMSE_fit1<-0
RMSE_fit2<-0
RMSE_fit3<-0
cof<-matrix(rep(1,9),3,3)
for (i in 1:1000){
simulate_gibbs()
fit<-rq(y1~x[,2]+x[,3],tau=t)
fit$coefficients/1000+mean_fit1->mean_fit1
(fit$coefficients-cof)*(fit$coefficients-cof)/1000+
  RMSE_fit1->RMSE_fit1
### quantreg for y1
fit<-rq(y2~x[,2]+x[,3],tau=t)
fit$coefficients/1000+mean_fit2->mean_fit2
(fit$coefficients-cof)*(fit$coefficients-cof)/1000+
  RMSE_fit2->RMSE_fit2
###quantreg for y2
fit<-rq(y3~x[,2]+x[,3],tau=t)
fit$coefficients/1000+mean_fit3->mean_fit3
(fit$coefficients-cof)*(fit$coefficients-cof)/1000+
  RMSE_fit3->RMSE_fit3
###quantreg for y3
}
mean_fit1
mean_fit2
mean_fit3
sqrt(RMSE_fit1)
sqrt(RMSE_fit2)
sqrt(RMSE_fit3)

### question1

gibbs_quantreg<-function(x,y,p,n=100){
  if (!dim(y)[2]==1){
    stop("Wrong in y!")
  }
  if (!dim(y)[1]==dim(x)[1]){
    stop("Different length of x and y!")
  }
  theta<-(1-2*p)/p/(1-p)
  tau2<-2/p/(1-p)
  psi=2+theta^2/tau2
  data_beta<-data.frame(0,0,0)
  names(data_beta)<-c("b1","b2","b3")
  z<-0
  b0=c(0,0,0)
  B0<-matrix(c(100,0,0,0,100,0,0,0,100),3,3)
  B0_solve<-solve(B0)
  data_beta[1,]<-c(0,0,0)
    for (i in 2:12000){
      chi<-(y-x%*%t(data_beta[i-1,]))^2/tau2
      for (j in 1:n){
      z[j]<-rgig(1,c(1/2,chi[j],psi))
      }
      sum<-t(x)%*%(x/cbind(z,z,z))    
      B<-solve(sum/tau2+B0_solve)
      a2<-(y/z-theta)/tau2
      b<-B%*%(t(x)%*%a2)
      
      data_beta[i,]<-mvrnorm(1,b,B)
    
  }
  data_beta
}

simulate_gibbs()
gibbs_quantreg(x,y1,0.9)->beta1
colMeans(beta1[2000:12000,])






















