library(survival)
library(MASS)
library(matrixcalc)
library(openxlsx)
library(nloptr)
library(nlme)
library(joineRML)
data(heart.valve)
data(pbc)

data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$grad) & !is.na(heart.valve$lvmi) & !is.na(heart.valve$ef), ]
s=0
Ni=NULL
ID=NULL
ab=NULL
Idx=unique(hvd$num)
for (i in 1:length(Idx) ) {
  idx=which(hvd$num ==Idx[i])
  s=s+length(idx)
  Ni[i]=length(idx)
  ID[i]=s
  ab=append(ab,rep(i,length(idx)))
}
Hvd=hvd[ID,]
l1=Hvd[which(Hvd$hs=="Homograft"),]
l2=Hvd[which(Hvd$hs=="Stentless valve"),]
fit1 = survfit(Surv(l1$fuyrs,l1$status)~1)
fit2 = survfit(Surv(l2$fuyrs,l2$status)~1)
plot(fit1$time,fit1$surv, ylim=c(0,1),main="Kaplan-Meier estimate", xlab = "Time (years)", ylab = "Survival",type="l",col=4,lwd=2)
lines(fit2$time,fit2$surv, col=2,type="b",cex=0.2,lwd=2)
legend("bottomleft",inset=.05, c("Homograft","Stentless valve"),col=c(4,2),lty=c(1,2),cex=1,bty="n")

## two patients
pt1<-hvd[which(hvd$num==227),]## survive 
pt2<-hvd[which(hvd$num==250),]## died
plot(pt1$time,log(pt1$lvmi), ylim=c(3,6),xlab = "Time (years)", ylab = "log(lvmi)",type="l",col=4,lwd=2)
lines(pt2$time,log(pt2$lvmi), col=2,type="b",cex=0.2,lwd=2)

plot(pt1$time,log(pt1$grad), ylim=c(1,5),xlab = "Time (years)", ylab = "log(grad)",type="l",col=4,lwd=2)
lines(pt2$time,log(pt2$grad), col=2,type="b",cex=0.2,lwd=2)

plot(pt1$time,(pt1$ef), ylim=c(40,100),xlab = "Time (years)", ylab = "ef",type="l",col=4,lwd=2)
lines(pt2$time,(pt2$ef), col=2,type="b",cex=0.2,lwd=2)


## parameter estimation
J1=4
J=3
n=length(Idx)
dim(hvd)
Ni ## repeated observation
ID ## last observation of every patient
id=hvd$num 
## data
t<-hvd$time
X<-hvd[,c(2,3,19,25)]
ho=which(X$hs=="Homograft")
X$hs=0
X$hs[ho]=1
X=as.matrix(X)
Z<-matrix(0,nrow=sum(Ni),ncol=J1*J)
for (j in 1:J) {
  Z[,2*j-1]=rep(1,sum(Ni))
  Z[,2*j]=t
}
Y=cbind(hvd$log.lvmi,hvd$ef/10, hvd$log.grad)
Y<-vec(Y)
## h0
Xt=X[ID,]
delta<-hvd$status[ID]
ttoevent<-hvd$time[ID]
data.y <-as.data.frame(cbind(id=rep(rep(1:n, Ni),J), y.name= rep(1:J, sum(Ni)), y.t=rep(t, J), x.1=rep(X[,1],J), x.2=rep(X[,2],J), x.3=rep(X[,3],J),  x.4=rep(X[,4],J)))
data.y$y= Y
data.x <- data.y[!duplicated(data.y$id),]	
data.t= as.data.frame(cbind(id=c(1:n), time=ttoevent, delta=delta, x.1=data.x$x.1, x.2=data.x$x.2,x.3=data.x$x.3,x.4=data.x$x.4)) 

n.break=5
t.break=quantile(hvd$fuyrs,seq(0, 1,length=6))
k=4

alpha.ini= matrix(0, nrow=k,ncol=J)
beta.ini= rep(0,J1*J)
Sigma.ini<-0.5*diag(2*J)
b.ini=mvrnorm(n,rep(0,2*J),Sigma.ini)
sigma.e2.ini=1
coxfit1 = coxph(Surv(ttoevent,delta) ~ Xt[,1] + Xt[,2] + Xt[,3]+Xt[,4],model = TRUE)
gamma0.ini= as.vector(coxfit1$coef[1:4])
gamma.ini = as.vector(coxfit1$coef[5:10])
data.t= cbind(data.t, b.ini)

base <- numeric(n.break)
SS= nrow(data.t)
idRange <- 1:nrow(data.t) 
data.Y <- data.y[data.y$id %in% idRange,]
data.T <- data.t[data.t$id %in% idRange,]

## function
{
  GS <- function(){
    b.out <- matrix(rep(0,2*J*SS),ncol=2*J)
    rownames(b.out) <- idRange
    
    for (id in idRange){
      
      data <- data.Y[data.Y$id==id,]
      N <- nrow(data)
      X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
      time <- data.T[data.T$id==id,]$time
      delta <- data.T[data.T$id==id,]$delta
      
      t <- data[data$y.name==1,]$y.t
      n <- length(t)
      Zi <- matrix(c(rep(1,n),t),n,2)
      for (j in 2:J) {
        t <- data[data$y.name==j,]$y.t
        n <- length(t)
        z <- matrix(c(rep(1,n),t),n,2)
        Zi <- superMatrix(list(Zi,z))
      }
      
      b.est <- rep(0,2*J)
      kbi.drv <- rep(1,2*J)
      
      for (int in 1:n.break) {
        if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0[int] }
      }  
      
      GS.iter <- 0
      while (GS.iter < 10 & sum(kbi.drv^2) > 1e-5){
        b.est.old <- b.est
        
        y.name.ind.all <- data$y.name
        alpha.sel.left.all <- J1*(y.name.ind.all-1)+1
        beta.sel.all <- 2*y.name.ind.all-1
        
        alpha.matrix <- matrix(rep(0,N*J1),nrow=N)
        for (i in 1:N) alpha.matrix[i,] <- alpha[(alpha.sel.left.all[i]):(alpha.sel.left.all[i]+J1-1)]
        
        mu <- apply((X)*alpha.matrix,1,sum) + beta[beta.sel.all] + beta[beta.sel.all+1]*data$y.t + Zi%*%b.est.old
        
        eta <- X[1,]%*%gamma0 + t(gamma)%*%b.est.old
        
        kbi.drv <- t(data$y-mu)%*%Zi/sigma.e2 + delta*gamma - as.numeric(h.base*time*exp(eta))*gamma - t(b.est.old)%*%solve(Sigma)
        Hi <- -t(Zi)%*%Zi/sigma.e2 - as.numeric(h.base*time*exp(eta))*gamma%*%t(gamma) - solve(Sigma)
        
        b.est <- b.est.old - solve(Hi)%*%t(kbi.drv)
        
        GS.iter <- GS.iter+1
      }
      
      b.out[rownames(b.out)==as.character(id),] <- as.vector(b.est)     
    }
    
    #appr.lik <- sum(f-log(H)/2) 
    
    b.out
    #list(appr.lik=appr.lik,b.est=b.est,H=H,f=f)
  }
  negPL.fun <- function(alpha.est=alpha, beta.est=beta, gamma0.est=gamma0, gamma.est=gamma, h0.est=h0, sigma.e2.est=sigma.e2, Sigma.est=Sigma, b.est=b){
    PL.neg <- 0
    for (id in idRange){
      data <- data.Y[data.Y$id==id,]
      N <- nrow(data)
      X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
      time <- data.T[data.T$id==id,]$time
      delta <- data.T[data.T$id==id,]$delta
      
      t <- data[data$y.name==1,]$y.t
      n <- length(t)
      Zi <- matrix(c(rep(1,n),t),n,2)
      for (j in 2:J) {
        t <- data[data$y.name==j,]$y.t
        n <- length(t)
        z <- matrix(c(rep(1,n),t),n,2)
        Zi <- superMatrix(list(Zi,z))
      }
      mu <- sapply(1:N, function(i) c(X[i,]%*%alpha.est[(J1*(data$y.name[i]-1)+1):(J1*data$y.name[i])]) + beta.est[2*data$y.name[i]-1] + beta.est[2*data$y.name[i]]*data$y.t[i] + b.est[rownames(b.est)==as.character(id),2*data$y.name[i]-1] + b.est[rownames(b.est)==as.character(id),2*data$y.name[i]]*data$y.t[i])
      eta <- as.numeric(X[1,]%*%gamma0.est+gamma.est%*%b.est[rownames(b.est)==as.character(id),])
      
      for (int in 1:n.break) {
        if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0.est[int] }
      }      
      
      
      kbi.neg <- N*log(sigma.e2.est)/2 + c(t(data$y-mu)%*%(data$y-mu)/sigma.e2.est/2)-delta*log(h.base) - delta*eta + h.base*time*exp(eta) + log(det(Sigma.est)+0.0000001)/2 + t(b.est[rownames(b.est)==as.character(id),])%*%solve(Sigma.est)%*%b.est[rownames(b.est)==as.character(id),]/2 
      Hi <- -t(Zi)%*%Zi/sigma.e2.est - h.base*time*exp(eta)*gamma.est%*%t(gamma.est) - solve(Sigma.est)
      PLi.neg <- kbi.neg+log(det(Hi)+0.00000001)/2
      PL.neg <- PL.neg + PLi.neg
    }
    as.numeric(PL.neg)/length(idRange) + sqrt(2)*lambda*(sqrt(sum((gamma.est[1:2])^2))+sqrt(sum((gamma.est[3:4])^2))+sqrt(sum((gamma.est[5:6])^2))+sqrt(sum((gamma.est[7:8])^2))+sqrt(sum((gamma.est[9:10])^2))) + lambda2*sum(abs(Sigma.est[upper.tri(Sigma.est)]))
  }
  data.X <- data.Y[!duplicated(data.Y$id),]
  Z <- NULL
  for (id in idRange){
    data <- data.Y[data.Y$id==id,]
    N <- nrow(data)
    
    t <- data[data$y.name==1,]$y.t
    n <- length(t)
    Zi <- matrix(c(rep(1,n),t),n,2)
    for (j in 2:J) {
      t <- data[data$y.name==j,]$y.t
      n <- length(t)
      z <- matrix(c(rep(1,n),t),n,2)
      Zi <- superMatrix(list(Zi,z))
    }
    Z <- rbind(Z, Zi)
  }
  
  # ptm <- proc.time()
  # calculate inital estimate by marginal models
  
  alpha.hist <- alpha.old <- alpha <- alpha.ini
  beta.hist <- beta.old <- beta <- beta.ini
  
  sigma.e2.hist <- sigma.e2.old <- sigma.e2 <- sigma.e2.ini
  
  Sigma <- Sigma.ini
  
  
  gamma0.hist <- gamma0.old <- gamma0 <- gamma0.ini
  gamma.hist <- gamma.old <- gamma <- gamma.ini
  
  h0.hist <- h0.old <- h0 <- c(0.1,0.5,1.25,3)
  
  b <- GS()

  negPL.hist <- negPL <- negPL.fun()
  
  k <- 0
  rdiff <- 1
  scale <- 1/5
  maxiter.opt <- 10
  
  while (k < maxiter & rdiff > rtol){
    
    #ptm <- proc.time()
    b <- GS() 
    
    alpha.old <- alpha
    beta.old <- beta
    gamma0.old <- gamma0
    gamma.old <- gamma
    h0.old <- h0
    Sigma.old <- Sigma
    sigma.e2.old <- sigma.e2
    negPL.old <- negPL
    
    f <- f.new <- negPL.fun()
    
    # update h0
    
    
    iter <- 0
    
    
    for (i in 1:n.break) {
      data.sub.id <- data.T$id[t.break[i] < data.T$time & data.T$time <= t.break[i+1]]
      data.T.sub <- data.T[data.T$id %in% data.sub.id, ]
      data.Y.sub <- data.Y[data.Y$id %in% data.sub.id, ]
      
      data.X.sub <- data.Y.sub[!duplicated(data.Y.sub$id),]
      
      
      drv <- 0
      
      for (id in data.sub.id) {
        data.Y.row.sel <- data.Y.sub$id==id
        data <- data.Y.sub[data.Y.row.sel,]
        N <- nrow(data)
        X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
        
        data.T.row.sel <- data.T.sub$id==id
        time <- data.T.sub[data.T.row.sel,]$time
        delta <- data.T.sub[data.T.row.sel,]$delta
        
        t <- data[data$y.name==1,]$y.t
        n <- length(t)
        Zi <- matrix(c(rep(1,n),t),n,2)
        for (j in 2:J) {
          t <- data[data$y.name==j,]$y.t
          n <- length(t)
          z <- matrix(c(rep(1,n),t),n,2)
          Zi <- superMatrix(list(Zi,z))
        }
        eta <- as.numeric(X[1,]%*%gamma0+gamma%*%b[rownames(b)==as.character(id),])
        H <- -t(Zi)%*%Zi/sigma.e2 - h0[i]*time*exp(eta)*gamma%*%t(gamma) - solve(Sigma)
        drv <- drv - sum(diag(solve(H)%*%gamma%*%t(gamma)))*h0[i]*time*exp(eta)/2 - delta/h0[i] + time*exp(c(as.matrix(data.X.sub[data.X.sub$id==id,substr(colnames(data.X.sub),1,2)=="x."])%*%gamma0)+c(t(b[rownames(b)==as.character(id),])%*%gamma))
      }
      
      
      if (abs(drv) > 1e-6) {
        iter <- 0
        step <- 0.5*0.9^k
        h0.update <- h0[i] - step*drv
        while (h0.update < 0) {
          step <- scale*step
          iter <- iter + 1
          h0.update <- h0[i] - step*drv
        } 
        if (abs(step*drv) > 1e-5){
          f.new <- negPL.fun(h0.est = replace(h0, i, h0.update))
          while ( abs(step*drv) >1e-5 & f - f.new  < step*drv^2/4) {
            step <- scale*step
            iter <- iter+1
            h0.update <- h0[i] - step*drv
            f.new <- negPL.fun(h0.est = replace(h0, i, h0.update))
          }
        }
        
        h0[i] <- ifelse(abs(step*drv)<1e-5, h0[i], h0[i] - step*drv)
        f <- ifelse(abs(step*drv)<1e-5, f, f.new)
      }
    }
    
    print(paste("h0",h0,",","iter",iter))   
    
    
    # update gamma
    iter <- 0
    
    for (i in 1:J){
      ###gamma.old=gamma[2*i-1: 2*i]
      drv1 <- drv2 <- H11<- H22 <- 0
      for (id in idRange){
        data <- data.Y[data.Y$id==id,]
        N <- nrow(data)
        X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
        time <- data.T[data.T$id==id,]$time
        delta <- data.T[data.T$id==id,]$delta
        
        t <- data[data$y.name==1,]$y.t
        n <- length(t)
        Zi <- matrix(c(rep(1,n),t),n,2)
        for (j in 2:J) {
          t <- data[data$y.name==j,]$y.t
          n <- length(t)
          z <- matrix(c(rep(1,n),t),n,2)
          Zi <- superMatrix(list(Zi,z))
        }
        
        for (int in 1:n.break) {
          if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0[int] }
        }   
        
        eta <- as.numeric(X[1,]%*%gamma0+gamma%*%b[rownames(b)==as.character(id),])
        H <- -t(Zi)%*%Zi/sigma.e2 - h.base*time*exp(eta)*gamma%*%t(gamma) - solve(Sigma)
        ind1<-ind2 <- numeric(length(gamma))
        ind1[2*i-1] <- 1
        ind2[2*i] <- 1
        
        drv1i <- -delta*b[rownames(b)==as.character(id),(2*i-1)] + h.base*time*exp(eta)*b[rownames(b)==as.character(id),(2*i-1)]-sum(diag(solve(H)%*%(gamma%*%t(gamma)*c(b[rownames(b)==as.character(id),(2*i-1)])+gamma%*%t(ind1)+ind1%*%t(gamma))))*h.base*time*exp(eta)/2
        drv2i <- -delta*b[rownames(b)==as.character(id),(2*i)] + h.base*time*exp(eta)*b[rownames(b)==as.character(id),(2*i)]-sum(diag(solve(H)%*%(gamma%*%t(gamma)*c(b[rownames(b)==as.character(id),(2*i)])+gamma%*%t(ind2)+ind2%*%t(gamma))))*h.base*time*exp(eta)/2
        
        h11 <- h.base*time*exp(eta)*(b[rownames(b)==as.character(id),(2*i-1)])^2-sum(diag(solve(H)%*%(gamma%*%t(gamma)*c(b[rownames(b)==as.character(id),(2*i-1)])+gamma%*%t(ind1)+ind1%*%t(gamma))))*h.base*time*exp(eta)*b[rownames(b)==as.character(id),(2*i-1)]/2-
          h.base*time*exp(eta)*((b[rownames(b)==as.character(id),(2*i-1)]*solve(H)%*%gamma+ solve(H)%*%ind1)[2*i-1])
        h22 <- h.base*time*exp(eta)*(b[rownames(b)==as.character(id),(2*i)])^2-sum(diag(solve(H)%*%(gamma%*%t(gamma)*c(b[rownames(b)==as.character(id),(2*i)])+gamma%*%t(ind2)+ind2%*%t(gamma))))*h.base*time*exp(eta)*b[rownames(b)==as.character(id),(2*i)]/2-
          h.base*time*exp(eta)*((b[rownames(b)==as.character(id),(2*i)]*solve(H)%*%gamma+ solve(H)%*%ind2)[2*i])       
        
        
        drv1 <- as.numeric(drv1 + drv1i)
        drv2 <- as.numeric(drv2 + drv2i)
        H11 <- as.numeric(H11 + h11)
        H22 <- as.numeric(H22 + h22)
      }
      
      drv= c(drv1, drv2)
      ### drv <- drv1 + sign(gamma[i])*lambda*length(idRange)
      hg= min(max(H11, H22, 0.01), 100) ###-max(-1*H11, -1*H22, 1)
      
      if ( (sqrt(sum((drv-hg*gamma[(2*i-1):(2*i)])^2))<=sqrt(2)*lambda) ||(sum((gamma[(2*i-1):(2*i)])^2)<=2*lambda^2) ){
        dg= gamma[(2*i-1):(2*i)]
        ### f.new <- negPL.fun(gamma.est = replace(gamma,c(2*i-1, 2*i),c(0,0)))
        gamma[(2*i-1)]= 0   ### ifelse(f.new <= f, 0, gamma[(2*i-1)])
        gamma[(2*i)]= 0     ### ifelse(f.new<= f,  0, gamma[(2*i)])
        ### f <- min(f, f.new)
      }else{
        dg= 1/hg*(drv- sqrt(2)*lambda*(drv-hg*gamma[(2*i-1):(2*i)])/ sqrt(sum((drv-hg*gamma[(2*i-1):(2*i)])^2)))
        
        iter <- 0
        step <- 1*(0.9)^k
        
        if (step*min(abs(dg))>1e-9){
          
          ###     while (step*min(abs(dg))> 10){ ###abs(gamma[i])/10
          ###        step <- scale*step
          ###        iter <- iter+1
          ###      }
          f.new <- negPL.fun(gamma.est = replace(gamma,c(2*i-1, 2*i),gamma[(2*i-1):(2*i)]-step*dg))
          delta.drv= sum(dg*drv)- sqrt(2)*lambda*sqrt(sum((gamma[(2*i-1):(2*i)]-hg )^2))+ sqrt(2)*lambda* sqrt(sum((gamma[(2*i-1):(2*i)])^2))
          
          while ( (step*min(abs(dg)) >1e-5) & (f - f.new < step*delta.drv/500)){
            step <- scale*step
            iter <- iter+1
            ###gamma[i]= sign(gamma[i]-step*drv)*max(abs(gamma[i]-step*drv)-lambda, 0)  
            f.new <- negPL.fun(gamma.est = replace(gamma,c(2*i-1, 2*i),gamma[(2*i-1):(2*i)]-step*dg))
          }
          
          gamma[(2*i-1)]= ifelse(f.new <= f, (gamma[(2*i-1):(2*i)]-step*dg)[1], gamma[(2*i-1)])
          gamma[(2*i)]= ifelse(f.new<= f,  (gamma[(2*i-1):(2*i)]-step*dg)[2], gamma[(2*i)])
          f <- min(f, f.new)
        }
      }
    }
    
    
    print(paste("gamma",gamma,",","iter",iter))    
    
    
    # update alpha
    iter <- 0
    mu <- sapply(1:nrow(data.Y), function(n){
      y.name.ind <- data.Y$y.name[n]
      alpha.sel <- (J1*(y.name.ind-1)+1):(J1*y.name.ind)
      beta.sel <- 2*y.name.ind-1
      b.row.sel <- rownames(b)==as.character(data.Y$id[n])
      as.matrix(data.Y[n,substr(colnames(data.Y),1,2)=="x."])%*%alpha[alpha.sel] + beta[beta.sel] + beta[beta.sel+1]*data.Y$y.t[n] + b[b.row.sel,beta.sel] + b[b.row.sel,beta.sel+1]*data.Y$y.t[n]
    })
    if (sum(abs(data.Y$y-mu)) > 1e-6){
      for (idx1 in 1:(length(alpha))){
        x.idx1 <- ifelse(idx1%%J1==0, J1, idx1%%J1)
        y.idx1 <- (idx1-x.idx1)/J1+1
        
        iter <- 0
        #mu <- sapply(1:nrow(data.Y), function(n) as.matrix(data.Y$x[n])%*%alpha[(J1*(data.Y$y.name[n]-1)+1):(J1*data.Y$y.name[n])] + beta[2*data.Y$y.name[n]-1] + beta[2*data.Y$y.name[n]]*data.Y$y.t[n] + b[rownames(b)==as.character(data.Y$id[n]),2*data.Y$y.name[n]-1] + b[rownames(b)==as.character(data.Y$id[n]),2*data.Y$y.name[n]]*data.Y$y.t[n])
        #drv <- -sum(((data.Y$y-mu)*as.matrix(data.Y$x)[,x.idx1]/(sigma.e2))[data.Y$y.name==y.idx1])/sigma.e2
        #drvdrv <- sum(as.matrix(data.Y$x)[data.Y$y.name==y.idx1,x.idx1]^2)/sigma.e2 
        
        drvdrv <- sum(data.Y[,substr(colnames(data.Y),1,2)=="x."][data.Y$y.name==y.idx1,x.idx1]^2)/sigma.e2  
        
        alpha.est <- alpha
        drv <- 1
        while (iter < 10 & sum(drv^2) > 1e-5){
          mu <- sapply(1:nrow(data.Y), function(n){
            y.name.ind <- data.Y$y.name[n]
            alpha.sel <- (J1*(y.name.ind-1)+1):(J1*y.name.ind)
            beta.sel <- 2*y.name.ind-1
            b.row.sel <- rownames(b)==as.character(data.Y$id[n])
            as.matrix(data.Y[n,substr(colnames(data.Y),1,2)=="x."])%*%alpha.est[alpha.sel] + beta[beta.sel] + beta[beta.sel+1]*data.Y$y.t[n] + b[b.row.sel,beta.sel] + b[b.row.sel,beta.sel+1]*data.Y$y.t[n]
          })
          
          drv <- -sum(((data.Y$y-mu)*data.Y[,substr(colnames(data.Y),1,2)=="x."][,x.idx1])[data.Y$y.name==y.idx1])/sigma.e2
          
          alpha.est.old <- alpha.est      
          alpha.est[idx1] <- alpha.est.old[idx1] - drv/drvdrv
          
          iter <- iter+1
        }
        alpha <- alpha.est
      }
    }
    
    print(paste("alpha",alpha,",","iter",iter))
    
    
    # update beta
    
    
    for (i in 1:length(beta)){  
      iter <- 0
      mu <- sapply(1:nrow(data.Y), function(n){
        y.name.ind <- data.Y$y.name[n]
        alpha.sel <- (J1*(y.name.ind-1)+1):(J1*y.name.ind)
        beta.sel <- 2*y.name.ind-1
        b.row.sel <- rownames(b)==as.character(data.Y$id[n])
        as.matrix(data.Y[n,substr(colnames(data.Y),1,2)=="x."])%*%alpha[alpha.sel] + beta[beta.sel] + beta[beta.sel+1]*data.Y$y.t[n] + b[b.row.sel,beta.sel] + b[b.row.sel,beta.sel+1]*data.Y$y.t[n]
      })
      
      drv <- (-t(data.Y$y-mu)%*%Z/sigma.e2)[i]
      
      step<- 0.1*(0.9)^k
      f.new <- negPL.fun(beta.est = replace(beta, i, beta[i]-step*drv))
      while ( (iter < 10) & (step*abs(drv) >1e-5) & (f - f.new < step*drv^2/10)){
        step <- scale*step
        iter <- iter+1  
        f.new <- negPL.fun(beta.est = replace(beta, i, beta[i]-step*drv))
      }
      beta[i]= ifelse(f.new<= f, beta[i]-step*drv, beta[i])
      f <- min(f, f.new)
      iter <- iter+1
    }
    print(paste("beta",beta,",","iter",iter))
    
    # update sigma.e2
    
    iter <- 0
    
    drv1 <- sum(sapply(idRange, function(id){
      data.Y.row.sel <- data.Y$id==id
      data <- data.Y[data.Y.row.sel,]
      N <- nrow(data)
      X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
      data.T.row.sel <- data.T$id==id
      time <- data.T[data.T.row.sel,]$time
      delta <- data.T[data.T.row.sel,]$delta
      
      t <- data[data$y.name==1,]$y.t
      n <- length(t)
      Zi <- matrix(c(rep(1,n),t),n,2)
      for (j in 2:J) {
        t <- data[data$y.name==j,]$y.t
        n <- length(t)
        z <- matrix(c(rep(1,n),t),n,2)
        Zi <- superMatrix(list(Zi,z))
      }
      
      for (int in 1:n.break){
        if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0[int] }
      } 
      
      eta <- as.numeric(X[1,]%*%gamma0+gamma%*%b[rownames(b)==as.character(id),])
      H <- -t(Zi)%*%Zi/sigma.e2 - h.base*time*exp(eta)*gamma%*%t(gamma) - solve(Sigma)
      drv1i <- sum(diag(solve(H)%*%t(Zi)%*%Zi))/(sigma.e2^2)/2
      as.numeric(drv1i)
    }))
    
    mu <- sapply(1:nrow(data.Y), function(n){
      y.name.ind <- data.Y$y.name[n]
      alpha.sel <- (J1*(y.name.ind-1)+1):(J1*y.name.ind)
      beta.sel <- 2*y.name.ind-1
      b.row.sel <- rownames(b)==as.character(data.Y$id[n])
      as.matrix(data.Y[n, substr(colnames(data.Y),1,2)=="x."])%*%alpha[alpha.sel] + beta[beta.sel] + beta[beta.sel+1]*data.Y$y.t[n] + b[b.row.sel,beta.sel] + b[b.row.sel,beta.sel+1]*data.Y$y.t[n]
    } )
    drv2 <- nrow(data.Y)*1/sigma.e2/2 - c(t(data.Y$y-mu)%*%(data.Y$y-mu)/(sigma.e2^2*2))
    drv <- drv1 + drv2
    
    
    if (f - negPL.fun(sigma.e2.est = sigma.e2-sign(drv)*1e-5) > (1e-5)^2/4){
      
      iter <- 0
      step <- 0.01*0.9^k
      
      sigma.e2.update <- sigma.e2 - step*drv
      while (sigma.e2.update < 0){
        step <- scale*step
        iter <- iter + 1
        sigma.e2.update <- sigma.e2 - step*drv
      } 
      if (abs(step*drv) >1e-5){
        f.new <- negPL.fun(sigma.e2.est=sigma.e2.update)
        while ( abs(step*drv) >1e-5 & f - f.new < step*drv^2/4){
          step <- scale*step
          iter <- iter+1
          sigma.e2.update <- sigma.e2 - step*drv
          f.new <- negPL.fun(sigma.e2.est = sigma.e2.update) 
        }
      }
      
      sigma.e2 <- ifelse(abs(step*drv)<1e-5, sigma.e2, sigma.e2.update)
      f <- ifelse(abs(step*drv)<1e-5, f, f.new)
    }
    
    
    print(paste("sigma.e2",sigma.e2,",","iter",iter))
    
    
    
    # update gamma0
    iter <- 0
    
    drv <- rep(0,length(gamma0))
    drv.list <- lapply(idRange, function(id){
      data <- data.Y[data.Y$id==id,]
      N <- nrow(data)
      X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
      time <- data.T[data.T$id==id,]$time
      delta <- data.T[data.T$id==id,]$delta
      
      t <- data[data$y.name==1,]$y.t
      n <- length(t)
      Zi <- matrix(c(rep(1,n),t),n,2)
      for (j in 2:J) {
        t <- data[data$y.name==j,]$y.t
        n <- length(t)
        z <- matrix(c(rep(1,n),t),n,2)
        Zi <- superMatrix(list(Zi,z))
      }
      
      for (int in 1:n.break) {
        if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0[int] }
      } 
      
      eta <- as.numeric(X[1,]%*%gamma0+gamma%*%b[rownames(b)==as.character(id),])
      H <- -t(Zi)%*%Zi/sigma.e2 - h.base*time*exp(eta)*gamma%*%t(gamma) - solve(Sigma)
      drv1i <- -sum(diag(solve(H)%*%gamma%*%t(gamma)))*h.base*time*exp(eta)*X[1,]/2
      
      drv2i <- -delta%*%as.matrix(data.X[data.X$id==id, substr(colnames(data.X),1,2)=="x."]) + as.numeric(h.base*time*exp(as.matrix(data.X[data.X$id==id, substr(colnames(data.X),1,2)=="x."])%*%gamma0+b[rownames(b)==as.character(id),]%*%gamma))*as.matrix(data.X[data.X$id==id, substr(colnames(data.X),1,2)=="x."])
      
      drv1i+drv2i
    })
    drv <- apply(matrix(unlist(drv.list), nrow=length(gamma0), byrow=TRUE), 2, sum)
    
    for (i in 1:length(gamma0)){
      if (abs(drv[i]) > 1e-5) {
        
        if (f - negPL.fun(gamma0.est = replace(gamma0,i,gamma0[i]-sign(drv[i])*1e-5)) > (1e-5)^2/4){
          
          step <- 0.01*0.9^k
          iter <- 0
          
          while ( abs(step*drv[i])>1e-5 & f - f.new  < step*drv[i]^2/4 ){
            step <- scale*step
            iter <- iter+1
            f.new <- negPL.fun(gamma0.est = replace(gamma0,i,gamma0[i]-step*drv[i]))
          }
          
          gamma0[i] <- ifelse(abs(step*drv[i])<1e-5, gamma0[i], gamma0[i] - step*drv[i])
          f <- ifelse(abs(step*drv[i])<1e-5, f, f.new)
        }
      }
    }
    
    print(paste("gamma0",gamma0,",","iter",iter))
    
    
    # update Sigma
    iter <- 0
    for (i in 1:(2*J)){
      for (j in i:(2*J)){  
        if (Sigma[i,j]!=0){
          Sigma.inv <- solve(Sigma)
          Sigma.drv <- matrix(numeric(ncol(Sigma)*ncol(Sigma)), ncol=ncol(Sigma))
          Sigma.drv[i,j] <- Sigma.drv[j,i] <- 1
          
          drv <- sum(sapply(idRange, function(id){
            data.Y.row.sel <- data.Y$id==id 
            data <- data.Y[data.Y.row.sel,]
            N <- nrow(data)
            X <- as.matrix(data[,substr(colnames(data),1,2)=="x."])
            data.T.row.sel <- data.T$id==id
            time <- data.T[data.T.row.sel,]$time
            delta <- data.T[data.T.row.sel,]$delta
            
            t <- data[data$y.name==1,]$y.t
            n <- length(t)
            Zi <- matrix(c(rep(1,n),t),n,2)
            for (jj in 2:J) {
              t <- data[data$y.name==jj,]$y.t
              n <- length(t)
              z <- matrix(c(rep(1,n),t),n,2)
              Zi <- superMatrix(list(Zi,z))
            }
            
            for (int in 1:n.break) {
              if (t.break[int] < time & time <= t.break[int+1]) { h.base <- h0[int] }
            } 
            
            eta <- as.numeric(X[1,]%*%gamma0+gamma%*%b[rownames(b)==as.character(id),])
            H <- -t(Zi)%*%Zi/sigma.e2 - h.base*time*exp(eta)*gamma%*%t(gamma) - Sigma.inv
            
            drv1i <- sum(diag(Sigma.inv%*%Sigma.drv))/2 - t(b[rownames(b)==as.character(id),])%*%Sigma.inv%*%Sigma.drv%*%Sigma.inv%*%b[rownames(b)==as.character(id),]/2
            drv2i <- sum(diag(solve(H)%*%Sigma.inv%*%Sigma.drv%*%Sigma.inv))/2
            as.numeric(drv1i + drv2i)
          }))
          
          ###drv <- drv + (i!=j)*lambda2*sign(Sigma[i,j])*length(idRange)*abs(Sigma[i,j])
          if (abs(drv) > 1e-5){
            ###if (f - negPL.fun(Sigma.est = replace(Sigma,c((j-1)*ncol(b)+i,(i-1)*ncol(b)+j),Sigma[i,j]-sign(drv)*1e-5)) > (1e-5)^2/4){
            iter <- 0
            step <- 0.01*0.9^k
            if (abs(step*drv) > 1e-5 ){
              Sigma.update <- replace(Sigma,c((j-1)*ncol(b)+i,(i-1)*ncol(b)+j),Sigma[i,j]-step*drv)
              while (det(Sigma.update) < 0){
                step <- scale*step
                iter <- iter + 1
                Sigma.update <- replace(Sigma,c((j-1)*ncol(b)+i,(i-1)*ncol(b)+j),Sigma[i,j]-step*drv)
              } 
              
              f.new <- negPL.fun(Sigma.est = Sigma.update)
              while ( abs(step*drv) > 1e-5 & f - f.new < step*drv^2/4){
                step <- scale*step
                iter <- iter+1
                Sigma.update <- replace(Sigma,c((j-1)*ncol(b)+i,(i-1)*ncol(b)+j),Sigma[i,j]-step*drv)
                f.new <- negPL.fun(Sigma.est = Sigma.update)
              }
            }
            Sigma[i,j] <- Sigma[j,i] <- ifelse(f.new>f, Sigma[i,j], Sigma[i,j]-step*drv)
            f= min(f, f.new)
            if (i!=j){
              Sigma.lasso= sign(Sigma[i,j])*max(abs(Sigma[i,j])-lambda2, 0)
              Sigma.update1 <- replace(Sigma,c((j-1)*ncol(b)+i,(i-1)*ncol(b)+j),Sigma.lasso)
              f.new <- negPL.fun(Sigma.est = Sigma.update1)
              Sigma[i,j] <- Sigma[j,i] <- ifelse(f.new< f, Sigma.lasso, Sigma[i,j])
              f <- min(f, f.new)###ifelse(abs(step*drv)<1e-5, f, f.new)
            }#####sign(Sigma[i,j]-step*drv)*max(abs(Sigma[i,j]-step*drv)-lambda2, 0)
            ###}
          }        
        }
      }
    }
    print(Sigma)    
    
    
    gamma.hist <- rbind(gamma.hist, gamma)
    gamma0.hist <- rbind(gamma0.hist, gamma0)
    alpha.hist <- rbind(alpha.hist, alpha)
    beta.hist <- rbind(beta.hist, beta)
    sigma.e2.hist <- rbind(sigma.e2.hist, sigma.e2)
    h0.hist <- rbind(h0.hist, h0)
    
    negPL <- f
    negPL.hist <- c(negPL.hist, negPL)
    rdiff <- c(abs(negPL - negPL.old)/(abs(negPL.old)+0.000001))
    k <- k+1
    
    print(paste("rdiff",rdiff,",","iter",k-1))
    
    print(list(alpha=alpha.hist,beta=beta.hist,gamma=gamma.hist,gamma0=gamma0.hist,
               Sigma=Sigma,sigma.e2=sigma.e2.hist,h0=h0.hist,negPL=negPL.hist, rdiff=rdiff,iter=k))
    
    print(list(alpha=alpha,beta=beta,gamma=gamma,gamma0=gamma0,
               Sigma=Sigma,sigma.e2=sigma.e2,h0=h0,negPL=negPL, rdiff=rdiff,iter=k))
    
  }
  out1 <- list(alpha=alpha, beta=beta, gamma=gamma, gamma0=gamma0, sigma.e2=sigma.e2, Sigma=Sigma, negPL=negPL)
  
}



