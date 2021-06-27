# R Code for 'A Bayesian phase I/II platform design for co-developing drug combination therapies for multiple indications'
# Contact the first author: Rongji Mu, rjmu@sjtu.edu.cn.

library(dlm);library(msm);library(mvtnorm)

effdose <- function(prior.pt.mat, prior.pe.mat, ngroup, vd, type){
  drange = range(vd)
  # coefficients of linear transformation
  coef.t = c((drange[2]+ drange[1])/2/(drange[1]-drange[2]), 1/(drange[2]-drange[1])) 
  nd = NULL
  u.mat = array(list(),dim=c(1,1,ngroup))
  vbeta = vgamma = matrix(0,ngroup,2)
  for(k in 1:ngroup){
    nd[k] = length(prior.pt.mat[[k]])
    for(j in 1:nd[k]){
      u.mat[[k]][j] = round(coef.t[1]+coef.t[2]*vd[[k]][j],2)
    }
    vbeta[k,] = lm(qnorm(prior.pt.mat[[k]])~u.mat[[k]])$coefficients
    if(type=='binary'){
      vgamma[k,] = lm(qnorm(prior.pe.mat[[k]])~u.mat[[k]])$coefficients
    }else{
      vgamma[k,] = lm(prior.pe.mat[[k]]~u.mat[[k]])$coefficients
    }
  }
  list(umat=u.mat, vmat=u.mat, nd=nd, vbeta=apply(vbeta,2,mean), vgamma = apply(vgamma,2,mean))
}

outcome <- function(dose.ind,  cohortsize, pt, vpara, u, v, type) {
  # 
  gamma0.true = vpara[1]; 
  gamma1.true = vpara[2]; 
  alpha.true = vpara[3] # efficacy model (4)
  rho.true = vpara[4] # used in the joint bivariate model 
  
  Y_T <- rbinom(cohortsize, 1, pt[dose.ind])
  n.y = 10000
  Y_T.temp<-rbinom(n.y,1,pt[dose.ind])
  
  ztv <- ifelse(pt[dose.ind]==0.5, 0.5/dnorm(pt[dose.ind]) ,qnorm(pt[dose.ind])/(2*pt[dose.ind]-1))
  Z_T <- ifelse(Y_T.temp==0, -ztv, ztv)
  Z_T_mean <- mean(Z_T)
  Z_E_mean <- gamma0.true+gamma1.true*ifelse(v[dose.ind]<= alpha.true, 
                                             v[dose.ind], alpha.true)
  
  Z <- rmvnorm(cohortsize, mean = c(Z_T_mean, Z_E_mean),
               sigma = matrix(c(1,rho.true, rho.true,1),nrow=2)) 
  #print(Z)
  if(type=='binary'){
    Y_E <- ifelse(Z[,2]>0, 1, 0)
  }else{
    Y_E <- Z[,2]
  }
  return(cbind(Y_T, Y_E, rep(u[dose.ind], cohortsize), rep(v[dose.ind], cohortsize)))
}

tox.prob <- function(vp,u) {
  beta0 = vp[1]
  beta1 = vp[2]
  Z_T_mean <- beta0+beta1*u
  return(pnorm(Z_T_mean))
}

eff.prob <- function(vp,v,type) {
  c_Y = 2
  gamma0 = vp[3]
  gamma1 = vp[4]
  alpha = vp[5]
  Z_E_mean <- gamma0+gamma1*ifelse(v<=alpha, v, alpha)
  if(type=='binary'){
    pe = pnorm(Z_E_mean)
  }else{
    Z_E.sample <- rnorm(10000, Z_E_mean, sd=1)
    pe = mean(Z_E.sample>c_Y)
  }
  return(pe)
}


norm2d <- function(v,rho, n) rmvnorm(n=n,mean=v,sigma=matrix(c(1,rho,rho,1),nrow=2))

get.utility <- function(vp,u,v,type) {
  y_E.cut = 2
  beta0 = vp[1]
  beta1 = vp[2]
  gamma0 = vp[3]
  gamma1 = vp[4]
  alpha = vp[5]
  rho = vp[6]
  # 
  Uti <- rbind(c(0, 50),
               c(25, 100))
  # V1
  # Uti <- rbind(c(0, 55),
  #              c(20, 100))
  n.y <- 100000
  prob <- matrix(0, nrow = 2, ncol = 2)
  mu.T <- beta0+beta1*u
  mu.E <- gamma0+gamma1*ifelse(v<=alpha, v, alpha)
  sample2d <- norm2d(v=c(mu.T, mu.E), rho=rho, n=n.y)
  
  c_Y = ifelse(type=='binary', 0, y_E.cut)
  
  prob[1,1] <- sum((sample2d[,1]>0) & (sample2d[,2]<=c_Y))/n.y
  prob[1,2] <- sum((sample2d[,1]>0) & (sample2d[,2]>c_Y))/n.y
  prob[2,1] <- sum((sample2d[,1]<=0) & (sample2d[,2]<=c_Y))/n.y
  prob[2,2] <- sum((sample2d[,1]<=0) & (sample2d[,2]>c_Y))/n.y
  uti <- sum(Uti*prob)
  return(uti)
}		

summary.mcmc <- function(matrix.uti,n.dose,u,v,type) {
  
  mean.matrix.uti <- apply(matrix.uti,2,mean)
  uti <- rep(0,n.dose)
  tox.mcmc <- rep(0,n.dose)
  eff.mcmc <- rep(0,n.dose)
  for (i in 1:n.dose) {
    tox.prob.mcmc <- apply(matrix.uti,1,tox.prob,u=u[i])
    eff.prob.mcmc <- apply(matrix.uti,1,eff.prob,v=v[i],type)
    tox.mcmc[i] <- sum(tox.prob.mcmc<phi.T)/nrow(matrix.uti)
    eff.mcmc[i] <- sum(eff.prob.mcmc>phi.E)/nrow(matrix.uti)
    uti[i] <- get.utility(mean.matrix.uti,u=u[i],v=v[i],type)
  }
  list(tox=tox.mcmc, eff=eff.mcmc, uti=uti)
}

mcmc.arms <- function(dat, vbeta=vbeta, vgamma=vgamma, ngroup, type) {
  N.post <- 10000
  N.burnin <- 5000
  ngroup = ngroup   # ngroup: number of indications or trials
  # set initial values
  # parameter appears in model (2)
  rho <- 0 
  
  beta0.hat <- vbeta[1] # 0.01, 0.3 
  beta1.hat <- vbeta[2]
  gamma0.hat <- vgamma[1]
  gamma1.hat <- vgamma[2]
  alpha.hat <- 0
  
  # parameters appear in toxicity model (3)
  beta0 <- rep(beta0.hat, ngroup)
  beta1 <- rep(beta1.hat, ngroup)
  beta1_share <- log(beta1.hat)  # shared among ngroup groups
  mu.beta0 <- beta0.hat
  tau2.beta0 <- (4*mu.beta0)^2
  mu.beta1 <- log(beta1.hat)
  
  
  tau2.beta1 <- ((log(3)-mu.beta1)/qnorm(0.9))^2 
  
  a.beta1 <- (2*mu.beta1)^2
  sigma2.beta1 <- 0.5*a.beta1 # \sigma_{\beta_1}^2 in manuscript
  
  
  # parameters appear in efficacy model (4)
  gamma0 <- rep(gamma0.hat, ngroup)
  gamma1 <- rep(gamma1.hat, ngroup)
  alpha <- rep(alpha.hat, ngroup)
  gamma1_share <- gamma1.hat  # shared among ngroup groups
  mu.gamma0 <- gamma0.hat
  tau2.gamma0 <- (4*mu.gamma0)^2
  
  mu.gamma1 <- log(gamma1.hat) 
  
  if(type=='binary'){
    tau2.gamma1 <- ((log(3)-mu.gamma1)/qnorm(0.9))^2
    a.gamma1 <- 4*mu.gamma1^2  
    sigma2.gamma1 <- 0.5*a.gamma1 # \sigma_{\gamma_1}^2 in manuscript
  }else{
    tau2.gamma1 <- 4*mu.gamma1^2
    a.gamma1 <- 4*mu.gamma1^2 
    sigma2.gamma1 <-  0.5*a.gamma1# \sigma_{\gamma_1}^2 in manuscript
  }
  
  #sigma2.gamma1 <- 0.5*a.gamma1 # \sigma_{\gamma_1}^2 in manuscript
  #tau2.gamma1 <- 0.5*a.gamma1
  
  
  n <- NULL
  for(k in 1:ngroup){
    n[k] = nrow(dat[[k]])
  }
  
  # N.post:
  rho_t <- rep(0, N.post) 
  beta0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  beta1_share_t <- rep(0, N.post)  
  sigma2.beta1_t <- rep(0, N.post) 
  
  gamma0_t <- matrix(0, nrow=N.post, ncol = ngroup)
  gamma1_t <- matrix(0, nrow=N.post, ncol = ngroup)
  alpha_t <- matrix(0, nrow=N.post, ncol = ngroup) 
  gamma1_share_t <- rep(0, N.post)  
  sigma2.gamma1_t <- rep(0, N.post) 
  
  
  # latent variables 
  Z_T = array(list(),c(1,1,ngroup))
  Z_T_t <- array(list(), dim = c(1, 1, ngroup))
  
  
  Z_E = array(list(),c(1,1,ngroup))
  Z_E_t <- array(list(), dim = c(1, 1, ngroup)) # Y_E for normal
  
  for(k in 1:ngroup){
    Z_T[[k]] = ifelse(dat[[k]][,1]==0, -0.5, 0.5) 
    if(type=='binary'){
      Z_E[[k]] = ifelse(dat[[k]][,2]==0, -0.5, 0.5)
    }else{
      Z_E[[k]] = dat[[k]][,2]  # Y_E
    }
  }
  
  for(ite in 1:N.post){
    #print(ite)
    # generate the latent variables
    xi_cut <- c(-Inf, 0, Inf)
    zeta_cut <- c(-Inf, 0, Inf)
    log.likelihood <- function(rho, beta0, beta1, gamma0, gamma1, alpha){
      Z_T.C <- Z_E.C <- NULL
      for(j in 1:ngroup){
        mean.Z_T <- beta0[j] + beta1[j]*dat[[j]][,3]
        mean.Z_E <- gamma0[j] + gamma1[j]*ifelse(dat[[j]][,4]<=alpha[j], dat[[j]][,4], alpha[j])
        Z_T.C <- c(Z_T.C, Z_T[[j]]-mean.Z_T)
        Z_E.C <- c(Z_E.C, Z_E[[j]]-mean.Z_E)
      }
      return(-length(Z_T.C)/2*log(1-rho^2)-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C+Z_T.C^2))
    }
    
    # rho
    logden <- function(x) log.likelihood(x, beta0, beta1, gamma0, gamma1, alpha)
    rho <- arms(rho, logden, function(x) ((x>-0.8)*(x<0.8)), 1)
    rho_t[ite] <- rho 
    var.z <- 1-rho^2 # rho: sigma_{12}^2
    
    for(k in 1:ngroup){
      
      mean.Z.E <- gamma0[k] + gamma1[k]*ifelse(dat[[k]][,4]<=alpha[k], dat[[k]][,4], alpha[k]) - rho*(Z_T[[k]]-beta0[k] - beta1[k]*dat[[k]][,3]) 
      if(type=='binary'){
        Z_E[[k]] <- rtnorm(n[k], mean=mean.Z.E, sd = sqrt(var.z), lower = xi_cut[dat[[k]][,2]+1], upper = xi_cut[dat[[k]][,2]+2])
      }	  
      # 
      mean.Z.T <- beta0[k] + beta1[k]*dat[[k]][,3] - rho*(Z_E[[k]]-gamma0[k] - gamma1[k]*ifelse(dat[[k]][,4]<=alpha[k], dat[[k]][,4], alpha[k]))
      Z_T[[k]] <- rtnorm(n[k], mean=mean.Z.T, sd = sqrt(var.z), lower = xi_cut[dat[[k]][,1]+1], upper = xi_cut[dat[[k]][,1]+2])
      
      
      
      # Z_T_t[ite,,k] <- Z_T[[k]] # no used
      # Z_E_t[ite,,k] <- Z_E[[k]]
      
      # beta0k
      mean.Z_E <- gamma0[k] + gamma1[k]*ifelse(dat[[k]][,4]<=alpha[k], dat[[k]][,4], alpha[k])    
      mu.n.beta0 <- mean(Z_T[[k]]-beta1[k]*dat[[k]][,3]-rho*(Z_E[[k]]-mean.Z_E))
      sigma2.n.beta0 <- (1-rho^2)/n[k]
      mean.beta0 <- (tau2.beta0*mu.n.beta0+sigma2.n.beta0*mu.beta0)/(tau2.beta0+sigma2.n.beta0)
      var.beta0 <- tau2.beta0*sigma2.n.beta0/(tau2.beta0+sigma2.n.beta0)
      beta0[k]  <- rnorm(1, mean = mean.beta0, sd=sqrt(var.beta0))
      beta0_t[ite,k] <- beta0[k]
      
      # beta1k
      log.likelihood.beta1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.beta1, beta1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3]
        mean.Z_E <- gamma0 + gamma1*ifelse(dat[,4]<=alpha, dat[,4], alpha)
        Z_T.C <- Z_T-mean.Z_T
        Z_E.C <- Z_E-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_T.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.beta1*(log(beta1)-beta1_share)^2)
      }
      logden <- function(x) log.likelihood.beta1(beta0[k],x,gamma0[k], gamma1[k], alpha[k], rho, sigma2.beta1, beta1_share, 
                                                 dat[[k]], Z_T[[k]], Z_E[[k]])
      beta1[k]  <- arms(beta1[k], logden, function(x) ((x>0)*(x<3)), 1)
      beta1_t[ite,k] <- beta1[k]
      
      # gamma0k
      mean.Z_T <- beta0[k] + beta1[k]*dat[[k]][,3]
      vstar <- ifelse(dat[[k]][,4]<=alpha[k], dat[[k]][,4], alpha[k])
      mu.n.gamma0 <- mean(Z_E[[k]]-gamma1[k]*vstar-rho*(Z_T[[k]]-mean.Z_T))
      sigma2.n.gamma0 <-(1-rho^2)/n[k]
      
      mean.gamma0 <- (tau2.gamma0*mu.n.gamma0+sigma2.n.gamma0*mu.gamma0)/(tau2.gamma0+sigma2.n.gamma0)
      var.n.gamma0 <-tau2.gamma0*sigma2.n.gamma0/(tau2.gamma0+sigma2.n.gamma0)
      var.gamma0 <- tau2.gamma0*sigma2.n.gamma0 /(tau2.gamma0+sigma2.n.gamma0)
      gamma0[k]  <- rnorm(1, mean = mean.gamma0, sd=sqrt(var.gamma0))
      gamma0_t[ite,k] <- gamma0[k]
      
      # gamma1k
      log.likelihood.gamma1 <- function(beta0, beta1, gamma0, gamma1, alpha, rho, sigma2.gamma1, gamma1_share, dat, Z_T, Z_E){
        mean.Z_T <- beta0 + beta1*dat[,3]
        mean.Z_E <- gamma0 + gamma1*ifelse(dat[,4]<=alpha, dat[,4], alpha)
        Z_T.C <- Z_T-mean.Z_T
        Z_E.C <- Z_E-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C)-1/2/sigma2.gamma1*(log(gamma1)-gamma1_share)^2)
      }
      logden <- function(x) log.likelihood.gamma1(beta0[k],beta1[k],gamma0[k], x, alpha[k], rho, sigma2.gamma1, gamma1_share, 
                                                  dat[[k]], Z_T[[k]], Z_E[[k]])
      gamma1.Star <- ifelse(type=='binary', 4, 6) 
      gamma1[k]  <- arms(gamma1[k], logden, function(x) ((x>0)*(x<gamma1.Star)), 1)
      gamma1_t[ite,k] <- gamma1[k]
      
      # alphak
      log.likelihood.alpha <- function(alphak, rho, beta0k, beta1k, gamma0k, gamma1k, Z_Tk, Z_Ek, datk) {
        u = datk[,3]
        v = datk[,4]
        mean.Z_T <- beta0k + beta1k*u
        mean.Z_E <- gamma0k + gamma1k*ifelse(v<=alphak, v, alphak)
        Z_T.C <- Z_Tk-mean.Z_T
        Z_E.C <- Z_Ek-mean.Z_E
        return(-1/2/(1-rho^2)*sum(Z_E.C^2-2*rho*Z_T.C*Z_E.C))
      }
      logden <- function(x) log.likelihood.alpha(x, rho, beta0[k], beta1[k], gamma0[k], gamma1[k], Z_T[[k]], Z_E[[k]], dat[[k]])
      alpha[k] <- arms(alpha[k],logden, function(x) ((x>-1)*(x<1)), 1)
      alpha_t[ite,k] <- alpha[k]
    }
    # beta1_share
    mu.n.beta1_share <- mean(log(beta1))
    sigma2.n.beta1_share <- sigma2.beta1/ngroup  
    
    mean.beta1_share <- (tau2.beta1*mu.n.beta1_share+sigma2.n.beta1_share*mu.beta1)/(tau2.beta1+sigma2.n.beta1_share)
    var.beta1_share <- (tau2.beta1*sigma2.n.beta1_share)/(tau2.beta1+sigma2.n.beta1_share)
    
    beta1_share <- rnorm(1, mean=mean.beta1_share, sd = sqrt(var.beta1_share))
    beta1_share_t[ite] <- beta1_share
    
    # sigma2.beta1
    logden <- function(x) -ngroup/2*log(x)- sum((log(beta1)-beta1_share)^2)/(2*x)
    sigma2.beta1 <- arms(sigma2.beta1, logden, function(x) (x>0)*(x<a.beta1), 1)
    sigma2.beta1_t[ite] <- sigma2.beta1
    
    # gamma1_share
    mu.n.gamma1_share <- mean(log(gamma1))
    sigma2.n.gamma1_share <- sigma2.gamma1/ngroup
    
    mean.gamma1_share <- (tau2.gamma1*mu.n.gamma1_share + sigma2.n.gamma1_share*mu.gamma1)/(tau2.gamma1+sigma2.n.gamma1_share)
    var.gamma1_share <- (tau2.gamma1*sigma2.n.gamma1_share)/(tau2.gamma1+sigma2.n.gamma1_share)
    gamma1_share <- rnorm(1, mean = mean.gamma1_share, sd=sqrt(var.gamma1_share))
    gamma1_share_t[ite] <- gamma1_share
    
    # sigma2.gamma1
    logden <- function(x) -ngroup/2*log(x) - sum((log(gamma1)-gamma1_share)^2)/(2*x)
    sigma2.gamma1 <- arms(sigma2.gamma1, logden, function(x) (x>0)*(x<a.gamma1), 1)
    sigma2.gamma1_t[ite] <- sigma2.gamma1
  }
  ind <- seq((N.burnin+1),N.post)
  # matrix.post <- cbind()[ind, ]
  matrix.uti <- array(list(),dim=c(1,1,ngroup))
  for(k in 1:ngroup){
    matrix.uti[[k]] = cbind(beta0_t[,k], beta1_t[,k], gamma0_t[,k], gamma1_t[,k], alpha_t[,k], rho_t)[ind, ]
  }  
  return(matrix.uti) 
}

main <- function(n.cohortsize, cohortsize, C_T, C_E, pt.mat,
                 prior.pt.mat, prior.pe.mat, vpt, vd, type){
  
  ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
  umat <- ed$umat
  vmat <- ed$vmat
  nd <- ed$nd
  vbeta <- ed$vbeta  # beta0.hat, beta1.hat
  vgamma <- ed$vgamma # gamma0.hat, gamma1.hat

  k.status <- rep(1, ngroup) # the initial status for each group, 1: continue going; 
  # 0:stop for too toxic or reached the max sample size
  dose.ind <- rep(1, ngroup)  # the starting dose levels for each group
  h.dose <- rep(1, ngroup)   # the current highest dose level has been tried for each group
  effd <- array(list(),dim=c(1,1,ngroup))
  dat <- array(list(),dim=c(1,1,ngroup))
  for(ich in 1:n.cohortsize){
    for(k in 1:ngroup){
      if(k.status[k]==1){
        u=umat[[k]]
        v=vmat[[k]]
        out.temp <- outcome(dose.ind[k], cohortsize, pt.mat[[k]], vpt[k,], u, v, type) # list(Y_T, Y_E)
        dat[[k]] <- rbind(dat[[k]], out.temp)
        # if(k==4){
        #   cat('data', k, '\n')
        #   print(dat[[k]])
        # }
      } 
    }
    mcmc.temp <- mcmc.arms(dat, vbeta=vbeta, vgamma=vgamma, ngroup=ngroup, type) # the return value of mcmc.arms should be a list like para[[k]]: which include 5 columns
    
    for(k in 1:ngroup){
      # op<-par(mfrow=c(2,2))
      # plot(mcmc.temp[[k]][,1], main = paste('beta0 = ', round(mean(mcmc.temp[[k]][,1]),2), ' k=', k))
      # plot(mcmc.temp[[k]][,2], main = paste('beta1 = ', round(mean(mcmc.temp[[k]][,2]),2)))        
      # plot(mcmc.temp[[k]][,3], main = paste('gamma0 = ', round(mean(mcmc.temp[[k]][,3]),2))) 
      # plot(mcmc.temp[[k]][,4], main = paste('gamma1 = ', round(mean(mcmc.temp[[k]][,4]),2)))
      # par(op)
      if(k.status[k]==1){
        u=umat[[k]]
        v=vmat[[k]]
        summ.temp<-summary.mcmc(matrix.uti=mcmc.temp[[k]],n.dose=nd[k], u=u, v=v, type)
        tox.summ <- summ.temp$tox
        cat('k:',k,'\n')
        cat('uti:', summ.temp$uti,'\n')
        if ((tox.summ[h.dose[k]]>.5) & (h.dose[k]!=nd[k])){ # escalation conditions
          dose.ind[k] <- h.dose[k]+1
          h.dose[k] <- dose.ind[k]
        }else{
          eff.summ <- summ.temp$eff
          capA.ind <- which((tox.summ>C_T)&(eff.summ>C_E))
          cat('tox:', tox.summ,'\n')
          cat('eff:', eff.summ,'\n')
          cat('admissible dose set:', capA.ind, '\n')
          if(length(capA.ind)==0){ # empty
            dose.ind[k] = 0
            k.status[k] = 0
          }else{# non empty
            uti <- summ.temp$uti[capA.ind]
            if(ich==n.cohortsize){ # reached the maximum sample size
              dose.ind[k] = capA.ind[which.max(uti)]
              k.status[k] = 0
            }else{
              # randomize patients within admissible doses
              uti.normalize <- uti/sum(uti)
              cumu.uti <- uti
              for (i in 1:length(uti)) cumu.uti[i] <- sum(uti.normalize[1:i])
              r <- runif(1,0,1)
              dose.ind[k] <- capA.ind[min(which(r<cumu.uti))]
              if(dose.ind[k]>h.dose[k]){
                dose.ind[k] <- h.dose[k]+1
                h.dose[k] <- dose.ind[k]
              }
            }
          }
        } 
      }
    }
    if(sum(k.status)==0) break
  }
  res = array(list(), dim=c(1,1,ngroup))
  for(k in 1:ngroup){
    res[[k]] = list(d.select = dose.ind[k], d=dat[[k]][,3])
  }
  return(res)
}

uti.fun <- function(pt, mu.E, rho, type){
  y_E.cut = 2
  n.y <- 10000
  y_T <- rbinom(n.y, 1, pt)
  Z_T <- ifelse(y_T==0, -1.3, 1.3)
  mu.T <- mean(Z_T)
  sample2d <- norm2d(v=c(mu.T, mu.E), rho=rho, n=n.y)
  Uti <- rbind(c(0, 50),
               c(25, 100))
  # V1
  # Uti <- rbind(c(0, 55), 
  #              c(20, 100))
  
  c_Y = ifelse(type=='binary', 0, y_E.cut)
  prob <- matrix(0,2,2)
  prob[1,1] <- sum((sample2d[,1]>0) & (sample2d[,2]<=c_Y))/n.y
  prob[1,2] <- sum((sample2d[,1]>0) & (sample2d[,2]>c_Y))/n.y
  prob[2,1] <- sum((sample2d[,1]<=0) & (sample2d[,2]<=c_Y))/n.y
  prob[2,2] <- sum((sample2d[,1]<=0) & (sample2d[,2]>c_Y))/n.y
  uti <- sum(Uti*prob)
  return(uti)
}

log_main<-function(nrep){
  pcs = pd = array(0,c(ngroup, ndoses, nrep, nscenario))
  for(isc in 1:nscenario){
    # true values of the parameters used in our models
    {# Scenarios used in simulation study
      if(type=='binary'){
        # Scenario 1
        #V0:      
        if(isc==1){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          vd[[1]] = c(0.50, 0.70, 0.80, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[3]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[4]] = c(0.10, 0.20, 0.30, 0.50)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[3]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[4]] = c(0.15, 0.30, 0.45, 0.55)
          
          pt.mat[[1]] = c(0.25, 0.45, 0.55, 0.60)
          pt.mat[[2]] = c(0.20, 0.30, 0.48, 0.55)
          pt.mat[[3]] = c(0.05, 0.20, 0.30, 0.48)
          pt.mat[[4]] = c(0.02, 0.05, 0.10, 0.20)
          
          vpt = rbind(c(-0.20, 1.9, 0, 0.5),
                      c(-0.30, 2.3, 0, 0.5),
                      c(-0.20, 2.3, 0, 0.5),
                      c(-0.20, 2.1, 0, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i])
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        # Scenario 2
        if(isc==2){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          vd[[1]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[3]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[4]] = c(0.10, 0.30, 0.50, 0.70)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[3]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[4]] = c(0.15, 0.30, 0.45, 0.55)
          
          
          pt.mat[[1]] = c(0.20, 0.30, 0.50, 0.60)
          pt.mat[[2]] = c(0.20, 0.30, 0.45, 0.60)
          pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.45)
          pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.50)
          
          vpt = rbind(c(-0.10, 2.5,0, 0.5),
                      c(-0.10, 2.5,0, 0.5),
                      c(-0.10, 2.5,0, 0.5),
                      c(-0.10, 2.5,0, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i])
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        # Scenario 3
        
        if(isc==3){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          vd[[1]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[2]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[3]] = c(0.50, 0.70, 0.90)
          vd[[4]] = c(0.10, 0.30, 0.50)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(0.15, 0.30, 0.45, 0.55)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[3]] = c(0.15, 0.30, 0.45)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[4]] = c(0.15, 0.30, 0.45)
          
          pt.mat[[1]] = c(0.20, 0.30, 0.45, 0.55) # too toxic
          pt.mat[[2]] = c(0.10, 0.15, 0.20, 0.25) # all futile
          pt.mat[[3]] = c(0.15, 0.30, 0.45) 
          pt.mat[[4]] = c(0.10, 0.15, 0.20) 
          
          vpt = rbind(c(-1.00, 1.7, 0.50, 0.5),
                      c(-1.00, 1.7,-0.30, 0.5),
                      c(-0.20, 1.7, 0.00, 0.5),
                      c(-0.20, 1.7, 0.50, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i])
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        
        # Scenario 4
        
        if(isc==4){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          # [3,]  0.1  0.15  0.20 0.25 0.30 0.35 0.40 0.45
          
          vd[[1]] = c(0.50, 0.70, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70)
          vd[[3]] = c(0.10, 0.30, 0.50)
          vd[[4]] = c(0.10, 0.30, 0.50)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[1]] = c(0.20, 0.30, 0.40)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[2]] = c(0.20, 0.30, 0.40)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[3]] = c(0.20, 0.30, 0.40)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[4]] = c(0.20, 0.30, 0.40)
          
          pt.mat[[1]] = c(0.25, 0.45, 0.55)
          pt.mat[[2]] = c(0.25, 0.30, 0.45) 
          pt.mat[[3]] = c(0.05, 0.10, 0.15) 
          pt.mat[[4]] = c(0.10, 0.20, 0.30)
          
          vpt = rbind(c(-0.20, 1.8, 0.00, 0.5),
                      c(-0.10, 1.5,-0.15, 0.5),
                      c(-0.30, 2.0, 0.00, 0.5),
                      c(-0.90, 1.5,-0.50, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i])
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        
      }else{
        y_E.cut = 2
        # Scenario 1
        if(isc==1){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd = array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          vd[[1]] = c(0.50, 0.70, 0.80, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[3]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[4]] = c(0.10, 0.20, 0.30, 0.50)
          
          # round(PI$vgamma,2)
          # [1] 2 3
          pt.mat[[1]] = c(0.30, 0.45, 0.50, 0.55)
          pt.mat[[2]] = c(0.20, 0.30, 0.45, 0.50)
          pt.mat[[3]] = c(0.15, 0.20, 0.30, 0.45)
          pt.mat[[4]] = c(0.05, 0.10, 0.15, 0.30)
          
          vpt = rbind(c(1.9, 2.0, 0.0, 0.5),
                      c(1.9, 2.5, 0.0, 0.5),
                      c(1.9, 2.5, 0.0, 0.5), 
                      c(1.9, 2.0, 0.0, 0.5))
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[3]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[4]] = c(1.50, 2.00, 2.50, 3.00)
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i]-y_E.cut)
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        # Scenario 2
        # Scenario 2
        if(isc==2){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd <- array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          vd[[1]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[3]] = c(0.10, 0.30, 0.50, 0.70)
          vd[[4]] = c(0.10, 0.30, 0.50, 0.70)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[3]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[4]] = c(1.50, 2.00, 2.50, 3.00)
          
          pt.mat[[1]] = c(0.20, 0.30, 0.50, 0.60)
          pt.mat[[2]] = c(0.20, 0.30, 0.50, 0.60)
          pt.mat[[3]] = c(0.10, 0.20, 0.30, 0.50)
          pt.mat[[4]] = c(0.10, 0.20, 0.30, 0.50)
          
          vpt = rbind(c(1.8, 3.0, -0.00, 0.5),
                      c(1.9, 3.0, -0.00, 0.5),
                      c(1.8, 2.5, -0.00, 0.5),
                      c(1.6, 2.5, -0.00, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i]-y_E.cut)
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        # Scenario 3
        if(isc==3){
          # beta: -1.23, 2.2 gamma: 2, 2  
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd <- array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          
          
          vd[[1]] = c(0.30, 0.50, 0.70, 0.90)
          vd[[2]] = c(0.10, 0.20, 0.30, 0.40)
          vd[[3]] = c(0.50, 0.70, 0.90)
          vd[[4]] = c(0.10, 0.30, 0.50)
          
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[1]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30, 0.40)
          prior.pe.mat[[2]] = c(1.50, 2.00, 2.50, 3.00)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[3]] = c(1.50, 2.00, 2.50)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[4]] = c(1.50, 2.00, 2.50)
          
          pt.mat[[1]] = c(0.20, 0.25, 0.30, 0.48) # too toxic
          pt.mat[[2]] = c(0.10, 0.15, 0.20, 0.25) # all futile
          pt.mat[[3]] = c(0.25, 0.30, 0.45) 
          pt.mat[[4]] = c(0.15, 0.20, 0.25)
          
          vpt = rbind(c(0.4, 2.0, 0.50, 0.5),
                      c(0.5, 2.0, -0.12, 0.5),
                      c(1.8, 2.0, -0.00, 0.5),
                      c(1.8, 2.0, -0.00, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i]-y_E.cut)
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
        
        # Scenario 4
        if(isc==4){
          ngroup = 4
          pt.mat = pe.mat = uti.mat = prior.pt.mat = prior.pe.mat = vd <- array(list(),dim=c(1,1,ngroup)) # initialize the prior set
          # [3,]  0.1  0.15  0.20 0.25 0.30 0.35 0.40 0.45
          prior.pt.mat[[1]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[1]] = c(1.50, 2.00, 2.50)
          
          prior.pt.mat[[2]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[2]] = c(1.50, 2.00, 2.50)
          
          prior.pt.mat[[3]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[3]] = c(1.50, 2.00, 2.50)
          
          prior.pt.mat[[4]] = c(0.10, 0.20, 0.30)
          prior.pe.mat[[4]] = c(1.50, 2.00, 2.50) 
          
          
          vd[[1]] = c(0.50, 0.70, 0.90)
          vd[[2]] = c(0.30, 0.50, 0.70)
          vd[[3]] = c(0.10, 0.30, 0.50)
          vd[[4]] = c(0.10, 0.30, 0.50)
          
          pt.mat[[1]] = c(0.30, 0.45, 0.55)
          pt.mat[[2]] = c(0.25, 0.35, 0.50) 
          pt.mat[[3]] = c(0.10, 0.15, 0.20) 
          pt.mat[[4]] = c(0.25, 0.30, 0.35)
          
          vpt = rbind(c(1.70, 2.5, -0.00, 0.5),
                      c(1.80, 1.0, -0.10, 0.5),
                      c(1.70, 2.0, -0.00, 0.5),
                      c(1.20, 2.0, -0.20, 0.5))
          
          ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
          for(k in 1:ngroup){
            u = ed$umat[[k]]
            v = ed$vmat[[k]]
            pe.mat[[k]] = uti.mat[[k]] <- rep(0, length(u))
            mu.E = vpt[k,1] + vpt[k,2]*ifelse(v<=vpt[k,3], v, vpt[k,3])
            for(i in 1:length(u)){
              pe.mat[[k]][i] = pnorm(mu.E[i]-y_E.cut)
              uti.mat[[k]][i]= uti.fun(pt.mat[[k]][i],mu.E[i],vpt[k,4],type)
            }
          }
        }
      }
      
      op<-par(mfrow=c(1,2),mar = c(4, 4, 4, 4) + 0.3, 
              cex=0.55)
      if(type=='binary'){
        c_urange = cbind(c(26, 49.1),
                         c(25, 45),
                         c(20, 45),
                         c(18, 51))
        c_erange = cbind(c(0, 0.50),
                         c(0, 0.45),
                         c(0, 0.50),
                         c(0, 0.50))
      }
      if(type=='continuous'){
        c_urange = cbind(c(25, 50),
                         c(20, 50),
                         c(15, 48),
                         c(18, 50))
        c_erange = cbind(c(0, 0.5),
                         c(0, 0.5),
                         c(0, 0.45),
                         c(0, 0.50))
      }
      prt_c = pre_c = uti_c<- NULL
      prt_c0 = pre_c0 = uti_c0<- NULL
      if(isc!=3){
        for(k in 1:4){
          prt_c <- cbind(prt_c,pt.mat[[k]])
          pre_c <- cbind(pre_c,pe.mat[[k]])
          uti_c <- cbind(uti_c,uti.mat[[k]])
        }
        matplot(prt_c, type = 'b', lty = 1, lwd = 0.5,
                col=1,
                ylim= c(min(c(unlist(pt.mat),unlist(pe.mat))), max(c(unlist(pt.mat), unlist(pe.mat)))),
                main=paste('Scenario', isc, sep=' '),
                xlab = 'Dose level', ylab = 'Probability',
                xaxt="n",mgp = c(2, 0.5, 0))
        abline(h=0.3, col='black', lty=3,lwd=2)
        # draw an axis
        if(isc==4){
          axis(1, at=1:3,labels=1:3, las=1)
        }else{
          axis(1, at=1:4,labels=1:4, las=1)
        }
        if(isc==1){
          legend("topleft",
                 legend = c('Toxicity', 'Utility'),
                 lty = 1,
                 col=c('black', 'chocolate3'),
                 lwd = 0.8,
                 seg.len = 1.3,
                 bty = 'n')
        }
        
        par(new=T)
        matplot(uti_c, type = 'b', ylim = c_urange[,isc],
                lty = 1, col = 'chocolate3',
                axes = FALSE, xlab = "", ylab = "", 
                lwd=0.5,
        )
        axis(side = 4, at = c(25, 30, 35, 40, 45, 50))# pretty(range(uti)))      # Add second axis
        mtext("Utility", side = 4, line = 2, cex = 0.55)
        
        ##############################
        matplot(pre_c, type = 'b', lty = 1, col = '1', 
                main=paste('Scenario', isc, sep=' '),
                xlab = 'Dose level',
                ylab = 'Probability',
                ylim = c_erange[,isc],
                mgp = c(2, 0.5, 0),
                xaxt = 'n',mgp = c(2, 0.5, 0))
        # draw an axis
        if(isc==4){
          axis(1, at=1:3,labels=1:3, las=1)
        }else{
          axis(1, at=1:4,labels=1:4, las=1)
        }
        abline(h=0.3, lty=3, lwd=2)
        if(isc==1){
          legend("topleft",
                 legend = 'Efficacy',
                 lty = 1,
                 col='black',
                 lwd = 0.8,
                 seg.len = 1.3,
                 bty = 'n')
        }
      }else{# isc = 3
        for(k in 1:2){
          prt_c <- cbind(prt_c,pt.mat[[k]])
          pre_c <- cbind(pre_c,pe.mat[[k]])
          uti_c <- cbind(uti_c,uti.mat[[k]])
        }
        for(k in 3:4){
          prt_c0 <- cbind(prt_c0,pt.mat[[k]])
          pre_c0 <- cbind(pre_c0,pe.mat[[k]])
          uti_c0 <- cbind(uti_c0,uti.mat[[k]])
        }
        matplot(prt_c, type = 'b', lty = 1, col=1, 
                ylim= c(min(c(unlist(pt.mat),unlist(pe.mat))), max(c(unlist(pt.mat), unlist(pe.mat)))),
                main=paste('Scenario', isc, sep=' '),
                xlab = 'Dose level',
                ylab = 'Probability',
                xaxt = 'n',mgp = c(2, 0.5, 0),
                lwd = 0.5)
        axis(1, at=1:4,labels=1:4, las=1)
        points(x=c(1,2,3),y=prt_c0[,1],type='b',pch='3', lwd=0.5)
        points(x=c(1,2,3),y=prt_c0[,2],type='b',pch='4', lwd=0.5)
        abline(h=0.3, col='black', lty=3, lwd=2)
        
        par(new=T)
        matplot(uti_c, type = 'b', lty = 1, col='chocolate3',
                main=paste('Scenario', isc, sep=' '), 
                ylim = c_urange[,isc], 
                axes = FALSE, xlab = "", 
                ylab = "", xaxt='n')
        points(x=c(1,2,3),y=uti_c0[,1],type='b',pch='3', col='chocolate3')
        points(x=c(1,2,3),y=uti_c0[,2],type='b',pch='4', col='chocolate3')
        axis(side = 4, at = c(25, 30, 35, 40, 45, 50))# pretty(range(uti)))      # Add second axis
        mtext("Utility", side = 4, line = 2, cex = 0.55)
        axis(1, at=1:4,labels=1:4, las=1)
        if(isc==3){
          legend("topleft",
                 legend = c('Toxicity', 'Utility'),
                 lty = 1,
                 col=c('black','chocolate3'),
                 lwd = 0.8,
                 seg.len = 1.3,
                 bty = 'n')
        }
        
        matplot(pre_c, type = 'b', lty = 1, col='black',
                main=paste('Scenario', isc, sep=' '),
                xlab = 'Dose level',
                ylab = 'Probability',
                xaxt = 'n',mgp = c(2, 0.5, 0),
                ylim = c_erange[,isc])
        points(x=c(1,2,3),y=pre_c0[,1],type='b',pch='3', col = 'black')
        points(x=c(1,2,3),y=pre_c0[,2],type='b',pch='4', col = 'black')
        axis(1, at=1:4,labels=1:4, las=1)
        abline(h=0.3, lty=3, lwd=2)
        
        if(isc==3){
          legend("topleft",
                 legend = 'Efficacy',
                 lty = 1,
                 col='black',
                 lwd = 0.8,
                 seg.len = 1.3,
                 bty = 'n')
        }
      }
      par(op)}
    
    ed <- effdose(prior.pt.mat, prior.pe.mat, ngroup, vd, type)
    umat <- ed$umat
    n.cohortsize = ifelse(isc==4, 14, 16)
    ####
    for(irep in 1:nrep){
      res = main(n.cohortsize, cohortsize=3, C_T=0.05, C_E=0.05, pt.mat, prior.pt.mat, prior.pe.mat, vpt, vd, type)
      for(k in 1:ngroup){
        u = umat[[k]]
        idselect <- res[[k]]$d.select
        d <- res[[k]]$d
        pcs[k,idselect,irep,isc] = pcs[k,idselect,irep,isc] + 1
        for(i in 1:length(d)){
          pd[k,which(u==d[i]),irep,isc] = pd[k,which(u==d[i]),irep,isc] + 1
        }
      }
    }
  }
  return(list(pcs=pcs,pd=pd))
}
# gloable settings
{
  phi.T = 0.3
  phi.E = 0.3
  nscenario = 4
  ngroup = 4  # K
  ndoses = 4  # J
  type = 'binary'
#  type = 'continuous'
  nrep = 200 # number of replication
}


tsart = Sys.time()
result = log_main(nrep)
save(result, file=paste('result',round(runif(1,0,10)*1000000,0),'.RData',sep = ''))



