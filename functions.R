#########################################################################
############################# generate data #############################
#########################################################################
gen_data <- function(n, # #of observations
                     g, # #of groups
                     g.size, # group size 
                     t, # first t groups of columns have nonzero coefficients
                     setting, case #, seeds ,alpha_1, alpha_2
){
  p = g*g.size # total columns
  t.p = t*g.size #total # of columns for nonzero coefficients
  
  ## true coefficients (beta)
  beta_0_gt = rep(0, p)
  if(setting == 1){
    #set.seed(seeds) 
    beta_0_gt[1:t.p] = runif(t.p, min=0.5, max=1.5)
  }
  
  if(setting == 2) beta_0_gt[1:t.p] = rep(1.5, t.p)
  
  if(setting == 3){
    #set.seed(seeds) 
    beta_0_gt[1:t.p] = runif(t.p, min=1.5, max=3)
  }
  
  if(setting == 4) beta_0_gt[1:t.p] = rep(3, t.p)
  
  ## true design matrix (X)
  if(case == 1){
    Sigma = diag(p)
    #set.seed(seeds) 
    X = mvrnorm(n, mu=as.vector(rep(0,p)), Sigma=diag(p))
  }
  
  if(case == 2){
    Sigma = matrix(0.5, p, p)
    diag(Sigma) = 1
    #set.seed(seeds) 
    X = mvrnorm(n, mu=as.vector(rep(0,p)), Sigma=Sigma)
  }
  
  if(case == 3){
    Sigma = matrix(NA, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] = 0.5^(abs(i-j))
      }
    }
    #set.seed(seeds)
    X = mvrnorm(n, mu=as.vector(rep(0,p)), Sigma=Sigma)
  }
  
  ## true Y
  # if (!require("pscl")) install.packages("pscl")
  # library(pscl)
  #sg = rigamma(1, alpha_1, alpha_2)
  #Y = mvrnorm(1, mu=as.vector(X%*%beta_0_gt), Sigma = sg*diag(n))
  
  #set.seed(seeds)
  Y = X%*%beta_0_gt + rnorm(n)*sqrt(1)
  #Y = X%*%beta_0_gt + rnorm(n,mean = 0, sd=sqrt(1/sigma_2))
  
  X2 = scale(X)
  Y2 = Y-mean(Y)
  Y2 = as.vector(Y2)
  
  #set.seed(seeds)
  n_test = sample(c(1:n), 40)
  x.train = X[-n_test,]
  y.train = Y[-n_test]
  
  x.test = X[n_test,]
  y.test = Y[n_test]
  
  #colnames(X) = rep(1:g, each=g.size)
  #data = list(Y = Y, beta_true = beta_0_gt, X = X, t=t, Sigma=Sigma)
  
  return(list(X.sd = X2, Y.sd = Y2, X = X, Y= Y
              , x.train=x.train, y.train=y.train, x.test=x.test, y.test=y.test
              , beta_true = beta_0_gt, t=t, Sigma=Sigma))

}

#########################################################################
############################## Evaluations ##############################
#########################################################################
Evaluation <- function(beta1, beta2){
  true.index <- which(beta1==1)
  false.index <- which(beta1==0)
  positive.index <- which(beta2==1)
  negative.index <- which(beta2==0)
  
  TP <- length(intersect(true.index,positive.index))
  FP <- length(intersect(false.index,positive.index))
  FN <- length(intersect(true.index,negative.index))
  TN <- length(intersect(false.index,negative.index))
  
  
  Precision <- TP/(TP+FP)
  if((TP+FP)==0) Precision <- 1
  Recall <- TP/(TP+FN)
  if((TP+FN)==0) Recall <- 1
  Sensitivity <- Recall
  Specific <- TN/(TN+FP)
  if((TN+FP)==0) Specific <- 1
  MCC.denom <- sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN)
  if(MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if((TN+FP)==0) MCC <- 1
  
  return(list(Precision=Precision,Recall=Recall,Sensitivity=Sensitivity,Specific=Specific,
              MCC=MCC,TP=TP,FP=FP,TN=TN,FN=FN))
}

result <-function(fit){
  #GAM = fit[-1,]; OBJ = fit[1,]
  
  GAM = fit$GAM; OBJ = fit$OBJ; 
  p = nrow(GAM)
  marg.gam = rep(0,p)
  for(u in 1:ncol(GAM)){
    marg.gam = marg.gam + GAM[,u]*exp(OBJ[u]-max(OBJ))
  }
  marg.gam = marg.gam / sum(exp(OBJ-max(OBJ)))
  gam0 = GAM[,which.max(OBJ)]
  #print(gam0)
  ind2 = which(gam0==1)
  post = exp(OBJ-max(OBJ))/sum(exp(OBJ-max(OBJ)))
  hppm = 1/sum(exp(OBJ-max(OBJ)))
  #print("# of Searched Models by S5");print(length(OBJ))
  #print("The MAP model is ")
  #print(which(gam0==1))
  #print(paste("with posterior probability",round(hppm,3) )) 
  return(list(hppm = which(gam0==1), hppm.prob = hppm, marg.prob = marg.gam,
              gam = GAM, obj = OBJ, post = post) )
}


result.evl.SSS = function(fit.sss, data, g, t
                          #r, alpha_1, alpha_2,delta
                          ){
  
  X = data$X.sd
  x.train = data$x.train
  x.test = data$x.test
  y.train = data$y.train
  y.test = data$y.test
  
  rt.sss = result(fit.sss)
  nonzer.gp.sss = rt.sss$hppm
  nonzer.clo.sss = which(colnames(X) %in% nonzer.gp.sss)
  
  train.all = cbind(y.train, x.train[,nonzer.clo.sss])
  colnames(train.all)[-1] = nonzer.clo.sss

  library(stats)
  fittedmodel.sss <- lm(y.train ~ ., data = as.data.frame(train.all))
  test.use = as.data.frame(x.test[, nonzer.clo.sss])
  colnames(test.use) = nonzer.clo.sss
  predfitted.sss <- predict(fittedmodel.sss, newdata = test.use)
  
  # n = nrow(x.train)
  # tau = min(1, n^(-1)*r^(2+2*delta))
  # colnames(x.train) = rep(1:g, each=g.size)
  # beta.est = rep(0, ncol(x.train))
  # 
  # ests = group_pmom_log_post(x.train, nonzer.gp.sss, y.train, 
  #                            n, alpha_1, alpha_2, r, tau, cal=FALSE) 
  # beta.est[nonzer.clo.sss] = ests[1:(length(ests)-1)]
  # sigma.est = ests[length(ests)]
  # 
  # predfitted.sss = x.test %*% beta.est + rep(sigma.est, nrow(x.test))
  #MSPE.sss = round(mean((y.test - predfitted.sss)^2), digits = 4)
  mspe_eval = mspe_eval(y.test, predfitted.sss)
  
  pred.beta.sss = rep(0, g)
  pred.beta.sss[nonzer.gp.sss] = 1
  true.beta = c(rep(1, t), rep(0, g-t))
  
  evl=Evaluation(true.beta, pred.beta.sss)
  
  return(list(eval.y=mspe_eval, beta.evl=evl))
}



result.evl.gplasso <- function(data, g, g.size, t, penalty){
  group=rep(1:g, each=g.size)
  
  #fit = grpreg(X=data$x.train, y=data$y.train, group, penalty=penalty)
  cvfit = cv.grpreg(data$x.train, data$y.train, group, penalty=penalty)
  
  beta = as.vector(coef(cvfit, lambda=cvfit$lambda.min))
  
  true_group = c(rep(1, t*g.size), rep(0, g*g.size-t*g.size))
  a = Evaluation(true_group, 1*(beta[-1] != 0))
  
  Y_hat = as.vector(predict(cvfit, data$x.test, type="response"
                            , lambda=cvfit$lambda.min))
  #MSPE = round(mean((data$y.test - Y_hat)^2), digits = 4)
  mspe_eval = mspe_eval(data$y.test, Y_hat)
  
  return(list(eval.y=mspe_eval, beta.evl=a))
}


result.evl.SSGL <- function(data, g, g.size, t
                            , lambda0seq=seq(1, 100, by=2)){
  group=rep(1:g, each=g.size)
  ## Now fit model for chosen lambda0 and lambda1 values
  ## This example relied on a well-chosen value of lambda0, which we won't know in general. 
  ## To solve this, one can use the cross-validation function as below:
  modSSGLcv = SSGLcv(Y=data$y.train, X=data$x.train
                     , lambda1=1, 
                     lambda0seq = lambda0seq,
                     groups = group,
                     nFolds = 10)
  ## In our experience, lambda0 sequences should go from around 1 to 100. 
  ## Now one can check what the chosen value of lambda0 is:
  modSSGL = SSGL(Y=data$y.train, X=data$x.train
                 , lambda1=1, lambda0=modSSGLcv$lambda0
                 , groups = group)
  
  true_group = c(rep(1, t*g.size), rep(0, g*g.size-t*g.size))
  a = Evaluation(true_group, 1*(modSSGL$beta != 0))
  
  Y_hat_ssgl =cbind(rep(1, nrow(data$x.test)), data$x.test)%*%as.vector(c(modSSGL$intercept, modSSGL$beta))
  #MSPE_ssgl = round(mean((data$y.test - Y_hat_ssgl)^2), digits = 4)
  mspe_eval = mspe_eval(data$y.test, Y_hat_ssgl)
  
  return(list(eval.y=mspe_eval, beta.evl=a))
}


## MSPE & normalized MSPE
mspe_eval <- function(y_true, y_pred){
  N = length(y_true)
  y_true_diff = y_true - mean(y_true)
  
  SSPE = sum((y_pred - y_true)^2)
  MSPE = SSPE / N
  RMSPE = sqrt(MSPE)
  
  STD_obs = sqrt(sum((y_true_diff)^2) / N)
  
  RSR = RMSPE / STD_obs
  
  return(list(MSPE=MSPE, RMSPE=RMSPE, RSR=RSR))
}


