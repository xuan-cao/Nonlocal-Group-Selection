##
group_pmom_log_post <- function(X, ind, Y, n, alpha_1, alpha_2, r, tau, cal=TRUE){
  
  library(phangorn)
  
  col = which(colnames(X) %in% ind)
  if(length(col) == 0) return(-100000)
  
  X_Gk = as.matrix(X[, col])
  Gk_car = ncol(X_Gk)
  beta_initial = rep(0.5, ncol(X_Gk))
  sigma_initial = 0.5
  
  obj_fun = function(beta_sigma){
    beta_Gk = beta_sigma[1:Gk_car]
    sigma_2 = beta_sigma[Gk_car+1]
    
    if(sigma_2 <= 0) sigma_2 = 10^-5
    obj = -(n+Gk_car)/2*log(2*pi) - ((n+Gk_car)/2+r*Gk_car+alpha_1+1)*log(sigma_2) -
      (r*Gk_car+Gk_car/2)*log(tau) - Gk_car*log(dfactorial(2*r-1)) -
      t(Y-X_Gk%*%beta_Gk)%*%(Y-X_Gk%*%beta_Gk)/(2*sigma_2) - t(beta_Gk)%*%beta_Gk/(2*tau*sigma_2) -
      alpha_2/sigma_2 + 2*r*sum(log(abs(beta_Gk)))
    return(-obj)
  }
  
  estimates = nlm(obj_fun, c(beta_initial, sigma_initial))
  
  ###get log Laplace approximation of f(LZj_j; dj)
  get_log_laplace = function(obj) {
    beta_Gk_hat = obj$estimate[1:Gk_car]
    sigma_2_hat = obj$estimate[Gk_car+1]
    
    xx = t(X_Gk) %*% X_Gk
    xy = t(X_Gk) %*% Y
    xy_xb = Y - X_Gk%*%beta_Gk_hat
    
    V11 = xx / sigma_2_hat + 1/(tau*sigma_2_hat) * diag(Gk_car) + diag(2*r/beta_Gk_hat^2)
    
    V12 = -(1/sigma_2_hat^2)*(xx %*% beta_Gk_hat - xy) - beta_Gk_hat / (tau*sigma_2_hat^2)
    
    V22 = -(1/sigma_2_hat^2)*((n+Gk_car)/2 + r*Gk_car + alpha_1 + 1) + t(xy_xb) %*% xy_xb / (4*sigma_2_hat^3) + 
      t(beta_Gk_hat)%*%beta_Gk_hat/(tau*sigma_2_hat^3) + (2*alpha_2/sigma_2_hat^3)
    
    V = rbind(cbind(V11, V12), cbind(t(V12), V22))
    
    result = (Gk_car/2+1) *log(2*pi) - obj$minimum - 1/2 * sum(log(abs(eigen(V)$values)))
    
    return(result)
  }
  
  logposterior = get_log_laplace(estimates)
  
  if(cal==TRUE)return(logposterior = logposterior)
  if(cal==FALSE)return(est = estimates$estimate)
  
}
