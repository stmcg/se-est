mln.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n){
  
  h <- floor(0.25*n+1)
  j <- floor(0.5*n+1)
  k <- floor(0.75*n+1)
  
  scenario <- estmeansd:::get.scenario(min.val = min.val, q1.val = q1.val, med.val = med.val, 
                                       q3.val = q3.val, max.val = max.val)
  
  #-------------------------------------------------------------------------------
  #Box-Cox transformation
  boxcoxtrans <- function(lambda, x)   
  {
    if (lambda == 0) x <- log(x, base = exp(1)) else x <- (x^lambda - 1)/lambda
    return(x)
  }
  
  if (scenario == 'S1'){
    logL.lambda <- function(lambda)
    {
      aa <- boxcoxtrans(lambda, min.val)
      mm <- boxcoxtrans(lambda, med.val)
      bb <- boxcoxtrans(lambda, max.val)
      mu <- (4/(4+n^0.75))*(aa+bb)/2 + (n^0.75/(4+n^0.75))*mm
      sigma <- (bb-aa)/(2*qnorm((n-0.375)/(n+0.25)))
      logL <- log(dnorm(aa,mu,sigma)) + log(dnorm(mm,mu,sigma)) + log(dnorm(bb,mu,sigma)) + (j-2)*log(pnorm(mm,mu,sigma)-pnorm(aa,mu,sigma)) + (n-j-1)*log(pnorm(bb,mu,sigma)-pnorm(mm,mu,sigma))
      return(logL)
    }
  } else if (scenario == 'S2'){
    logL.lambda <- function(lambda)
    {
      qq1 <- boxcoxtrans(lambda, q1.val)
      mm <- boxcoxtrans(lambda, med.val)
      qq3 <- boxcoxtrans(lambda, q3.val)
      mu <- (0.7+0.39/n)*(qq1+qq3)/2 + (0.3-0.39/n)*mm
      sigma <- (qq3-qq1)/(2*qnorm((0.75*n-0.125)/(n+0.25)))
      logL <- log(dnorm(qq1,mu,sigma)) + log(dnorm(mm,mu,sigma)) + log(dnorm(qq3,mu,sigma)) + (h-1)*log(pnorm(qq1,mu,sigma)) + (j-h-1)*log(pnorm(mm,mu,sigma)-pnorm(qq1,mu,sigma)) + (k-j-1)*log(pnorm(qq3,mu,sigma)-pnorm(mm,mu,sigma)) + (n-k)*log(1-pnorm(qq3,mu,sigma))
      return(logL)
    }
  } else if (scenario == 'S3'){
    logL.lambda <- function(lambda)
    {
      aa <- boxcoxtrans(lambda,min.val)
      qq1 <- boxcoxtrans(lambda,q1.val)
      mm <- boxcoxtrans(lambda,med.val)
      qq3 <- boxcoxtrans(lambda,q3.val)
      bb <- boxcoxtrans(lambda,max.val)
      mu <- (2.2/(2.2+n^0.75))*(aa+bb)/2+(0.7-0.72/n^0.55)*(qq1+qq3)/2+(0.3+0.72/n^0.55-2.2/(2.2+n^0.75))*mm
      sigma <- (bb-aa)/(4*qnorm((n-0.375)/(n+0.25)))+(qq3-qq1)/(4*qnorm((0.75*n-0.125)/(n+0.25)))
      logL <- log(dnorm(aa,mu,sigma)) + log(dnorm(qq1,mu,sigma)) + log(dnorm(mm,mu,sigma)) + log(dnorm(qq3,mu,sigma)) + log(dnorm(bb,mu,sigma)) + (h-2)*log(pnorm(qq1,mu,sigma)-pnorm(aa,mu,sigma)) + (j-h-1)*log(pnorm(mm,mu,sigma)-pnorm(qq1,mu,sigma)) + (k-j-1)*log(pnorm(qq3,mu,sigma)-pnorm(mm,mu,sigma)) + (n-k-1)*log(pnorm(bb,mu,sigma)-pnorm(qq3,mu,sigma))
      return(logL)
    }
  }
  
  
  #the MLE of lambda
  lambda.hat <- suppressWarnings(optimize(f = logL.lambda, interval = c(0,10), tol = 0.001,maximum = TRUE)$maximum)
  if (lambda.hat < 1e-08){
    lambda.hat <- 0
  }
  
  #-------------------------------------------------------------------------------
  if (scenario == 'S1'){
    min.val.lambda <- boxcoxtrans(lambda.hat, min.val)
    med.val.lambda <- boxcoxtrans(lambda.hat, med.val)
    max.val.lambda <- boxcoxtrans(lambda.hat, max.val)
  } else if (scenario == 'S2'){
    q1.val.lambda <- boxcoxtrans(lambda.hat, q1.val)
    med.val.lambda <- boxcoxtrans(lambda.hat, med.val)
    q3.val.lambda <- boxcoxtrans(lambda.hat, q3.val)
  } else if (scenario == 'S3'){
    min.val.lambda <- boxcoxtrans(lambda.hat, min.val)
    q1.val.lambda <- boxcoxtrans(lambda.hat, q1.val)
    med.val.lambda <- boxcoxtrans(lambda.hat, med.val)
    q3.val.lambda <- boxcoxtrans(lambda.hat, q3.val)
    max.val.lambda <- boxcoxtrans(lambda.hat, max.val)
  }
  
  
  #-------------------------------------------------------------------------------
  #estimating mu and sigma after Box-Cox transformation
  if (scenario == 'S1'){
    logL <- function(theta){
      mu <- theta[1]
      sigma <- theta[2]
      logL <- dnorm(min.val.lambda,mu,sigma, log = TRUE) + dnorm(med.val.lambda,mu,sigma, log = TRUE) + dnorm(max.val.lambda,mu,sigma, log = TRUE) + (j-2)*log(pnorm(med.val.lambda,mu,sigma)-pnorm(min.val.lambda,mu,sigma)) + (n-j-1)*log(pnorm(max.val.lambda,mu,sigma)-pnorm(med.val.lambda,mu,sigma))
      return(-logL)
    }
  } else if (scenario == 'S2'){
    logL <- function(theta){
      mu <- theta[1]
      sigma <- theta[2]
      logL <- dnorm(q1.val.lambda,mu,sigma, log = TRUE) + dnorm(med.val.lambda,mu,sigma, log = TRUE) + dnorm(q3.val.lambda,mu,sigma, log = TRUE) + (h-1)*pnorm(q1.val.lambda,mu,sigma, log = TRUE) + (j-h-1)*log(pnorm(med.val.lambda,mu,sigma)-pnorm(q1.val.lambda,mu,sigma)) + (k-j-1)*log(pnorm(q3.val.lambda,mu,sigma)-pnorm(med.val.lambda,mu,sigma)) + (n-k)*log(1-pnorm(q3.val.lambda,mu,sigma))
      return(-logL)
    }
  } else if (scenario == 'S3'){
    logL <- function(theta){
      mu <- theta[1]
      sigma <- theta[2]
      logL <- dnorm(min.val.lambda,mu,sigma, log = TRUE) + dnorm(q1.val.lambda,mu,sigma, log = TRUE) + dnorm(med.val.lambda,mu,sigma, log = TRUE) + dnorm(q3.val.lambda,mu,sigma, log = TRUE) + dnorm(max.val.lambda,mu,sigma, log = TRUE) + (h-2)*log(pnorm(q1.val.lambda,mu,sigma)-pnorm(min.val.lambda,mu,sigma)) + (j-h-1)*log(pnorm(med.val.lambda,mu,sigma)-pnorm(q1.val.lambda,mu,sigma)) + (k-j-1)*log(pnorm(q3.val.lambda,mu,sigma)-pnorm(med.val.lambda,mu,sigma)) + (n-k-1)*log(pnorm(max.val.lambda,mu,sigma)-pnorm(q3.val.lambda,mu,sigma))
      return(-logL)
    }
  }
  
  
  #-------------------------------------------------------------------------------
  if (scenario == 'S1'){
    mean.LW <- (4/(4+n^0.75))*(min.val.lambda+max.val.lambda)/2 + (n^0.75/(4+n^0.75))*med.val.lambda
    sd.LW <- (max.val.lambda-min.val.lambda)/(2*qnorm((n-0.375)/(n+0.25)))
  } else if (scenario == 'S2'){
    mean.LW <- (0.7+0.39/n)*(q1.val.lambda+q3.val.lambda)/2 + (0.3-0.39/n)*med.val.lambda
    sd.LW <- (q3.val.lambda-q1.val.lambda)/(2*qnorm((0.75*n-0.125)/(n+0.25)))
  } else if (scenario == 'S3'){
    mean.LW <- (2.2/(2.2+n^0.75))*(min.val.lambda+max.val.lambda)/2+(0.7-0.72/n^0.55)*(q1.val.lambda+q3.val.lambda)/2+(0.3+0.72/n^0.55-2.2/(2.2+n^0.75))*med.val.lambda
    sd.LW <- (max.val.lambda-min.val.lambda)/(4*qnorm((n-0.375)/(n+0.25)))+(q3.val.lambda-q1.val.lambda)/(4*qnorm((0.75*n-0.125)/(n+0.25)))
  }
  
  #-------------------------------------------------------------------------------
  est.MLE <- try(suppressWarnings(optim(par = c(mean.LW,sd.LW), fn = logL)), silent = TRUE)
  if ('try-error' %in% class(est.MLE)){
    print(paste('Lambda', lambda.hat, 'mean.LW', mean.LW, 'sd.LW', sd.LW))
    est.MLE <- list(par = c(mean.LW,sd.LW))
  }
  mu.lambda <- est.MLE$par[1]
  sigma.lambda <- est.MLE$par[2]
  
  data.lambda <- rnorm(10000, mu.lambda, sigma.lambda)
  if (lambda.hat !=0) data.lambda <- data.lambda[(data.lambda > -1/lambda.hat) & (data.lambda < 2*mu.lambda+1/lambda.hat)]
  if (lambda.hat == 0) {data.origin <- exp(data.lambda)} else {data.origin<-(lambda.hat*data.lambda+1)^(1/lambda.hat)}
  
  #estimates by applying the MLN method
  mean.MLN <- mean(data.origin)
  sd.MLN <- sd(data.origin)
  
  return(list(est.mean = mean.MLN, est.sd = sd.MLN, 
              location = mu.lambda, scale = sigma.lambda, shape = lambda.hat))
}


