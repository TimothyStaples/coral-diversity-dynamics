beta.greta.0r <- function(response, mm, pred.mm,
                          n = 1000){
  require(greta)
  require(arm)
  
  # set up beta parameters for mm
  beta <- normal(0, 10, dim=ncol(mm))
  
  # set up linear predictor
  linpred <- mm %*% beta
  
  # set up priors
  phi <- cauchy(0,5,truncation=c(0,Inf))
  
  # set up beta shape parameters
  mu <- invlogit(linpred)
  A <- mu * phi
  B <- (1 - mu) * phi
  
  # model
  distribution(response) <- beta(A, B)
  
  mod <- model(beta)
  draws <- mcmc(mod, n_samples = n)
  
  draw.summary <- summary(draws)

  pred.linpred <- pred.mm %*% beta
  spline.y <- calculate(pred.linpred, values=draws)
  pred.y <-calculate(linpred, values=draws)
  
  return(list(parameters=draw.summary, 
              spline.preds=spline.y, 
              data.preds=pred.y, 
              model=draws))
  
}