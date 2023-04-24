beta.greta.2r <- function(response, mm, r1, r2, pred.mm,
                          n = 1000){
  require(greta)
  require(arm)
  
  # set up beta parameters for mm
  beta <- normal(0, 10, dim=ncol(mm))
  
  # set up linear predictor
  linpred <- mm %*% beta
  
  # set up priors
  phi <- cauchy(0,5,truncation=c(0,Inf))
  
  # hierarchical model for site & plot effect
  # use the first level as the baseline like in lm()
  r1_sd <- lognormal(0, 10)
  r2_sd <- lognormal(0, 10)
  
  r1_offset <- normal(0, r1_sd, dim = length(unique(r1))-1)
  r1_effect <- rbind(0, r1_offset)
  r2_offset <- normal(0, r2_sd, dim = length(unique(r2))-1)
  r2_effect <- rbind(0, r2_offset)
  
  # set up beta shape parameters
  mu <- invlogit(linpred + r1_effect[r1] + r2_effect[r2])
  A <- mu * phi
  B <- (1 - mu) * phi
  
  # model
  distribution(response) <- beta(A, B)
  
  mod <- model(beta, r1_effect, r2_effect)
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