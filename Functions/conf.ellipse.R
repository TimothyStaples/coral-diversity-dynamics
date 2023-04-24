conf.ellipse <- function(covariance, center, conf=0.95, n.points=100){
  
  # Derived from ordiellipse() code in vegan package
  scale = sqrt(qchisq(conf, 2))
  theta <- (0:n.points) * 2 * pi/n.points
  Circle <- cbind(cos(theta), sin(theta))
  Q <- chol(covariance, pivot = TRUE)
  o <- attr(Q, "pivot")
  
  return(t(center + scale * t(Circle %*% Q[, o])))
  }
