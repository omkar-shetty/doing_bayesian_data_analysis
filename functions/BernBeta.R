## Find the HDI of a probability density function by ICDF 
## (inverse cumulative density function) 
HDIofICDF <- function(ICDFname, credMass=0.95, tol=1e-8, ...) {
  # Arguments:
  # ICDFname is Râ€™s name for the inverse cumulative density function
  # of the distribution (such as qnorm, qbeta, etc.).
  # credMass is the desired mass of the HDI region.
  # tol is passed to Râ€™s optimize function.
  # Return value:
  # Highest density iterval (HDI) limits in a vector, that is, 
  # the HDI is the narrowest interval of that mass.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  # HDIofICDF(qbeta, shape1 = 30 , shape2 = 12)
  # Notice that the parameters of the ICDFname must be explicitly named;
  # e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  incredMass = 1.0 - credMass
  intervalWidth = function(lowTailPr, ICDFname, credMass, ...){
    ICDFname(credMass + lowTailPr, ...) - ICDFname(lowTailPr, ...)
  }
  optInfo = optimize(intervalWidth, c(0, incredMass), ICDFname = ICDFname,
                     credMass = credMass, tol = tol, ...)
  HDIlowTailPr = optInfo$minimum
  return(c(ICDFname(HDIlowTailPr, ...),
           ICDFname(credMass + HDIlowTailPr, ...)))
}
## Bayesian updating for Bernoulli likelihood and beta prior. 
BernBeta <- function(priorShape, data, credMass = 0.95){
  # Input:
  # prior: vector of parameter values for the prior beta distribution;
  # data: vector of 1's and 0's;
  # credMass: the probability mass of equal tailed credible interval (default 0.95).
  # Output:
  # posterior beta distribution;
  # creats a three-panel graph of prior, likelihood, and posterior with 
  # highest posterior density distribution.
  # Example:
  # postShape <= BernBeta(priorShape = c(1, 1), data = c(1, 0, 0, 1, 1))
  # Check for errors in input  arguments:
  if (length(priorShape) != 2){
    stop('prior must have two components.')}
  if (any(priorShape <= 0)){
    stop('prior components must be positive.')}
  if (any(data != 1 & data != 0)){
    stop('data must be a vector of 1s and 0s.')}
  if (credMass <=0 | credMass >= 1.0){
    stop('credMass must be between 0 and 1.')}
  # Rename the prior shape parameters for convenience
  a = priorShape[1]; b = priorShape[2]
  z = sum(data == 1)  # number of 1's in data
  N = length(data)    # number of flips in data
  # Compute the posterior shape parameters
  postShape = c(a + z, b + N -z)
  # Compute the evidence p(D)
  pData = beta(a + z, b + N - z) / beta(a, b)
  # Determine the limits of the highest density interval
  hpdLim = HDIofICDF(qbeta, shape1 = postShape[1], shape2 = postShape[2])
  # Plotting graphs
  binwidth = 0.005  # Arbitrary small value for comb on theta 
  theta = seq(from = binwidth / 2, to = 1 - binwidth / 2, by = binwidth)
  # Compute the likelihood of the data at each value of theta
  pDataGivenTheta = theta^z * (1 - theta)^(N-z)
  # Compute the prior at each value of theta
  pTheta = dbeta(theta, a, b)
  # Compute the posterior at each value of theta
  pThetaGivenData = dbeta(theta, a + z, b + N -z)
  # Open a widow with three panels
  windows(7, 10)  # creat window of specified width, height inches.
  layout(matrix(c(1, 2, 3), nrow = 3))  # 3x1 panels
  par(mar = c(3, 3, 1, 0), 
      mgp = c(2, 1, 0), 
      mai = c(.5, .5, .3, .1))  # margin specs
  maxY = max(c(pTheta, pThetaGivenData))  # max y for plotting on prior & posterior
  # Plot the prior
  plot(theta, pTheta, type = 'l', lwd = 3,
       xlim = c(0, 1), ylim = c(0, maxY), cex.axis = 1.2,
       xlab = bquote(theta), ylab = bquote(p(theta)), cex.lab = 1.5,
       main = bquote("Prior: " * 'beta(' * theta * ';' * .(a) * ',' * .(b) * ')'), 
       cex.main = 1.5)
  # Plot the likelihood: p(data|theta)
  plot( theta, pDataGivenTheta, type = "l", lwd = 3,
        xlim = c(0, 1), ylim = c(0, 1.1 * max(pDataGivenTheta)), cex.axis = 1.2, 
        xlab = bquote(theta), ylab = bquote("p(D|" * theta * ")"), cex.lab = 1.5,
        main = bquote('Likelihood: ' * "Data(z=" * .(z) * ",N=" * .(N) * ')'), 
        cex.main = 1.5)
  # Plot the posterior.
  plot(theta, pThetaGivenData, type = "l", lwd = 3,
       xlim = c(0, 1), ylim = c(0, maxY), cex.axis = 1.2,
       xlab = bquote(theta), ylab = bquote("p(" * theta * "|D)"),
       cex.lab = 1.5, 
       main = bquote('Posterior: ' * 
                       "beta(" * theta * ";" * .(a + z) * 
                       "," * .(b + N - z) * ")" * '; P(D)=' * .(signif(pData,3))),
       cex.main = 1.5)
  hpdHt = mean(c(dbeta(hpdLim[1], a + z, b + N - z), 
                 dbeta(hpdLim[2], a + z, b + N - z)))
  lines(c(hpdLim[1], hpdLim[1]), c(-0.5, hpdHt), type = "l", lty = 2, lwd = 1.5, col = 'red')
  lines(c(hpdLim[2], hpdLim[2]), c(-0.5, hpdHt), type = "l", lty=2, lwd = 1.5, col = 'red')
  lines(hpdLim, c(hpdHt, hpdHt), type = "l", lwd = 2, col = 'red')
  text(hpdLim[1], hpdHt, bquote(.(round(hpdLim[1], 3))),
       adj = c(1.2, -0.1), cex = 1.2)
  text(hpdLim[2], hpdHt, bquote(.(round(hpdLim[2], 3))),
       adj = c(-0.1, -0.1), cex = 1.2)
  # Return posterior 
  return (postShape)
} 
# End of function


## Bayesian updating for Bernoulli likelihood and beta prior. 
BernBeta_plot <- function(priorShape ,data, credMass = 0.95){
  # Input:
  # prior: vector of parameter values for the prior beta distribution;
  # data: vector of 1's and 0's;
  # credMass: the probability mass of equal tailed credible interval (default 0.95).
  # Output:
  # posterior beta distribution;
  # creats a three-panel graph of prior, likelihood, and posterior with 
  # highest posterior density distribution.
  # Example:
  # postShape <= BernBeta(priorShape = c(1, 1), data = c(1, 0, 0, 1, 1))
  pTheta <- list()
  pThetaGivenData <- list()
  for(i in 1:length(priorShape)){
    
    ps <- priorShape[[i]]
    # Rename the prior shape parameters for convenience
    a = ps[1]; b = ps[2]
    z = sum(data == 1)  # number of 1's in data
    N = length(data)    # number of flips in data
    # Compute the posterior shape parameters
    postShape = c(a + z, b + N -z)
    # Compute the evidence p(D)
    pData = beta(a + z, b + N - z) / beta(a, b)
    # Determine the limits of the highest density interval
    hpdLim = HDIofICDF(qbeta, shape1 = postShape[1], shape2 = postShape[2])
    # Plotting graphs
    binwidth = 0.005  # Arbitrary small value for comb on theta 
    theta = seq(from = binwidth / 2, to = 1 - binwidth / 2, by = binwidth)
    # Compute the likelihood of the data at each value of theta
    pDataGivenTheta = theta^z * (1 - theta)^(N-z)
    # Compute the prior at each value of theta
    pTheta[[i]] = data.frame(prob =dbeta(theta, a, b))
    # Compute the posterior at each value of theta
    pThetaGivenData[[i]] = data.frame(prob = dbeta(theta, a + z, b + N -z))
    # Open a widow with three panels
  }

  pThetaGivenData <- rbindlist(pThetaGivenData, idcol = T)
  setnames(pThetaGivenData, c('prior','prob'))
  pThetaGivenData[,theta := rep(theta, length(priorShape))]
  
  pThetaGivenData[,prior2 := ifelse(prior == 1,2,1)]
  
  post_plot <- ggplot(data = pThetaGivenData) + 
    geom_line(aes(x = theta ,y = prob, 
                  color = as.factor(prior2), group = as.factor(prior2))) +
    theme_bw() +
    ggtitle(paste0('Posterior based on ',length(data), ' data points')) +
    labs(x = 'theta', y = 'p(theta |data)')
  
  
  pTheta <- rbindlist(pTheta, idcol = T)
  setnames(pTheta, c('pr','prob'))
  pTheta[,prior := ifelse(pr == 1,'non-informative prior','informative prior')]
  pTheta[,theta := rep(theta, length(priorShape))]
  
  prior_plot <- ggplot(data = pTheta) + 
    geom_line(aes(x = theta ,y = prob, 
                  color = prior, group = prior),size = 1) +
    theme_bw() +
    ggtitle(paste0('Prior Distributions')) +
    labs(x = 'theta', y = 'p(theta | a,b)')
  
  
  return(list(prior_plot = prior_plot, post_plot = post_plot))
 
} 
# End of function