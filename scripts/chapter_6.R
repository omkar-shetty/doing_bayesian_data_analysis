source('functions/DBDA2E-utilities.R')
source('functions/BernBeta.R')


library(data.table)
library(ggplot2)

# Ex 6.1
post = BernBeta_plot( priorShape =list(c(1,1),c(1.5,8.5)) , data=c(rep(1,4),rep(0,5)) )
post$prior_plot
post$post_plot


likelihood <- function(n,y,theta){
  return((theta^y)*((1-theta)^(n-y)))
}
theta = seq(from = 0.01, to = 0.99, by = 0.01)

lik <- data.frame(theta = theta, likelihood = likelihood(100,30,theta))
ggplot(data = lik, aes(x = theta, y = likelihood)) + geom_line(col = 'black') +
  labs(x = 'theta', y = 'P(D|theta)') + theme_bw()

plot(theta, likelihood(9,4,theta), type = 'l')

success <- 0:100

plot(success, dbinom(success, size=100, prob=.444),type='l')
