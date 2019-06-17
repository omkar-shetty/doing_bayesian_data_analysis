

# Ex 5.1 ----------------------------------------------------------------------

# data from table 5.4
# tpos/tneg : event that the test result is positive/negative for the disease
# dpos/dneg : event that the subject is actually positive/negative for the disease

#Based on the example data
pr_dpos <- 0.001 #prior prob that the a random person has the disease

p_tpos_dpos <- 0.99 # P(test = + | disease = +)
p_tpos_dneg <- 0.05 # P(test = + | disease = -)

# 1st round - Prob (dpos | tpos) - Using Bayes rule

p_dpos_tpos <- pr_dpos * p_tpos_dpos /(p_tpos_dpos*pr_dpos + 
                                         p_tpos_dneg*(1-pr_dpos))
# 0.01943463

# 2nd round - Prob (dpos |tpos,tneg)

p_tneg_dpos <- 1 - p_tpos_dpos
p_tneg_dneg <- 1 - p_tpos_dneg
p_dpos_tneg <- p_dpos_tpos * p_tneg_dpos/(p_tneg_dpos*p_dpos_tpos +
                                            p_tneg_dneg*(1 - p_dpos_tpos))
#0.0002085862
