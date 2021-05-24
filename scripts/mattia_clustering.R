#Matti Pirinen 13.5.2020
# rm(list=ls())

# Generate data
ind = rep(c(1,2),c(100,50)) #indicator of cluster; 1 = mean is 1.0, 2 = mean is 0.5.
# ind = rep(c(1,2),c(50,100))
# ind = rep(c(1,2),c(75,75))

n = length(ind)
mus = c(1, 0.5)
sigmas = c(0.15, 0.15) #SDs of the two clusters
x = rnorm(n, mus[ind], sigmas[ind]) #generate data (individual trajectory with 100 points around 1 and 50 around .5)


# Estimation starts from here
#Set parameter values
mus = c(1, 0.5) #two groups have fixed mean: the two doses we assume for the drug
sigmas = c(0.15, 0.15) #and fixed SD
n.iter = 1000 #Gibbs sampler iterations
n.burnin = 50 #Burnin period will be discarded from final results (iterations 1,...,n.burnin)

#Initialize variables for the sampler
gr = sample(c(1,2), prob = c(0.8, 0.2), size = n, replace = T) #initialize group memberships
prop = mean(gr == 1) #initialize proportion to their empirical estimate
sum.gr = matrix(0, ncol = 2, nrow = n) #col1 counts membership in group 1, col2 in group2
res.prop = rep(NA, n.iter - n.burnin) #posterior distribution of prop

#Gibbs sampler
for(ii in 1:n.iter){
    
    # theta ~ beta(0.5+sum(gr==1),0.5+sum(gr==2))
    prop = rbeta(1, 0.5 + sum(gr == 1), 0.5 + sum(gr == 2)) #prior for prop is Beta(0.5, 0.5)
    if(ii > n.burnin) res.prop[ii-n.burnin] = prop
    
    # likelihood = [ theta*p(x|N(1, 0.15)); theta*p(x|N(0.5,0.15)) ]
    loglk = cbind( log(prop) + dnorm(x, mus[1], sigmas[1], log = T), log(1-prop) + dnorm(x, mus[2], sigmas[2], log = T) )
    
    pr.1 = 1/(1 + exp(loglk[,2]-loglk[,1])) #probability of group 1 for each observation
    gr = 2 - rbinom(n, prob = pr.1, size =  1) #sample group membership, where pr.1 is probability of 1, and other outcome value is 2

    if(ii > n.burnin){
        sum.gr[gr == 1, 1] = sum.gr[gr == 1, 1] + 1
        sum.gr[gr == 2, 2] = sum.gr[gr == 2, 2] + 1
    }
}

layout(matrix(c(1,1,2,2), nrow = 1))
barplot(t(sum.gr)/(n.iter - n.burnin), col = c("blue","white"), main = "Membership in group 1 (in blue)" )
plot(density(x))
hist(res.prop, breaks = 40, col = "limegreen", main ="posterior of proportion of 1")

