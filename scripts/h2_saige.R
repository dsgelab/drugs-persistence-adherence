load("/home/cordioli/drugs/results/R5_GRM_V1-statins.rda")
names(modglmm)

tau <- modglmm$theta[2]

h2 <- tau/(pi^2/3+tau)
h2
