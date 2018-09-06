likelihoods = read.table("likelihoods.tab", header = TRUE)  # read tab file
LRT <- -2*(likelihoods$lnull - likelihoods$lalt)
pval <- 1-pchisq(LRT, 1)
l_pvals <- cbind(likelihoods, pval, -log(pval))
write.table(l_pvals, file ="likelihoods.pvals.tab", row.names=F, col.names=F, quote=F) 
