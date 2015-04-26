# Functions for printing genotypes during the MCMC analysis. Each individual gets their
# own file.

print_g<-function(k,geno.mat,tot.mat){
  rnames<-row.names(tot.mat)
  for(i in 1:nrow(tot.mat)){
    cat(k,geno.mat[i,],file=paste("./genotypes/",rnames[i],"_g-mcmc.out",sep=""),append=TRUE)
    cat("\n",file=paste("./genotypes/",rnames[i],"_g-mcmc.out",sep=""),append=TRUE)
  }
}
