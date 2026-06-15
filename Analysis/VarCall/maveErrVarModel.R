model
{
    for (m in 1:M) {
        ##f.bv[m] ~ dt(mu[m],tau[m],tdf) ## model for replicate values
        f.bv[m] ~ dnorm(mu[m],tau[m]) ## model for replicate values
        tau[m] <- prec.reps[m]*pow(gamma[batchM[m]],-2)
        mu[m] <- beta[batchM[m]] + eta[variantM[m]]*gamma[batchM[m]]
        ## measurement variability of replicates
        log(prec.reps[m]) <- inprod(delta[1:7],ns[m,1:7])
    }

    
   delta[1] ~ dnorm(0,0.01)  ## sd=10; var=100, prec=0.01
   for (n in 2:7) {
       ##delta[n] ~ dnorm(0,0.01)  ## sd=10; var=100, prec=0.01
       delta[n] ~ dnorm(0,0.0625)  ## sd=4; var=16, prec=0.0625
       ##delta[n] ~ ddexp(0,1)  ## sd=1; var=1, prec=1
       ##delta[n] ~ dgamma(1.0,0.1)  ## mu=10, sd=10
    }

    ## batch-specific random effects
    for (b in 1:B) {
        beta[b] ~ dnorm(mu.beta + slope.beta*gamma[b], prec.beta)
        gamma[b] ~ dlnorm(mu.gamma, prec.gamma)
    }

    for (v in 1:V) {
        ## variant-specific function effects:
        eta[v] ~ dnorm(alpha[del[v]], prec.eta[1])  ##LDA (see prec updates below)
    }

    for (v in 1:V){
        ##  Prior on Pathogenicity Status, Pr(D)
        pi.del[v] ~ dbeta(a.pp,b.pp)
        P[v,1] <- (1.0 - pi.del[v])
        P[v,2] <- pi.del[v]
        del[v] ~ dcat(P[v,1:2])
    }

    
    ## variability of batch coefficients in mean regression
    prec.beta <- pow(sigma2.beta, -1.0)
    sigma2.beta <- pow(sigma.beta, 2.0)
    sigma.beta ~ dnorm(0.0, cauch.prec.b)
    cauch.prec.b ~ dgamma(0.5, 0.5)
    ## variability of batch scaling in mean regression
    prec.gamma <- pow(sigma2.gamma, -1.0)
    sigma2.gamma <- pow(sigma.gamma, 2.0)
    sigma.gamma ~ dnorm(0.0, cauch.prec.g)
    cauch.prec.g ~ dgamma(0.5, 0.5)

    ## variability of variant coefficients in mean regression
    for (k in 1:2){
        prec.eta[k] <- prec.se1
        sigma.eta[k] <- pow(sigma2.eta[k], 0.5)
    }
 
## function component variances
    sigma2.eta[1] <-pow(prec.se1,-1)
    sigma2.eta[2] <- sigma2.eta[1] ##equal
    prec.se1 ~ dgamma(25.0, 1.0)
    ##prec.se2 ~ dgamma(0.5, 0.5)

    ## Prior on batch location and scale random effects means
    mu.beta ~ dunif(-10,10)
    slope.beta ~ dnorm(-1.0,1.0)
    mu.gamma ~ dnorm(0.0,1.0)
    
    prob <- mean(del[])-1
        
}

