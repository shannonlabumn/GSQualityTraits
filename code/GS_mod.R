GS_run <- function(p,g,t,model){
  g <- read_geno(filename = g, ploidy=4, map=T, min.minor.allele=5,
                 dominance=T)
  mod1 <- Stage1(file=p, traits=t, solver = "spats",spline=c("row","range"))
  broad_het <- mod1$fit

  mod2 <- Stage2(data=mod1$blue, vcov=mod1$vcov)
  gxev <- summary(mod2$vars)
  if(model=="additive"){
    nvh <- Stage2(data=mod1$blue, geno=g, vcov=mod1$vcov, non.add="none", silent=FALSE)

    hv <- nvh$vars
    h <- summary(nvh$vars)[2,2]
    prep1 <- blup_prep(data=mod1$blue, vcov=mod1$vcov, geno=g, vars=nvh$vars)
    gb1 <- blup(prep1, geno=g, what="BV")
    gv1 <- blup(prep1, geno=g, what="GV")

    prep2 <- blup_prep(data=mod1$blue, vcov=mod1$vcov, vars=mod2$vars)
    gv2 <- blup(prep2, what="GV")

    #gvplot <- merge()


    allR <- list(mod1=mod1, mod2=mod2,broadSenseHeritability=broad_het,Varcomp=gxev,Varcomph1=nvh,heritability=h, geno=g, prep=prep1, gebv=gb1, genoValueM=gv1,genoValueP=gv2)
    return(allR)

  }else{
    nvh1 <- Stage2(data=mod1$blue, geno=g, vcov=mod1$vcov, non.add="g.resid", silent=FALSE)# non dominance model
    nvh2 <- Stage2(data=mod1$blue, geno=g, vcov=mod1$vcov, non.add="dom", silent=FALSE)# dominance,non additive model


    na <- data.frame(non.add=c("g.resid","dom"), AIC=c(nvh1$aic,nvh2$aic))
    if(na[1,2] > na[2,2]){
      hv <- nvh2$vars
      h <- summary(nvh2$vars)[2,2]
      prep1 <- blup_prep(data=mod1$blue, vcov=mod1$vcov, geno=g, vars=nvh2$vars)

    }else{
      hv <- nvh1$vars
      h <- summary(nvh1$vars)[2,2]
      prep1 <- blup_prep(data=mod1$blue, vcov=mod1$vcov, geno=g, vars=nvh1$vars)
    }

    gb1 <- blup(prep1, geno=g, what="BV")
    gv1 <- blup(prep1, geno=g, what="GV")

    prep2 <- blup_prep(data=mod1$blue, vcov=mod1$vcov, vars=mod2$vars)
    gv2 <- blup(prep2, what="GV")

    #gvplot <- merge()


    allR <- list(mod1=mod1, mod2=mod2,broadSenseHeritability=broad_het,Varcomp=gxev,Varcomph1=nvh1,Varcomph2=nvh2,heritability=h, geno=g, prep=prep1, gebv=gb1, genoValueM=gv1,genoValueP=gv2)
    return(allR)
  }
}
