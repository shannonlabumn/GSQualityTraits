crossVal2 <- function(blue,vcov, geno,vars,gv2, indv, trait){
  MBS <- NULL
  predA <- data.frame()

  set.seed(0928761)
  for(Rep in 1:5){

    N <- length(indv)
    folds <- split(sample(indv),cut(1:N,10))

    for (i in 1:10) {
      prep <- blup_prep(data=blue,
                        vcov=vcov,
                        geno=geno,
                        vars=vars,
                        mask=data.frame(id=folds[[i]]),)
      pred <- blup(prep, geno=geno, what="GV")
      MBS <- rbind(MBS, pred[pred$id %in% folds[[i]],])

      data1 <- merge(MBS, gv2, by="id")
      names(data1) <- c("id","MBS","MBS_r2", "blups", "blups_r2")
      data1$key <- trait
    }
    #return(data1)

    predA <- rbind(predA, data.frame(Rep = Rep, predA=cor(data1$MBS, data1$blups),  Acc=sqrt(cor(data1$MBS, data1$blups)/mean(data1$blups_r2)), trait=trait))

    pred_Data <- list(dtab=data1, cval=predA)

  }
  return(pred_Data)
}
