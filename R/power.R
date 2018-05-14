#' @title power.STEPS
#' @description A function to calculate power for STEPS given a set of parameters
#' @param maf minor allele frequency of SNP, between 0 and 0.5.
#' @param qntl quantile to choose y1 and y2, between 0 and 0.5. Default value is 0.1, that is, subjects with primary phenotype of top 10\% and bottom 10\% are in cohort
#' @param N sample size of dataset
#' @param b1 parameter to characterize association between genotype and primary trait. See 'Details' for more information.
#' @param b3 parameter to characterize association between secondary trait and primary trait. See 'Details' for more information.
#' @param g1 parameter to characterize association between genotype and secondary trait. See 'Details' for more information.
#' @return An R data frame with powers for both continuous/binary traits and 4 significance levels from 1E-5 to 1E-8.
#' @examples
#' power.STEPS(maf=0.3,qntl=0.1,N=1000,b1=-0.4,b3=-0.7,g1=0.3)
#' @details
#' For continuous secondary traits, model
#' \deqn{Z = g0+g1*G+0.4*X+e1}
#' \deqn{Y = b0+b1*G+0.4*X+b3*Z+e2}
#' For binary secondary traits, model
#' \deqn{Z = g0+g1*G+0.4*X+e1}
#' \deqn{D = I(Z>cutoff)}
#' \deqn{Y = b0+b1*G+0.4*X+b3*Z+e2}
#' where 'Z'/'D' is continuous/binary secondary trait, 'Y' is primary trait, 'X' is covariate following a standard normal distribution,
#' 'G' is genotype following HWE with MAF of 'maf', error term 'e1'/'e2' follows a standard normal distribution,
#' only subjects with primary phenotype at top/bottom quantile of 'qntl' are retained as extreme phenotype sampling design.
#' @export
power.STEPS=function(maf,qntl=0.1,N,b1,b3,g1){
  power.tot=c()
  for(i in c("continuous","binary")){
    for(j in c("1e-05","1e-06","1e-07","1e-08")){
      coef=power.coef.list[[i]][[j]]
      for(k in 1:nrow(coef)){
        x=coef$X.1[k]
        est=coef$Estimate[k]
        if(x=="(Intercept)"){
          log.power=est
        }else{
          value=paste(unlist(strsplit(x,":")),collapse = "*")
          log.power=log.power+eval(parse(text=paste0(est,"*",value)))
        }
      }
      power=exp(log.power)
      power.tot=rbind(power.tot,c(i,j,power))
    }
  }
  power.tot=data.frame(power.tot,stringsAsFactors = F)
  colnames(power.tot)=c("st.type","alpha","power")
  power.tot$alpha=as.numeric(power.tot$alpha)
  power.tot$power=as.numeric(power.tot$power)
  power.tot$power=ifelse(power.tot$power>1,1,power.tot$power)
  return(power.tot)
}

