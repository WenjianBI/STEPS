#' @title Seoncdary Trait (ST) Analysis in Extreme Phenotype Sequencing (EPS)
#' @description A function to conduct genetic association analyses of secondary traits in extreme phenotype sequencing
#' @param data.mat An R matrix with each row for one subject.
#' Column 1 is primary phneotype, column 2 is secondary phenotype,
#' column 3 is genotype, the other columns are covariates. Missing data should be NA.
#' @param sec.type an R character to specify secodnary trait type: "binary" or "continuous"
#' @param ini.para initial parameters. Default is NULL, in which initial parameters would be estimated.
#' @param y0 lower cutoff value in EPS, more information can be seen in Details.
#' @param y1 lower cutoff value in EPS, more information can be seen in Details.
#' @param y2 upper cutoff value in EPS, more information can be seen in Details.
#' @param y3 upper cutoff value in EPS, more information can be seen in Details.
#' @return An R with the following components:
#' \item{pval}{p-value of Wald test to associate secondary trait with genotype}
#' \item{inv.fsh}{inverse matrix of fisher matrix}
#' \item{fnl.para}{an R list of a) res.opt: output of function optim(); b) res.lst: only estimated parameter}
#' @examples
#' ## First generate an parameter
#' par.ls = list(b0=0,b1=rnorm(2),b3=rnorm(1),g0=0,g1=rnorm(2))
#' ## Continuous Secondary Trait
#' data.cont = data.simu(par.ls,"continuous")
#' out=STEPS.snp(data.cont,"continuous")
#' out$pval
#' ## Binary Secondary Trait
#' data.bina = data.simu(par.ls,"binary")
#' out=STEPS.snp(data.bina,"binary")
#' out$pval
#' @details
#' Models to characterize Genotype, Primary Trait, Secondary Trait and Covariate can be seen in Details of help(data.simu).
#' Cutoff values y0<y1<y2<y3 are to specify the cutoff values in EPS. Subjects with primary trait Y between (y0,y1) or between (y2,y3) were retained to gentoype/sequence, and other subjects were removed.
#' If not specified, y0 is -Inf, y3 is Inf, y1 and y2 are estimated based on dataset. In the current version, we simply remove subjects with any missing data.
#' @export
STEPS.snp <- function(data.mat,
                      sec.type,
                      ini.para=NULL,
                      y0=-Inf,
                      y1=NULL,
                      y2=NULL,
                      y3=Inf)
{
  if(class(data.mat)!="matrix"|mode(data.mat)!="numeric") stop("class(data.mat) should be 'matrix' and mode(data.mat) should be 'numeric'.")
  if(sec.type!="binary"&sec.type!="continuous") stop("Argument 'sec.type' should be either 'binary' or 'continuous'.")
  row.na=which(apply(is.na(data.mat),1,any))
  if(length(row.na)!=0){
    data.mat=data.mat[-1*row.na,]
    warning(paste("Totally",length(row.na),"observations are removed because of data missing."))
  }

  if(sec.type=="continuous") res.snp = wald.snp.cont(data.mat,ini.para,y0,y1,y2,y3)
  if(sec.type=="binary"){
    D = data.mat[,2]   # column 2 is secodnary phenotype
    uD = sort(unique(D))
    if(length(uD)!=2) stop("Please check data input. If sec.type=='binary', secodnary phenotype (column 2) should be binary.")
    data.mat[,2] = ifelse(data.mat[,2]==uD[1],0,1)   # transform to 0 and 1
    res.snp = wald.snp.bina(data.mat,ini.para,y0,y1,y2,y3)
  }
  return(res.snp)
}


#################################################################
#
# -------- lower functions for binary secondary traits ----------
#
#################################################################

wald.snp.bina = function(data.mat,        # numeric matrix of data with one row per subject. Column 1 is primary phneotype, column 2 is secondary phenotype, column 3 is genotype, the other columns are covariates
                         ini.para=NULL,
                         y0=-Inf,
                         y1=NULL,
                         y2=NULL,
                         y3=Inf)
{
  fnl.para = para.est.bina(data.mat,ini.para,y0,y1,y2,y3)
  fsh.info = fisher.comp.bina(data.mat,y0,y1,y2,y3,fnl.para$res.lst)
  inv.fsh = solve(fsh.info)
  wald = fnl.para$res.opt$par[2]^2/inv.fsh[2,2]
  p.value = 1-pchisq(wald,df=1)
  names(p.value) = "pval"
  return(list(pval=p.value,inv.fsh=inv.fsh,fnl.para=fnl.para))
}

#####################################################################
#
# -------- lower functions for continuous secondary traits ----------
#
#####################################################################

wald.snp.cont = function(data.mat,
                         ini.para=NULL,
                         y0=-Inf,
                         y1=NULL,
                         y2=NULL,
                         y3=Inf)
{
  fnl.para = para.est.cont(data.mat,ini.para,y0,y1,y2,y3)
  fsh.info = fisher.comp.cont(data.mat,y0,y1,y2,y3,fnl.para$res.lst)
  inv.fsh = solve(fsh.info)
  wald = fnl.para$res.opt$par[2]^2/inv.fsh[2,2]
  p.value = 1-pchisq(wald,df=1)
  names(p.value) = "pval"
  return(list(pval=p.value,inv.fsh=inv.fsh,fnl.para=fnl.para))
}
