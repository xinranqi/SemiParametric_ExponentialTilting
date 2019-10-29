
library(haven)
library(lattice)
library(glmnet)
library(survival)
library(ggplot2)
library(gglasso)
library(grpreg)
library(stringr)

# data manipulation
#set.seed(123)
#aims_irgs_jan2018 <- read_sas("C:/Users/Administrator/Desktop/aims_irgs_whole.sas7bdat",NULL)
#irgs_indx <- read_sas("C:/Users/Administrator/Desktop/irgs_indx.sas7bdat",NULL)
#irgs_cov <- read_sas("C:/Users/Administrator/Desktop/covariates.sas7bdat",NULL)
#irgs_covP <- read_sas("C:/Users/Administrator/Desktop/pcovariates.sas7bdat",NULL)
#irgs_covD <- read_sas("C:/Users/Administrator/Desktop/dcovariates.sas7bdat",NULL)

################################################
# data manipulation for patient data subset
options(digits=4)
options(max.print = 2000)
aims_irgs_jan2018 <- read_sas("C:/Users/xinqi/Desktop/aims_irgs_whole.sas7bdat",NULL)
irgs_indx <- read_sas("C:/Users/xinqi/Desktop/irgs_indx.sas7bdat",NULL)
irgs_cov <- read_sas("C:/Users/xinqi/Desktop/covariates.sas7bdat",NULL)
irgs_covP <- read_sas("C:/Users/xinqi/Desktop/pcovariates.sas7bdat",NULL)
irgs_covD <- read_sas("C:/Users/xinqi/Desktop/dcovariates.sas7bdat",NULL)

agvhd2 <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
agvhd3 <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
anc500 <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
cgvhd <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
dead <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
lfs <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
plate20 <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
relapse <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)
trm <- read_sas("C:/Users/xinqi/Desktop/p.sas7bdat",NULL)


# data manipulation for paired data subsets
# options(digits=4)
# options(max.print = 2000)
# aims_irgs_jan2018 <- read_sas("C:/Users/xinqi/Desktop/aims_irgs_whole.sas7bdat",NULL)
# irgs_indx <- read_sas("C:/Users/xinqi/Desktop/irgs_indx.sas7bdat",NULL)
# irgs_cov <- read_sas("C:/Users/xinqi/Desktop/covariates.sas7bdat",NULL)
# irgs_covP <- read_sas("C:/Users/xinqi/Desktop/pcovariates.sas7bdat",NULL)
# irgs_covD <- read_sas("C:/Users/xinqi/Desktop/dcovariates.sas7bdat",NULL)

# agvhd2 <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# agvhd3 <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# anc500 <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# cgvhd <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# dead <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# lfs <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# plate20 <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# relapse <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)
# trm <- read_sas("C:/Users/xinqi/Desktop/paired.sas7bdat",NULL)


paired = (irgs_indx[which(irgs_indx$data_pair %in% c("Ponly","paired")),-c(1,2)])[, seq(from=(185*3 + 1), to=(185*5), by=1)]
pairedCompleted = irgs_indx[which(irgs_indx$data_pair %in% c("Ponly","paired")),-c(1,2)]
pairedcov = irgs_covP[which(irgs_covP$data_pair %in% c("Ponly","paired")),-1]
pairedFullInfo = aims_irgs_jan2018[which(aims_irgs_jan2018$data_pair %in% c("Ponly","paired")),]
OutcObsN = dim(pairedFullInfo)[1] - c(sum(is.na(pairedFullInfo$agvhd2)),
                                      sum(is.na(pairedFullInfo$agvhd3)),
                                      sum(is.na(pairedFullInfo$anc500)),
                                      sum(is.na(pairedFullInfo$cgvhd)),
                                      sum(is.na(pairedFullInfo$dead)),
                                      sum(is.na(pairedFullInfo$lfs)),
                                      sum(is.na(pairedFullInfo$plate20)),
                                      sum(is.na(pairedFullInfo$relapse)),
                                      sum(is.na(pairedFullInfo$trm)))

# Function to remove zero indicator variables
Exclude = function(DATASET,cutoff,seq) {
  Nobs = OutcObsN[seq]
  DATASET = DATASET[which(as.numeric(DATASET$Frequency) < (Nobs - cutoff)),1]
  DATASET
}
agvhd2Minor = Exclude(agvhd2,50,1)
agvhd3Minor = Exclude(agvhd3,50,2)
anc500Minor = Exclude(anc500,50,3)
cgvhdMinor = Exclude(cgvhd,50,4)
deadMinor = Exclude(dead,50,5)
lfsMinor = Exclude(lfs,50,6)
plate20Minor = Exclude(plate20,50,7)
relapseMinor = Exclude(relapse,50,8)
trmMinor = Exclude(trm,50,9)

IrgLabels = colnames(irgs_indx)[-c(1,2)]
IrgLTrans = NULL
for (i in 1:length(IrgLabels)) {
  IrgLTrans[i] = paste("Table", IrgLabels[i], sep=" ")
}
Indxagvhd2 = which(IrgLTrans %in% agvhd2$Table)
Indxagvhd3 = which(IrgLTrans %in% agvhd3$Table)
Indxanc500 = which(IrgLTrans %in% anc500$Table)
Indxcgvhd = which(IrgLTrans %in% cgvhd$Table)
Indxdead = which(IrgLTrans %in% dead$Table)
Indxlfs = which(IrgLTrans %in% lfs$Table)
Indxplate20 = which(IrgLTrans %in% plate20$Table)
Indxrelapse = which(IrgLTrans %in% relapse$Table)
Indxtrm = which(IrgLTrans %in% trm$Table)

# W dataframes
Wmatrx_paired = paired

# plot density plot of NA's for irg_p, irg_d and mismatch in the paired dataset
NAcount_paired = numeric(dim(Wmatrx_paired)[1])
for(i in 1:(dim(Wmatrx_paired)[1])) {
  NAcount_paired[i] = sum(is.na(Wmatrx_paired[i,]))
}
summary(NAcount_paired)
qplot(NAcount_paired, geom="histogram",binwidth=5,alpha=I(1),
      main="Histogram of Paired Data Subset's irg NA",
      xlab="Number of NA",xlim=c(0,120),ylim=c(0,150))



# ############# Imputation #############
# frequency of "0,1" of each mismatch w.r.t irg
for (i in 1:(dim(Wmatrx_paired)[2])) {
  irg_misi = Wmatrx_paired[,i]
  obs = (dim(Wmatrx_paired)[1]) - sum(is.na(irg_misi))
  zeros = sum(irg_misi==0,na.rm=TRUE) / obs
  ones  = sum(irg_misi==1,na.rm=TRUE) / obs
  irg_misi[is.na(irg_misi),1] = sample(c(0,1),size=sum(is.na(irg_misi)),
                                       prob=c(zeros,ones),replace=TRUE)
  Wmatrx_paired[,i] = irg_misi
}
# complete data set using imputation
FullInfo_paired = cbind.data.frame(pairedFullInfo,Wmatrx_paired)


# ########## Standardization ###########
Standardization = function(DATA) {
  STANDARDIZED = matrix(0,nrow=(dim(DATA)[1]),ncol=(dim(DATA)[2]))
  for (i in 1:(dim(DATA)[2])) {
    dummy_i = DATA[,i]
    mean_i = sum(dummy_i) / (dim(DATA)[1])
    if (mean_i>0) {
      STANDARDIZED[,i] = data.matrix((dummy_i - mean_i) / sqrt(mean_i*(1-mean_i)))
    }
  }
  STANDARDIZED = data.frame(STANDARDIZED)
  names(STANDARDIZED) = colnames(DATA)
  STANDARDIZED
}
# for group LASSO, we do not use standardization for the SNPs index matrix
# Wmatrx_paired_st = Standardization(Wmatrx_paired)

# complete data set using imputation & standardization data frame
FullInfo_paired_st = cbind.data.frame(pairedFullInfo, Wmatrx_paired)


# ############### LASSO ################
LASSOFUNC1 = function(XCovariates,TIMES,EVENTS,PENALTY,LAMBDAseq) {
  set.seed(123)
  irg_time = TIMES
  irg_event = EVENTS
  XCovariates = data.matrix(XCovariates)
  IndxNA = unique(c(which(is.na(irg_time)),which(is.na(irg_event))))
  cvfit = cv.glmnet(XCovariates[-IndxNA,],
                    Surv(time=(irg_time[-IndxNA]),event=(irg_event[-IndxNA])), 
                    family="cox",standardize=TRUE,
                    lambda=LAMBDAseq,
                    penalty.factor=PENALTY,
                    maxit=100000)
  plot(cvfit)
  
  # coefficients using the min rule
  Coefficients = coef(cvfit,s=cvfit$lambda.min)
  Active.Index = which(Coefficients!=0)
  Active.Coefficients = Coefficients[Active.Index]
  cat("The active indexes of irg are",Active.Index,"\n",
      "With corresponding coefficients:",Active.Coefficients,"\n",
      "The tuning parameter is:",cvfit$lambda.min,"\n")
  Active.Index
}



# ############### Group LASSO (Logistic Regression Model) ################
GroupLASSO = function(XCovariates,TIMES,EVENTS,PENALTY,GROUP) {
  set.seed(123)
  
  irg_time = TIMES
  irg_event = ifelse(EVENTS==0, -1, 1)
  XCovariates = data.matrix(XCovariates)
  IndxNA = unique(c(which(is.na(irg_time)),which(is.na(irg_event))))
  
  gglassofit = cv.gglasso(x=XCovariates[-IndxNA,], y=irg_event[-IndxNA], group=GROUP,
                          loss = "logit", pred.loss = "misclass",
                          nlambda = 1000,
                          lambda.factor = ifelse(nrow(XCovariates) < ncol(XCovariates), 0.05, 0.001),
                          lambda = NULL,
                          pf = PENALTY, nfolds=5,
                          eps = 1e-05, maxit = 1e+06, delta=1, intercept=TRUE)
                        
  plot(gglassofit) # plots the coefficients against the log-lambda sequence

  Coefficients = coef(gglassofit$gglasso.fit, gglassofit$lambda.1se)
  Active.Index = which((Coefficients[-1])!=0)
  Active.Coefficients = Coefficients[Active.Index]
  cat("The active indexes of irg are", Active.Index, ".", "\n",
      "With corresponding coefficients:", Active.Coefficients, "." ,"\n",
      "The tuning parameter is:", gglassofit$lambda.1se, ".", "\n")
  Active.Index
}

GroupLASSOCox(COVARIATE7,intxcgvh,FullInfo_paired$cgvhd,
              penaltyfactor7,group7,TRUE,cgvhdMinor)

# ############### Group LASSO (Cox Regression Model) ################
GroupLASSOCox = function(XCovariates,TIMES,EVENTS,PENALTY,GROUP, SubGroupLASSO,MinorIRGs) {
  set.seed(123)
  
  irg_time = TIMES
  irg_event = EVENTS
  # irg_event = ifelse(EVENTS==0, -1, 1)
  XCovariates = data.matrix(XCovariates)
  IndxNA = unique(c(which(is.na(irg_time)),which(is.na(irg_event))))

  gglassofit = cv.grpsurv(X=XCovariates[-IndxNA,], y=Surv(time=(irg_time[-IndxNA]),event=(irg_event[-IndxNA])),
                          group=GROUP,
                          penalty="grLasso", 
                          gamma=3, 
                          alpha=1,
                          nlambda=1000, 
                          lambda = seq(from=0.005, to=0.5, by=0.005),
                          lambda.min={if (nrow(XCovariates[-IndxNA,]) > ncol(XCovariates[-IndxNA,])) 0.001 else .05},
                          eps=.001, 
                          max.iter=10000, 
                          gmax=length(unique(GROUP)), 
                          tau=1/3,
                          group.multiplier=PENALTY,
                          nfolds=10, seed=123,
                          se='bootstrap', returnY=FALSE, trace=FALSE)
                       
  # plots the coefficients against the log-lambda sequence
  plot(gglassofit) 
  print(summary(gglassofit))
  
  # Beta at minimum CVE
  Coefficients = coef(gglassofit) 
  Active.Index = which(Coefficients!=0)
  Active.Coefficients = Coefficients[Active.Index]
  cat("The tuning parameter (lambda.min) is:", gglassofit$lambda.min, ".", "\n",
      "With corresponding coefficients:", "\n")
  print(Active.Coefficients)
  

  # performing sub-group LASSO to reduce the number of covariates
  # if there are more than 3 IRG SNPs selected in the first step group LASSO, then perform the second step group LASSO
  if(SubGroupLASSO==TRUE && (length(Active.Index) > (sum(GROUP <= sum(1-PENALTY)) + 6))) {
    
    # find indexes for IRG SNPs with minor frequency (<50)
    IrgLabels.Active = rownames(as.data.frame(Active.Coefficients))
    IrgLTrans.Active = NULL
    for (i in 1:length(IrgLabels.Active)) {
      IrgLTrans.Active[i] = paste("Table", IrgLabels.Active[i], sep=" ")
    }
    MinorIndx.Active = which(IrgLTrans.Active %in% MinorIRGs$Table)
    cat("Based on the group LASSO output (among all IRG SNPs), IRG SNPs with minor frequency are", 
        rownames(as.data.frame(Active.Coefficients))[MinorIndx.Active], ".","\n")
    
    # if there is IRG SNPs with minor frequency, then combine levels in each group; if NOt, just print out group LASSO results in the first step
    if (length(MinorIndx.Active) > 0) {
      # since SNPs with minor frequency are NOT grouped (only some levels from the same group are of minor frquency), we then map them to the overall IRG SNPs variables 
      MinorIRGinAllIndex = which(colnames(Wmatrx_paired) %in% rownames(as.data.frame(Active.Coefficients))[MinorIndx.Active])
      
      # after mapping to the overall IRG SNPs, we can then complete all the group levels for those incomplete group levels (this is neccessary for the second group LASSO)
      completePatient = NULL

      # complete IRG SNPs groups using binary trees
      for (i in 1:length(MinorIRGinAllIndex)) {
        incompleteIndx = MinorIRGinAllIndex[i] 
        IRGmodulusIndex = ifelse((incompleteIndx %% 185) == 0, 185, incompleteIndx %% 185)
        completePatient = c(completePatient, IRGmodulusIndex)
        completePatient = unique(sort(completePatient, decreasing=FALSE))
      }
      
      # complete IRG SNPs indexes
      CompleteIRGIndexs =  sum(GROUP <= sum(1-PENALTY)) + c(completePatient, (185 + completePatient))
      
      # IRG SNPs prepared to be combined
      MinorIRGforPatient = XCovariates[ , sum(GROUP <= sum(1-PENALTY)) + completePatient] + 
                           XCovariates[ , sum(GROUP <= sum(1-PENALTY)) + 185 + completePatient]
      CombinedMinorIRG = MinorIRGforPatient
      
      if (dim(CombinedMinorIRG)[2] == 1) {
        CombinedMinorIRG = as.matrix(as.vector(CombinedMinorIRG), ncol=1)
        colnames(CombinedMinorIRG) = paste0("irg_pc",completePatient)
      } else {
        colnames(CombinedMinorIRG) = paste0("irg_pc",completePatient)
      }

      `%notin%` <- Negate(`%in%`)
      # indexes of IRG SNPs with higher frequency
      HighFreqIRGsIndx = (Active.Index[-seq(from=1, to=sum(GROUP <= sum(1-PENALTY)), by=1)])[(Active.Index[-seq(from=1, to=sum(GROUP <= sum(1-PENALTY)), by=1)] %notin% CompleteIRGIndexs)]
      
      # combining IRG SNPs with minor frequency
      MinorIRGCombinedXCovariates = cbind.data.frame(XCovariates[ , c(seq(from=1, to=sum(GROUP <= sum(1-PENALTY)), by=1), HighFreqIRGsIndx)], CombinedMinorIRG)
      
      # group vector for sub-group LASSO
      MinorIRGCombinedGROUP = as.numeric(as.factor(c(GROUP[c(seq(from=1, to=sum(GROUP <= sum(1-PENALTY)), by=1), HighFreqIRGsIndx)], seq(from=9999, to=(9998+ncol(CombinedMinorIRG)), by=1))))
      
      # penalty vector for sub-group LASSO
      MinorIRGCombinedPENALTY = c(numeric(sum(1-PENALTY)), 
                                  rep(1, times=sum(unique(MinorIRGCombinedGROUP[MinorIRGCombinedGROUP>sum(1-PENALTY)])>0)))
      set.seed(123)
      SubGGLasso = cv.grpsurv(X=MinorIRGCombinedXCovariates[-IndxNA, ], y=Surv(time=(irg_time[-IndxNA]),event=(irg_event[-IndxNA])),
                              group=MinorIRGCombinedGROUP,
                              penalty="grLasso", 
                              gamma=3, 
                              alpha=1,
                              nlambda=1000, 
                              lambda = seq(from=0.005, to=0.5, by=0.005),
                              lambda.min={if (nrow(XCovariates[-IndxNA,]) > ncol(XCovariates[-IndxNA,])) 0.001 else .05},
                              eps=.001, 
                              max.iter=10000, 
                              gmax=length(unique(GROUP)), 
                              tau=1/3,
                              group.multiplier=MinorIRGCombinedPENALTY,
                              nfolds=10, seed=123,
                              se='bootstrap', returnY=FALSE, trace=FALSE)
      
      plot(SubGGLasso) # plots the coefficients against the log-lambda sequence
      print(summary(SubGGLasso))
      
      # Beta at minimum CVE
      SubCoefficients = coef(SubGGLasso) 
      SubActive.Index = which(SubCoefficients!=0)
      SubActive.Coefficients = SubCoefficients[SubActive.Index]
      cat("The tuning parameter (lambda.min) is:", SubGGLasso$lambda.min, ".", "\n",
          "With corresponding coefficients:", "\n")
      print(SubActive.Coefficients)
      
    } else {
      print("There is NOT any IRG SNP with minor frequency. The group LASSO output is kept the same as the first step.")
    }
    
  }
}



# ##############################################
# for paired data subset: mis-match vs irg_p vs irg_d
attach(FullInfo_paired)

# first outcome
COVARIATE1 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$disease10,pairedcov$disease20,pairedcov$disease30,pairedcov$disease40,
  pairedcov$disease50,pairedcov$disease100,
  pairedcov$ragecat21,pairedcov$ragecat22,
  pairedcov$cmvmatch1,pairedcov$cmvmatch2,pairedcov$cmvmatch3,pairedcov$cmvmatch4,
  pairedcov$condint1,pairedcov$condint2,pairedcov$condint3,
  pairedcov$tbidosegp1,pairedcov$tbidosegp2,pairedcov$tbidosegp3,
  pairedcov$tbidosegp4,
  pairedcov$firsttx1,
  pairedcov$Pcluster1, pairedcov$Pcluster2, pairedcov$Pcluster3, pairedcov$Pcluster4,
  pairedcov$Pcluster5, pairedcov$Pcluster6, pairedcov$Pcluster7,
  Wmatrx_paired)
group1 = c(rep(1, times=5), rep(2, times=4), rep(3, times=6), rep(4, times=2),
           rep(5, times=4), rep(6, times=3), rep(7, times=4), 8, rep(9, times=7),
           rep(seq(from=(9+1), to=(9+185), by=1), times=2))
penaltyfactor1 = c(rep(0, times=9), rep(1, times=185))

Pairedinx_1 = GroupLASSOCox(COVARIATE1,intxsurv,FullInfo_paired$dead,
                            penaltyfactor1,group1,TRUE,deadMinor)


# second outcome
COVARIATE2 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$cmvmatch1,pairedcov$cmvmatch2,pairedcov$cmvmatch3,pairedcov$cmvmatch4,
  pairedcov$graftype21,pairedcov$graftype22,
  pairedcov$tbidosegp1,pairedcov$tbidosegp2,pairedcov$tbidosegp3,pairedcov$tbidosegp4,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$t_deplete1,pairedcov$t_deplete0,
  pairedcov$cmldxtx1,pairedcov$cmldxtx2,
  pairedcov$disease10,pairedcov$disease20,pairedcov$disease30,pairedcov$disease40,
  pairedcov$disease50,pairedcov$disease100,
  pairedcov$gvhdprop21,pairedcov$gvhdprop22,pairedcov$gvhdprop23,pairedcov$gvhdprop24,
  pairedcov$gvhdprop27,pairedcov$gvhdprop28,
  pairedcov$indxtxgp1,pairedcov$indxtxgp2,pairedcov$indxtxgp3,pairedcov$indxtxgp4,
  Wmatrx_paired)
group2 = c(rep(1, times=5), rep(2, times=4), rep(3, times=2), rep(4, times=4),
           rep(5, times=4), rep(6, times=2), rep(7, times=2), rep(8, times=6),
           rep(9, times=6), rep(10, times=4),
           rep(seq(from=(10+1), to=(10+185), by=1), times=2) )
penaltyfactor2 = c(rep(0, times=10), rep(1, times=185))

Pairedinx_2 = GroupLASSOCox(COVARIATE2,intxrel,FullInfo_paired$lfs,
                            penaltyfactor2,group2,TRUE,lfsMinor)


# third outcome
COVARIATE3 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$disstat1,pairedcov$disstat2,pairedcov$disstat3,
  pairedcov$kps0,pairedcov$kps1,
  pairedcov$cmvmatch1,pairedcov$cmvmatch2,pairedcov$cmvmatch3,pairedcov$cmvmatch4,
  pairedcov$graftype21,pairedcov$graftype22,
  pairedcov$tbidosegp1,pairedcov$tbidosegp2,pairedcov$tbidosegp3,pairedcov$tbidosegp4,
  pairedcov$gvhdprop21,pairedcov$gvhdprop22,pairedcov$gvhdprop23,pairedcov$gvhdprop24,
  pairedcov$gvhdprop27,pairedcov$gvhdprop28,
  pairedcov$Pcluster1,pairedcov$Pcluster2,pairedcov$Pcluster3,pairedcov$Pcluster4,
  pairedcov$Pcluster5,pairedcov$Pcluster6,pairedcov$Pcluster7,
  Wmatrx_paired)
group3 = c(rep(1, times=5), rep(2, times=4), rep(3, times=3), rep(4, times=2),
           rep(5, times=4), rep(6, times=2), rep(7, times=4), rep(8, times=6),
           rep(9, times=7),
           rep(seq(from=(9+1), to=(9+185), by=1), times=2) )
penaltyfactor3 = c(rep(0, times=9), rep(1, times=185))

Pairedinx_3 = GroupLASSOCox(COVARIATE3,intxrel,FullInfo_paired$trm,
                            penaltyfactor3,group3,TRUE,trmMinor)


# fourth outcome
COVARIATE4 = cbind.data.frame(
  pairedcov$indxtxgp1,pairedcov$indxtxgp2,pairedcov$indxtxgp3,pairedcov$indxtxgp4,
  pairedcov$disease10,pairedcov$disease20,pairedcov$disease30,pairedcov$disease40,
  pairedcov$disease50,pairedcov$disease100,
  pairedcov$condint1,pairedcov$condint2,pairedcov$condint3,
  pairedcov$cmldxtx1,pairedcov$cmldxtx2,
  Wmatrx_paired)
group4 = c(rep(1, times=4), rep(2, times=6), rep(3, times=3), rep(4, times=2),
           rep(seq(from=(4+1), to=(4+185), by=1), times=2) )
penaltyfactor4 = c(rep(0, times=4), rep(1, times=185))

Pairedinx_4 = GroupLASSOCox(COVARIATE4,intxrel,FullInfo_paired$relapse,
                            penaltyfactor4,group4,TRUE,relapseMinor)


# fifth outcome
COVARIATE5 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$tbi0,pairedcov$tbi1,
  pairedcov$cmldxtx1,pairedcov$cmldxtx2,
  pairedcov$graftype21,pairedcov$graftype22,
  pairedcov$yeartx20,pairedcov$yeartx21,pairedcov$yeartx22,pairedcov$yeartx23,
  pairedcov$abomatch0,pairedcov$abomatch1,pairedcov$abomatch2,pairedcov$abomatch3,
  pairedcov$firsttx1,
  Wmatrx_paired)
group5 = c(rep(1, times=5), rep(2, times=4), rep(3, times=2), rep(4, times=2),
           rep(5, times=2), rep(6, times=4), rep(7, times=4), 8,
           rep(seq(from=(8+1), to=(8+185), by=1), times=2) )
penaltyfactor5 = c(rep(0, times=8), rep(1, times=185))

Pairedinx_5 = GroupLASSOCox(COVARIATE5,dytxagvh2,FullInfo_paired$agvhd2,
                            penaltyfactor5,group5,TRUE,agvhd2Minor)


# sixth outcome
COVARIATE6 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$tbi0,pairedcov$tbi1,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$kps0,pairedcov$kps1,
  pairedcov$graftype21,pairedcov$graftype22,
  pairedcov$yeartx20,pairedcov$yeartx21,pairedcov$yeartx22,pairedcov$yeartx23,
  pairedcov$gvhdprop21,pairedcov$gvhdprop22,pairedcov$gvhdprop23,pairedcov$gvhdprop24,
  pairedcov$gvhdprop27,pairedcov$gvhdprop28,
  pairedcov$abotyped1,pairedcov$abotyped2,pairedcov$abotyped3,
  pairedcov$Pcluster1,pairedcov$Pcluster2,pairedcov$Pcluster3,pairedcov$Pcluster4,
  pairedcov$Pcluster5,pairedcov$Pcluster6,pairedcov$Pcluster7,
  Wmatrx_paired)
group6 = c(rep(1, times=5), rep(2, times=2), rep(3, times=4), rep(4, times=2),
           rep(5, times=2), rep(6, times=4), rep(7, times=6), rep(8, times=3),
           rep(9, times=7),
           rep(seq(from=(9+1), to=(9+185), by=1), times=2) )
penaltyfactor6 = c(rep(0, times=9), rep(1, times=185))

Pairedinx_6 = GroupLASSOCox(COVARIATE6,dytxagvh3,FullInfo_paired$agvhd3,
                            penaltyfactor6,group6,TRUE,agvhd3Minor)


# seventh outcome
COVARIATE7 = cbind.data.frame(
  pairedcov$yeartx20,pairedcov$yeartx21,pairedcov$yeartx22,pairedcov$yeartx23,
  pairedcov$gvhdprop21,pairedcov$gvhdprop22,pairedcov$gvhdprop23,pairedcov$gvhdprop24,
  pairedcov$gvhdprop27,pairedcov$gvhdprop28,
  pairedcov$graftype21,pairedcov$graftype22,
  pairedcov$dnrsex1,
  pairedcov$condint1,pairedcov$condint2,pairedcov$condint3,
  pairedcov$abotyped1,pairedcov$abotyped2,pairedcov$abotyped3,
  Wmatrx_paired)
group7 = c(rep(1, times=4), rep(2, times=6), rep(3, times=2), 4,
           rep(5, times=3), rep(6, times=3),
           rep(seq(from=(6+1), to=(6+185), by=1), times=2) )
penaltyfactor7 = c(rep(0, times=6), rep(1, times=185))

Pairedinx_7 = GroupLASSOCox(COVARIATE7,intxcgvh,FullInfo_paired$cgvhd,
                            penaltyfactor7,group7,TRUE,cgvhdMinor)


# eighth outcome
COVARIATE8 = cbind.data.frame(
  pairedcov$yeartx20,pairedcov$yeartx21,pairedcov$yeartx22,pairedcov$yeartx23,
  pairedcov$disease10,pairedcov$disease20,pairedcov$disease30,pairedcov$disease40,
  pairedcov$disease50,pairedcov$disease100,
  pairedcov$tbi0,pairedcov$tbi1,
  pairedcov$dnrsex1,
  pairedcov$disstat1,pairedcov$disstat2,pairedcov$disstat3,
  pairedcov$abomatch0,pairedcov$abomatch1,pairedcov$abomatch2,pairedcov$abomatch3,
  pairedcov$kps0,pairedcov$kps1,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$htpr21,pairedcov$htpr22,pairedcov$htpr23,pairedcov$htpr24,
  Wmatrx_paired)
group8 = c(rep(1, times=4), rep(2, times=6), rep(3, times=2), 4,
           rep(5, times=3), rep(6, times=4), rep(7, times=2), rep(8, times=4),
           rep(9, times=4),
           rep(seq(from=(9+1), to=(9+185), by=1), times=2) )
penaltyfactor8 = c(rep(0, times=9), rep(1, times=185))

Pairedinx_8 = GroupLASSOCox(COVARIATE8,dytxanc,FullInfo_paired$anc500,
                            penaltyfactor8,group8,TRUE,anc500Minor)


# ninth outcome
COVARIATE9 = cbind.data.frame(
  pairedcov$hla_match5,pairedcov$hla_match6,pairedcov$hla_match7,pairedcov$hla_match8,
  pairedcov$hla_match9,
  pairedcov$yeartx20,pairedcov$yeartx21,pairedcov$yeartx22,pairedcov$yeartx23,
  pairedcov$disease10,pairedcov$disease20,pairedcov$disease30,pairedcov$disease40,
  pairedcov$disease50,pairedcov$disease100,
  pairedcov$cmvmatch1,pairedcov$cmvmatch2,pairedcov$cmvmatch3,pairedcov$cmvmatch4,
  pairedcov$kps0,pairedcov$kps1,
  pairedcov$gvhdprop21,pairedcov$gvhdprop22,pairedcov$gvhdprop23,pairedcov$gvhdprop24,
  pairedcov$gvhdprop27,pairedcov$gvhdprop28,
  pairedcov$disstat1,pairedcov$disstat2,pairedcov$disstat3,
  pairedcov$abomatch0,pairedcov$abomatch1,pairedcov$abomatch2,pairedcov$abomatch3,
  pairedcov$ragecat21,pairedcov$ragecat22,
  pairedcov$wtpr21,pairedcov$wtpr22,pairedcov$wtpr23,pairedcov$wtpr24,
  pairedcov$htpr21,pairedcov$htpr22,pairedcov$htpr23,pairedcov$htpr24,
  Wmatrx_paired)
group9 = c(rep(1, times=5), rep(2, times=4), rep(3, times=6), rep(4, times=4),
           rep(5, times=2), rep(6, times=6), rep(7, times=3), rep(8, times=4),
           rep(9, times=2), rep(10, times=4), rep(11, times=4), 
           rep(seq(from=(11+1), to=(11+185), by=1), times=2) )
penaltyfactor9 = c(rep(0, times=11), rep(1, times=185))

Pairedinx_9 = GroupLASSOCox(COVARIATE9,dytxpl20,FullInfo_paired$plate20,
                            penaltyfactor9, group9,TRUE,plate20Minor)


detach(FullInfo_paired)


