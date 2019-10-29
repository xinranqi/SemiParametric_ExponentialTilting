
#####################################################################################################################################
### Conditional Expectation Maximization Maximize Minorize + Numerical Methods for generalized marginal compatibility (GMC) model ###
#####################################################################################################################################
# GMC form of univariate marginal probability for baseline dose level with cluster size r, i.e., q_y,r_(0)
# what is in the numerator of the GMC form of q_y,r_(0)
Q_yr_0_numer = function(yResp, rCluSize, qMaxClus, Rmax, wTilt) {
  ts = seq(from=yResp, to=min(Rmax-rCluSize+yResp, Rmax), by=1)
  hypergeoms = dhyper(yResp, ts, Rmax-ts, rCluSize, log=FALSE)
  sum(hypergeoms * qMaxClus[ts+1] * exp(wTilt*(yResp-ts)))
}



# return a vector of GMC form of q_y,r_(0) with corresponding support yResp between 0 ~ rClusterSize
Q_r_0 = function(rCluSize, qMaxClus, Rmax, wTilt) {
  Qr0unscaled = numeric(rCluSize + 1)
  
  for (Yresps in 0:rCluSize) {
    Qr0unscaled[Yresps+1] = Q_yr_0_numer(Yresps, rCluSize, qMaxClus, Rmax, wTilt)
  }
  
  Qr0 = if(sum(Qr0unscaled)==0) {
    rep(0, times=rCluSize+1)
  } else {
    Qr0unscaled/sum(Qr0unscaled)
  }
  
  Qr0
}
# to be honest, once we have changed the EM algorithm estimation path, the above two functions are not applied in neither E nor M steps



# Conditional Expectation Maximization Maximize Minorize + Numerical Methods Algorithm
# expectation step with respect to the observed layer
# a vector of univariate marginal probabilities for clusters of size Rmax in the latent layer
Q_R_latent = function(Rmax, wTilt, qMaxClus) {
  tprime_ = seq(from=0, to=Rmax, by=1)
  q_R_tilted = qMaxClus * exp(-wTilt*tprime_)
  # normalizing constant for the maximum cluster size unscaled exponential tilted univariate marginal probabilities
  Normalizer_ = sum(q_R_tilted)
  
  # diagnosis of normalizing constant ==0 or not
  q_R_tilted = if(Normalizer_==0) {
    rep(0, times=Rmax+1)
  } else {
    q_R_tilted / Normalizer_
  }
  
  list(Normalizer_ = Normalizer_,
       q_R_tilted = q_R_tilted)
}



# function to calculate conditional probability p_r,s,t_(0) in the latent layer: for the baseline dose level, given that s responses out of cluster size r, probability of having t responses out of maximum cluster size R
P_rs_tLatent = function(Rmax, wTilt, qMaxClus) {
  # an array of conditional probability p_r,s,t_(0) in the latent layer
  prst.array = array(0, dim=c(Rmax+1, Rmax, Rmax+1), dimnames=list(paste0("ResponsesFromr", 0:Rmax),
                                                                   paste0("ClusterSize", 1:Rmax),
                                                                   paste0("ResponsesFromR", 0:Rmax)))
  qR_latent_ = Q_R_latent(Rmax, wTilt, qMaxClus)$q_R_tilted
  
  # fulfillment of entries in prst.array, unscaled for now
  for (rs in 1:Rmax) {
    for (Ss in 0:rs) {
      Ts = seq(from=Ss, to=min(Rmax, Rmax-rs+Ss), by=1)
      # conditional probability Pr(S=s|T=t): given that t responses out of maximum cluster size R, probability of having s responses out of cluster size r
      # Pr(S=s|T=t, latent) calculated from hypergeometric distribution since the latent layer is subject to marginal compatibility
      prst.array[Ss+1, rs, Ts+1] = dhyper(Ss, Ts, Rmax-Ts, rs, log=FALSE) * qR_latent_[Ts+1]
    }
  }
  
  # a matrix of normalizing constants to normalize each vector of conditonal probabilities p_r,s,t_(0) corresponding to (rs, Ss), i.e., (cluster size, # of responses)
  prst.SumOverT = apply(prst.array, 1:2, sum)
  
  # normalize to pmf scale
  for (rs in 1:Rmax) {
    for (Ss in 0:rs) {
      Ts = seq(from=Ss, to=min(Rmax, Rmax-rs+Ss), by=1)
      
      prst.array[Ss+1, rs, Ts+1] = if(prst.SumOverT[Ss+1, rs]==0) {
        rep(0, times=min(Rmax, Rmax-rs+Ss)-Ss+1)
      } else {
        prst.array[Ss+1, rs, Ts+1]/prst.SumOverT[Ss+1, rs]
      }
      
    }
  }
  
  # correction of trivial situation: for maximum cluster size, when t==s, conditional probability p_r,s,t_(0)==1 because given that s responses out of maximum cluster size R, with probability 1 there are t responses out of maximum cluster size R
  diag(prst.array[, Rmax, ]) = rep(1, times=Rmax+1)
  prst.array[, Rmax, ][lower.tri(prst.array[, Rmax, ])] = rep(0, times=Rmax*(Rmax+1)/2)
  prst.array[, Rmax, ][upper.tri(prst.array[, Rmax, ])] = rep(0, times=Rmax*(Rmax+1)/2)
  
  prst.array
}



# matrix of nrs in the observed layer: number of observed clusters with s reponses out of cluster size r
# few words about input Data structure: two columns, i.e., nTotal & nResp
N_rs = function(Data, Rmax) {
  nrs.array = array(0, dim=c(Rmax, Rmax+1, dim(Data)[1]), dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                        paste0("ResponsesFromr", 0:Rmax),
                                                                        paste0("Observation", 1:dim(Data)[1])))
  nrs.matrix = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                               paste0("ResponsesFromr", 0:Rmax)))
  for (rs in 1:Rmax) {
    for (Ss in 0:rs) {
      for (obs in 1:dim(Data)[1]) {
        nrs.array[rs, Ss+1, obs] = Data$nTotal[obs]==rs & Data$nResp[obs]==Ss
      }
      nrs.matrix[rs, Ss+1] = sum(nrs.array[rs, Ss+1, ])
    }
  }
  
  nrs.matrix
}



# a vector of univariate marginal probabilities out of cluster of size r (<=R) in the latent layer, calculated from the marginal compatibility assumption (satisfied in the latent layer)
Q_r_latent = function(Rmax, wTilt, qMaxClus, rClusterSize) {
  qR_latent_ = Q_R_latent(Rmax, wTilt, qMaxClus)$q_R_tilted
  qr_latent_ = numeric(rClusterSize+1)
  
  # fulfillment of entries in the qr_latent_ vector
  for(Ss in 0:rClusterSize) {
    Ts = seq(from=Ss, to=min(Rmax, Rmax-rClusterSize+Ss), by=1)
    qr_latent_[Ss+1] = sum(dhyper(Ss, Ts, Rmax-Ts, rClusterSize, log=FALSE) * qR_latent_[Ts+1])
  }
  
  qr_latent_
}



# matrix of expectation of nrs in the latent layer, based on the observed nrs
N_rsLatent = function(Data, Rmax, wTilt, qMaxClus) {
  # nrs maxtrix in the observed layer
  Nrs_ = N_rs(Data, Rmax)
  # expectation of nrs matrix for the latent layer
  nrsLatent.matrix = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                     paste0("ResponsesFromr", 0:Rmax)))
  for(rs in 1:Rmax) {
    nrsLatent.matrix[rs, seq(from=1, to=(rs+1), by=1)] = sum(Nrs_[rs, ]) * Q_r_latent(Rmax, wTilt, qMaxClus, rs)
  }
  
  nrsLatent.matrix
}



# expectation of N_t in the latent layer: under complete data (cluster size being R), number of clusters having t responses
N_t_primeLatent = function(Data, Rmax, wTilt, qMaxClus) {
  NrsLatent_ = N_rsLatent(Data, Rmax, wTilt, qMaxClus)
  PrstLatent_ = P_rs_tLatent(Rmax, wTilt, qMaxClus)
  # a vector of expectation of N_t in the latent layer
  N_tPrimeLatent = numeric(Rmax+1)
  
  # fulfillment of entries in the vector of expectation of N_t in the latent layer
  for (Ts in 0:Rmax) {
    # array PrstLatent_ is arranged as [Ss, rs, Ts], i.e., [# responses from r, cluster size r, # responses from R], thus transpose is need in matrix elemnentwise multiplication
    N_tPrimeLatent[Ts+1] = sum(NrsLatent_ * t(PrstLatent_[, , Ts+1]))
  }
  
  N_tPrimeLatent
}



# expectation of N_t in the observed layer: under complete data (cluster size being R), number of clusters having t responses; calculated via the linkage equation between the observed and latent layers
N_t_prime = function(Data, Rmax, wTilt, qMaxClus) {
  NormalizeC_ = Q_R_latent(Rmax, wTilt, qMaxClus)$Normalizer_
  NtLatent_ = N_t_primeLatent(Data, Rmax, wTilt, qMaxClus)
  # expectation of N_t in the observed layer
  NtObserved_ = numeric(Rmax+1)
  
  # fulfillment of expectation of N_t in the observed layer
  for (Ts in 0:Rmax) {
    NtObserved_[Ts+1] = NormalizeC_ * NtLatent_[Ts+1] * exp(wTilt*Ts)
  }
  
  NtObserved_
}



# maximization step
Q_t_maximize = function(Data, Rmax, wTilt, qMaxClus) {
  N_tprime_ = N_t_prime(Data, Rmax, wTilt, qMaxClus)
  
  qt_ = if(sum(N_tprime_)==0) {
    rep(0, times=Rmax+1)
  } else {
    N_tprime_/sum(N_tprime_)
  }

  qt_
}



# Newton-Raphson algorithm to updata parameter wTilt
# return a vector of unscaled GMC form q_y,r_(0) with support 0 ~ rClusterSize
MomentsHelper = function(yResp, rCluSize, qMaxClus, Rmax, wTilt) {
  ts = seq(from=yResp, to=min(Rmax-rCluSize+yResp, Rmax), by=1)
  hypergeoms = dhyper(yResp, ts, Rmax-ts, rCluSize, log=FALSE)
  hypergeoms * qMaxClus[ts+1] * exp(wTilt*(yResp-ts))
}



# moment sums of unscaled GMC form q_y,r_(0)
Moments = function(yResp, rCluSize, qMaxClus, Rmax, wTilt) {
  ts = seq(from=yResp, to=min(Rmax-rCluSize+yResp, Rmax), by=1)
  pmf = MomentsHelper(yResp, rCluSize, qMaxClus, Rmax, wTilt)
  firstM = pmf * (yResp - ts)
  secondM = pmf * ((yResp - ts)^2)
  
  list(pmfsum = sum(pmf),
       firstMsum = sum(firstM),
       secondMsum = sum(secondM))
}



# function to calculate first derivative + quadratic (w^2) penality included
FirsSecOrder_derivatives = function(Data, Rmax, wTilt, qMaxClus) {
  Nrs_ = N_rs(Data, Rmax)
  
  # matrix of moment sums stored in advance for quick sort
  pmfsumMatrix = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                 paste0("ResponsesFromr", 0:Rmax)))
  firstMsumMatrix = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                    paste0("ResponsesFromr", 0:Rmax)))
  secondMsumMatrix = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                     paste0("ResponsesFromr", 0:Rmax)))
  for (rs in 1:Rmax) {
    for (Ss in 0:rs) {
      Sums_ = Moments(Ss, rs, qMaxClus, Rmax, wTilt)
      pmfsumMatrix[rs, Ss+1] = Sums_$pmfsum
      firstMsumMatrix[rs, Ss+1] = Sums_$firstMsum
      secondMsumMatrix[rs, Ss+1] = Sums_$secondMsum
    }
  }
  
  # normalizing constant: sum of moment sums according to each row <-> 0:r; normalizing constants is invariant to # of responses, but depend on the cluster size. Thus the function returns a vector of normalizing constants w.r.t different cluster sizes (1:Rmax)
  NormalizePmf = apply(pmfsumMatrix, 1, sum)
  NormalizeFirstM = apply(firstMsumMatrix, 1, sum)
  NormalizeSecondM = apply(secondMsumMatrix, 1, sum)
  
  # matrix of expectation form, prepared for multiplication with Nrs_
  # comp1: first order derivative components
  # comp2: second order derivative components
  comp1 = matrix(0, nrow=Rmax, ncol=Rmax+1)
  comp2 = matrix(0, nrow=Rmax, ncol=Rmax+1)
  
  for (rs in 1:Rmax) {
    for (Ss in 0:rs) {
      comp1[rs, Ss+1] = ifelse(pmfsumMatrix[rs, Ss+1]==0, 0, firstMsumMatrix[rs, Ss+1]/pmfsumMatrix[rs, Ss+1]) - ifelse(NormalizePmf[rs]==0, 0, NormalizeFirstM[rs]/NormalizePmf[rs])
      comp2[rs, Ss+1] = ifelse(pmfsumMatrix[rs, Ss+1]==0, 0, (secondMsumMatrix[rs, Ss+1]*pmfsumMatrix[rs, Ss+1] - (firstMsumMatrix[rs, Ss+1])^2) / ((pmfsumMatrix[rs, Ss+1])^2)) - ifelse(NormalizePmf[rs]==0, 0, (NormalizeSecondM[rs]*NormalizePmf[rs] - (NormalizeFirstM[rs]^2)) / (NormalizePmf[rs]^2))  
    }
  }
  
  list(NRS_ = Nrs_,
       comp1 = comp1,
       comp2 = comp2)
}



# use R package nloptr to solve for wTilt, i.e., the exponential tilting parameter
library('nloptr')
  
nlSolvewTilt = function(Data, Rmax, wTiltInitial, qMaxClus, penalty) {
  # plug in qMaxClus from previous updates and other parameters to get the loglikelihood function with uni-parameter wTilt
  loglikelihoods_ = function(WTilt) {
    LogLikelihood(Data, Rmax, WTilt, qMaxClus, penalty)
  }
  
  # partial derivative w.r.t wTilt after plugging in qMaxClus from previous updates and other parameters
  partialwTilt_ = function(WTilt) {
    componentMatrices = FirsSecOrder_derivatives(Data, Rmax, WTilt, qMaxClus)
    sum((componentMatrices$NRS_) * (componentMatrices$comp1)) - 2*penalty*WTilt
  }
  
  # objective function to be minimized, nonlinear problem
  eval_f = function(WTilting) {
    return(list("objective" = -loglikelihoods_(WTilting),
                "gradient"  = -partialwTilt_(WTilting)))
  }
  
  local_opts = list("algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1e-4)
  opts = list("algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1e-4,
              "maxeval" = 1000,
              "local_opts" = local_opts)
  
  res = nloptr(x0 = wTiltInitial,
               eval_f = eval_f,
               opts = opts)
  res$x0
}



# use R function optim to solve for wTilt, i.e., the exponential tilting parameter
OptimSolvewTilt = function(Data, Rmax, wTiltInitial, qMaxClus, penalty) {
  # objective function to be minimized, nonlinear programming problem
  # plug in qMaxClus from previous updates and other parameters to get the loglikelihood function with uni-parameter wTilt
  minusloglikelihoods_ = function(WTilt) {
    -LogLikelihood(Data, Rmax, WTilt, qMaxClus, penalty)
  }
  
  # partial derivative w.r.t wTilt after plugging in qMaxClus from previous updates and other parameters
  minuspartialwTilt_ = function(WTilt) {
    componentMatrices = FirsSecOrder_derivatives(Data, Rmax, WTilt, qMaxClus)
    -sum((componentMatrices$NRS_) * (componentMatrices$comp1)) + 2*penalty*WTilt
  }

  # 1-D minimization: "Brent" or optimize() being preferred
  res = optim(wTiltInitial, minusloglikelihoods_, minuspartialwTilt_, method = "Brent", lower = -50, upper = 50)
  res$par
  
  # using optimize function
  # res = optimize(minusloglikelihoods_, interval=c(-50, 50), maximum = FALSE, tol=1e-4)
  # res$minimum
  
  # using nlm function
  # res = nlm(minusloglikelihoods_, wTiltInitial)
  # res$estimate
  
  # using uniroot function
  # res = uniroot(minuspartialwTilt_, interval=c(-50, 50), tol=1e-4)
  # res$root
}


 
# simulations have shown that using NR to solve for wTilt leads to unstable estimates, thus ruled out
# Newtown Raphson algorithm to update wTilt, i.e., exponential tilting parameter
NRsolvewTilt = function(Data, Rmax, wTilt, qMaxClus, penalty, itermax=100, tol=1e-4) {
  iter = 0
  difference = 10
  
  while (iter<itermax & difference>tol) {
    iter = iter + 1
    wTiltPre = wTilt
    qMaxClusPre = qMaxClus
    
    compMatrices = FirsSecOrder_derivatives(Data, Rmax, wTiltPre, qMaxClusPre)
    denomStepsize =  sum((compMatrices$NRS_) * (compMatrices$comp2)) - 2*penalty
    
    stepsize = if(denomStepsize==0) {
      0
    } else {
      (sum((compMatrices$NRS_) * (compMatrices$comp1)) - 2*penalty*wTiltPre) / denomStepsize
    }
    
    wTilt = wTiltPre - stepsize
    difference = abs(stepsize)
  }
  
  wTilt
}



# estimation of parameters involved in the GMC model via conditional expectation maximization maximize minorize + Newton Raphson algorithm
GMC.est = function(Data, Rmax, wTiltInit, qMaxClusInit, penalty, itermax=500, tol=1e-4) {
  Iter = 0
  delta = 10
  wTiltingParameter = wTiltInit
  qCompleteData = qMaxClusInit
  cat("log likelihood before parameter update is:", LogLikelihood(Data, Rmax, wTiltingParameter, qCompleteData, penalty), "\n")
  
  while (Iter<itermax & delta>tol) {
    Iter = Iter + 1
    wTiltingParameterPre = wTiltingParameter
    qCompleteDataPre = qCompleteData
    
    # wTiltingParameter = nlSolvewTilt(Data, Rmax, wTiltingParameterPre, qCompleteDataPre, penalty)
    wTiltingParameter = OptimSolvewTilt(Data, Rmax, wTiltingParameterPre, qCompleteDataPre, penalty)
    # wTiltingParameter = NRsolvewTilt(Data, Rmax, wTiltingParameterPre, qCompleteDataPre, penalty)
    cat("In iteration",Iter,", the log-likelihood after wTilt but before qMaxClus update is:", LogLikelihood(Data, Rmax, wTiltingParameter, qCompleteDataPre, penalty), "\n")
    
    # use the updated w estimate to update the univariate marginal probabilities for the maximum cluster size
    qCompleteData = Q_t_maximize(Data, Rmax, wTiltingParameter, qCompleteDataPre)
    cat("In iteration",Iter,", the log-likelihood after both wTilt and qMaxClus updates is:", LogLikelihood(Data, Rmax, wTiltingParameter, qCompleteData, penalty), "\n")
    
    delta = abs(wTiltingParameter-wTiltingParameterPre) + sum(abs(qCompleteData-qCompleteDataPre))
  }
  
  list(qCompleteData = qCompleteData,
       wTiltingParameter = wTiltingParameter,
       Iter = Iter)
}



# log-likelihood function, used for likelihood ratio test & non-linear optimization programming; quadratic (w^2) penalty term included
LogLikelihood = function(Data, Rmax, wTilt, qMaxClus, penalty) {
  Nrs_ = N_rs(Data, Rmax)
  
  # matrix to store univariate marginal probabilities for varying cluster sizes
  q_r_s_ = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                           paste0("ResponsesFromr", 0:Rmax)))
  for(rs in 1:Rmax) {
    Ss = seq(from=0, to=rs, by=1)
    q_r_s_[rs, Ss+1] = Q_r_0(rs, qMaxClus, Rmax, wTilt)
  }
  
  # matrix to store elementwise product of nrs * log(q_s,r)
  elementProduct = matrix(0, nrow=Rmax, ncol=Rmax+1, dimnames=list(paste0("ClusterSize", 1:Rmax),
                                                                   paste0("ResponsesFromr", 0:Rmax)))
  for (rs in 1:Rmax) {
    for (ss in 0:rs) {
      elementProduct[rs, ss+1] = ifelse(Nrs_[rs, ss+1]>0 & q_r_s_[rs, ss+1]>0, Nrs_[rs, ss+1] * log(q_r_s_[rs, ss+1]), 0)
    }
  }
  
  sum(elementProduct) - penalty*(wTilt^2)
}






# data generation: given that the maximum cluster size being 10, simulation as follows
set.seed(123)
Rglobal = 10
w = -5

# matrix to store univariate marginal probabilities for varying cluster sizes, rows <-> cluster size; columns <-> # of responses
QvaryingClusSizes = matrix(0, nrow=10, ncol=Rglobal+1)

# univariate marginal probabilities for the maximum cluster size
# QvaryingClusSizes[10, ] = rep(1/(Rglobal+1), times=Rglobal+1)
QvaryingClusSizes[10, ] = dbinom(0:Rglobal, size=Rglobal, prob=0.5, log=FALSE)

for (ClusterSize in (Rglobal-9):(Rglobal-1)) {
  QvaryingClusSizes[ClusterSize, 1:(ClusterSize+1)] = Q_r_0(ClusterSize, QvaryingClusSizes[10, ], Rglobal, w)
}

# simulation dataset
RandomDraws = matrix(0, nrow=100, ncol=2)
# we assume that the distribution of clustersize is uniform between 3 ~ 5
RandomDraws[, 1] = sample(seq(from=Rglobal-9, to=Rglobal, by=1), size=nrow(RandomDraws), replace=TRUE, prob=rep(0.1, times=10))
for (observ in 1:nrow(RandomDraws)) {
  clustSize = RandomDraws[observ, 1]
  RandomDraws[observ, 2] = sample(seq(from=0, to=clustSize, by=1), size=1, prob=QvaryingClusSizes[clustSize, 1:(clustSize+1)])
}
colnames(RandomDraws) = c("nTotal", "nResp")
RandomDraws = as.data.frame(RandomDraws)



# initial values similar to corresponding true values
# the first output corresponds to w=0.5; whereas the second output corresponds to w=1.5
wTiltInits = 0
# qMaxClusInits = dbinom(0:Rglobal, size=Rglobal, prob=0.5, log=FALSE)
qMaxClusInits = rep(1/(Rglobal+1), times=Rglobal+1)
ests = GMC.est(RandomDraws, Rglobal, wTiltInits, qMaxClusInits, penalty = 0.05)
ests



ws = seq(from=-5, to=5, by=0.05)
loglikehoods1 = numeric(length(ws))
# loglikehoods2 = numeric(length(ws))
for (i in 1:length(ws)) {
  loglikehoods1[i] = LogLikelihood(RandomDraws, Rglobal, ws[i], QvaryingClusSizes[10, ], penalty = 0.05)
  # loglikehoods2[i] = LogLikelihood(RandomDraws, Rglobal, ws[i], ests$qCompleteData)
}
plot(x=ws, y=loglikehoods1)
# plot(x=ws, y=loglikehoods2)
ws[which.max(loglikehoods1)]
# ws[which.max(loglikehoods2)]



