# This is the code to create the plot for the confidence intervals.

suppressMessages({
 
  if (!require(writexl)) install.packages('writexl', repos =
                                            'https://cran.revolutionanalytics.com/')
  
  
  library(writexl)
  
})

setwd("/Users/daiguorong/Dropbox/CCSHD/Uploaded_Files")

source("~/Dropbox/CCSHD/Simulation/CCSHD_Functions.R")

dat = readRDS(file = "sepsis_data.rds")

D = dat $ D
S = dat $ S
X = dat $ X
Z = dat $ Z

n = length(D)
N1 = sum(S)
N = length(S)
N0 = N - N1

analysis = "primary"
method = c("LR", "KS", "RF")
M = 10
Y = NULL

tau = N1 / N
eta = NULL

output = matrix( 0, ncol(X) - 1, 2 * (length(method) + 1) )

initial.fit = prelim(tau, eta, D, Y, S, X, analysis)
thetahatIN = initial.fit $ initial.est
Omegahat = initial.fit $ Jacobian
psi = initial.fit $ psi
cov.thetahatIN = initial.fit $ initial.cov

std.initial = sqrt(diag(cov.thetahatIN))
CI.initial = cbind(thetahatIN - 1.96 * std.initial, thetahatIN + 1.96 * std.initial) 
# 95% confidence intervals of the coefficients based on thetahatIN, recorded in a p * 2 matrix
CIL.initial = CI.initial[, 2] - CI.initial[, 1]

output[ , 1 : 2] = CI.initial[-1, ]

Z.V = Z[1 : n, ]
Z.N = Z[(n + 1) : N1, ]

for(j in 1 : length(method))
{
  
  set.seed(1)
  muhatZ = muhat.function(D, Z.V, Z.N, M, method[j])
  
  result = CCS.M.Est(tau, eta, D, Y, S, X, muhatZ, thetahatIN, Omegahat, psi, analysis)
  thetahat = result $ est
  cov.thetahat.new = result $ cov.new
  # cov.thetahat = result $ cov
  
  std.new = sqrt(diag(cov.thetahat.new))
  CI.new = cbind(thetahat - 1.96 * std.new, thetahat + 1.96 * std.new) 
  # 95% confidence intervals of the coefficients based on thetahat, recorded in a p * 2 matrix
  
  output[ , 2 * j + (1 : 2)] = CI.new[-1, ]
  
}

dat = output

X.names = c("urineoutput", "lactate_min", "bun_mean", "sysbp_min", "age", "sodium_max",
            "aniongap_max", "creatinine_min", "spo2_mean", "inr_max", "metastatic_cancer")
uu = length(X.names)

lb = numeric(0)
ub = numeric(0)
for(i in 1 : uu)
{
  
  dat1 = dat[i, ]
  lb = c(lb, dat1[seq(1, ncol(dat) - 1, 2)])
  ub = c(ub, dat1[seq(2, ncol(dat), 2)])
  
}
covariates = rep(factor(X.names), each = 4)
# method = factor(1 : 4)
method = c("Zero", "LR", "KS", "RF")
method = factor(method, levels = c("Zero", "LR", "KS", "RF"))

ci = data.frame("Covariates" = covariates, "Imputation" = method, "Coefficients" = (lb + ub) / 2,
                "lb" = lb, "ub" = ub)

g1 = ggplot(data = ci, aes(x = Covariates, y = Coefficients, ymin = lb, 
                           ymax = ub, fill = Imputation))

g2 = geom_crossbar_pattern(stat = "identity", 
                           pattern = "none",
                           pattern_density = 0.2,
                           pattern_spacing = 0.007,
                           pattern_fill = 'black',
                           position = position_dodge(width = 0.62),  
                           width = 0.62, show.legend = T)
g1 + g2 + theme(
  legend.position = c(.95, .95),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
) 
