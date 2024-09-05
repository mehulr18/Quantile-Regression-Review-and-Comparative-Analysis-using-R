rm(list=ls()) # clear environment

# !! Packages needed (you need to comment-them-in if you don't have these installed already) !!
# install.packages("quantreg")
# install.packages("truncnorm")
# install.packages("boot")
# install.packages("tidyverse")

# Activation of the required packages
library(truncnorm)
library(quantreg)
library(boot)
library(tidyverse)
library(ggplot2)

set.seed(4321) #setting a seed to ensure replicability

# Data Generating Process (DGP) - Proposition:
# Using the model provided in the paper "Revisiting the income, energy consumption and carbon emission nexus: new 
# evidence from quantile regression for different country groups" Yaya Keho (2017)
# CO2 = beta0 + beta1 * income + beta2 * income^2 + epsilon
# - CO2: emissions resp. environmental degradation
# - income: GDP per capita
# - income^2: GDP per capita squared
# - epsilon: error term (heteroscedastic)
# - beta0: intercept
# - beta1: is expected to be positive
# - beta2: is expected to be negative (environmental Kuznets curve)

# Create a Data Generating Process (DGP) function, that simulates data 

# Proposition: A HETEROSCEDASTIC DGP
dgp_heteroscedastic <- function(n,beta0,beta1,beta2){
  x <- rtruncnorm(n, a=0, b=Inf, mean = 100, sd =30) # independent variable (Income)
  eps <- rtruncnorm(n, a=-1.5*beta0, b= Inf, mean =0, sd = x*70) # heteroscedastic error term
  y <- beta0 + beta1 * x + beta2 * x^2 + eps # the dgp (the true relationship of the regressor and the regressant; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}

################################################################################
######### GRAPHICAL VISUALIZATION OF THE DATA AND THE REGRESSION LINES #########
################################################################################

# Generating the data according to our dgp (here: heteroscedastic data)
model_coefficients <- c(5000,175,-0.75)
new_data <- dgp_heteroscedastic(n = 100000,
                 beta0 = 5000,
                 beta1 = 175,
                 beta2 = -0.75)
# Defining our quantiles of interest 
qs <- c(0.05,0.25,0.5,0.75,0.95) 

# Plotting the data as well as a linear regression line and the respective 
#     quantile regression lines
plottings <- ggplot(data = new_data,
                    mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Relationship between Carbon Emissions and Income (Environmental Kuznets Curve)") +
  ylab("Carbon Emissions (in tons)") +
  xlab("Income (scaled)") +
  geom_smooth(data = new_data,
              aes(color= "OLS Regression"),
              linetype = 1,
              method = "lm",
              formula = y ~ poly(x,2, raw=TRUE)) +
  geom_quantile(data = new_data,
                aes(color = "Quantile Regression"),
                linetype = 1,
                quantiles = qs, 
                formula = y ~ poly(x,2, raw=TRUE))+
  scale_fill_manual(values = c( "OLS Regression" = "darkgreen",
                                "Quantile_Regression" = "darkred"))+
  labs(color = "Regression Lines")

plottings # you need to execute this line to see the plot

################################################################################
########################## SIMULATION PART #####################################
################################################################################

#------------------------------------------------------------------------------#
#DIGRESSION:
#Now, before we can really start with the simulation part, we need to make a few
# steps.
#This is, because our model only determines the average conditional relationship
# (the beta-coefficients) between the dependent variable and the covariates.
#But, since we are interested in some particular quantiles, we need to "extract"
# the respective "true" values of the beta-coefficients first, in order to use 
# them later to make assessments of quantile regression in general.

#For this, I propose to create a thousand repeated draws using our dgp and then
# taking the average value, we obtain for the respective betas

#Defining a function, that calculates the "true" betas as averages from S repeated
# draws of the dgp with n observations

beta_averages <- function(S){
  beta0s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta0 coefficient estimates
  beta1s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta1 coefficient estimates
  beta2s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta2 coefficient estimates
  for(i in 1:S){ #repeating S times
    df <- dgp_heteroscedastic(n=10000, 
                              beta0=5000,
                              beta1=175,
                              beta2=-0.75) # every iteration we draw a new data set
    beta0s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                    data=df,
                    tau=qs)$coefficients)[,1] # stores the coefficient estimates for beta0 (for all resp. quantiless)
    beta1s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                          data=df,
                          tau=qs)$coefficients)[,2] # stores the coefficient estimates for beta1 (for all resp. quantiles)
    beta2s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                          data=df,
                          tau=qs)$coefficients)[,3] # stores the coefficient estimates for beta2 (for all resp. quantiles)
  }
  beta0_avgs <- rowMeans(beta0s) #calculates the average beta0 for all resp. quantiles
  beta1_avgs <- rowMeans(beta1s) #calculates the average beta0 for all resp. quantiles
  beta2_avgs <- rowMeans(beta2s) #calculates the average beta0 for all resp. quantiles
  quantiles <- qs #including the respective quantiles to make it more appealing
  averages <- t(data.frame(quantiles,beta0_avgs,beta1_avgs,beta2_avgs)) #creating a data frame (and transposing it)
  return(averages)
}
true_betas <- beta_averages(S=1000) # (NOTE: calculation time approximately 7.5 minutes on my laptop when using S=1000)
true_betas
# so this creates a (4x5) matrix such that:
# in the first row we have our quantiles of interest
# row 2,3, and 4 are the respective regression coefficients for these quantiles

#------------------------------------------------------------------------------#

#The following code is for using monte-carlo simulations in order to asses the 
# finite sample properties of quantile regression.
#Therefore, the following will be investigated below:
# 1) mean-squared errors (mse) of the beta estimates 
#     -> measure of the quality of the estimates.
# 2) consistency of the beta estimates
# 3) coverage probability

# MC simulates numerous random samples from an underlying distribution to obtain estimator realizations for each sample
# Then, it uses these realizations to approx the small samples distribution of the estimator and evaluates properties
# like coverage prob. of confidence intervals
# The estimator properties vary based on the chosen underlying distribution

S = 1000 #number of Monte-Carlo replications

## (1) CALCULATING THE MEAN SQUARED ERROR OVER DIFFERENT SAMPLE SIZES FOR OUR ESTIMATES ##

#### Creating a function that calculates the MEAN SQUARED ERROR (mse)
# note: MSE(Theta) = Mean((E(Theta_hat) - Theta)^2)

# This function can later be used to calculate the mse over a sequence of sample sizes using sapply()
mse_function_qr_n_seq <- function(n){
  beta0s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta0 coefficient estimates
  beta1s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta1 coefficient estimates
  beta2s <- matrix(NA,nrow=length(qs),ncol=S) #empty matrix for beta2 coefficient estimates
  mse_beta0 <- matrix(NA,nrow=length(qs),ncol=1) #empty matrix for the mses of beta0
  mse_beta1 <- matrix(NA,nrow=length(qs),ncol=1) #empty matrix for the mses of beta1
  mse_beta2 <- matrix(NA,nrow=length(qs),ncol=1) #empty matrix for the mses of beta2
  betas_ols <- matrix(NA,nrow=3,ncol=S) # empty matrix for the beta coefficients estimates for a ordinary least squares model
  mse_betas_ols <- matrix(NA,nrow=3,ncol=1) # empty matrix for the mses of the beta coefficient estimates for the ols model
  for(i in 1:S){ #repeating S times
    df <- dgp_heteroscedastic(n=n, 
                              beta0=5000,
                              beta1=175,
                              beta2=-0.75) # every iteration we draw a new data set
    fit_ols <- lm(Carbon_Emissions~poly(Income,2,raw=TRUE),
                  data=df)
    betas_ols[,i] <- fit_ols$coefficients
    beta0s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                       data=df,
                       tau=qs)$coefficients)[,1] # stores the coefficient estimates for beta0 (for all resp. quantiles)
    beta1s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                       data=df,
                       tau=qs)$coefficients)[,2] # stores the coefficient estimates for beta1 (for all resp. quantiles)
    beta2s[,i] <- t(rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                       data=df,
                       tau=qs)$coefficients)[,3] # stores the coefficient estimates for beta2 (for all resp. quantiles)
  }
  for(q in 1:3){
    mse_betas_ols[q,1] <- mean((betas_ols[q,]-model_coefficients[q])^2)
  }
  for(z in 1:length(qs)){
    mse_beta0[z,1] <- mean((beta0s - true_betas[2,z])^2) # calculates the mse for the estimate of beta0 for all resp. quantiles
    mse_beta1[z,1] <- mean((beta0s - true_betas[3,z])^2) # calculates the mse for the estimate of beta1 for all resp. quantiles
    mse_beta2[z,1] <- mean((beta0s - true_betas[4,z])^2) # calculates the mse for the estimate of beta2 for all resp. quantiles
    quantiles <- qs
  }
  df3 <- rbind(mse_beta0,mse_beta1,mse_beta2,mse_betas_ols)
  #df <- rbind(mse_beta0,mse_beta1,mse_beta2)
  #return(df)
  return(df3)
  #this function row binds our results for the mses for the resp. betas for the resp. quantiles
  # and returns a vector with 15 rows (if we use 5 quantiles and 3 coefficients)
}

#Now: calculating the mean-squared errors using this function over a range of 
# sample sizes to see, whether the estimates tend to be unbiased (decreasing mse)
n_seq <- c(10,15,20,30,50,75,100,250,500,1000,2000) # a sequence of different sample sizes 

MSE_seq <- sapply(n_seq,mse_function_qr_n_seq) #storing the results
MSE_seq

#PLOTTING THE RESULTS
#First: create data frames
# For the mses of beta 0 estimates for all resp. quantiles
mse_beta0_df <- data.frame(sample_size = n_seq,
                          log_mse_beta0_tau_0.05 = log(MSE_seq[1,]),
                          log_mse_beta0_tau_0.25 = log(MSE_seq[2,]),
                          log_mse_beta0_tau_0.5 = log(MSE_seq[3,]),
                          log_mse_beta0_tau_0.75 = log(MSE_seq[4,]),
                          log_mse_beta0_tau_0.95 = log(MSE_seq[5,]),
                          log_mse_beta0_ols = log(MSE_seq[16,])
                          )
# For the mses of beta 1 estimates for all resp. quantiles
mse_beta1_df <- data.frame(sample_size = n_seq,
                           log_mse_beta1_tau_0.05 = log(MSE_seq[6,]),
                           log_mse_beta1_tau_0.25 = log(MSE_seq[7,]),
                           log_mse_beta1_tau_0.5 = log(MSE_seq[8,]),
                           log_mse_beta1_tau_0.75 = log(MSE_seq[9,]),
                           log_mse_beta1_tau_0.95 = log(MSE_seq[10,]),
                           log_mse_beta1_ols = log(MSE_seq[17,])
)
# For the mses of beta 2 estimates for all resp. quantiles
mse_beta2_df <- data.frame(sample_size = n_seq,
                           log_mse_beta2_tau_0.05 = log(MSE_seq[11,]),
                           log_mse_beta2_tau_0.25 = log(MSE_seq[12,]),
                           log_mse_beta2_tau_0.5 = log(MSE_seq[13,]),
                           log_mse_beta2_tau_0.75 = log(MSE_seq[14,]),
                           log_mse_beta2_tau_0.95 = log(MSE_seq[15,]),
                           log_mse_beta2_ols = log(MSE_seq[18,])
)

#Second: create plots for these (one per beta estimate)
# Beta0:
plot_mse_beta0 <- ggplot(data=mse_beta0_df,aes(x=sample_size))+
  ggtitle("Mean-Squared Error of beta0 for increasing sample size") +
  ylab("log of MSE") +
  xlab("Sample Size") +
  geom_line(aes(y=log_mse_beta0_tau_0.05,color="0.05 Quantile")) +
  geom_line(aes(y=log_mse_beta0_tau_0.25,color="0.25 Quantile")) +
  geom_line(aes(y=log_mse_beta0_tau_0.5,color="0.5 Quantile")) +
  geom_line(aes(y=log_mse_beta0_tau_0.75,color="0.75 Quantile")) +
  geom_line(aes(y=log_mse_beta0_tau_0.95,color="0.95 Quantile")) +
  geom_line(aes(y=log_mse_beta0_ols,color="ols")) +
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))

plot_mse_beta0 #plots the results for beta0

# Beta1:
plot_mse_beta1 <- ggplot(data=mse_beta1_df,aes(x=sample_size))+
  ggtitle("Mean-Squared Error of beta1 for increasing sample size") +
  ylab("log of MSE") +
  xlab("Sample Size") +
  geom_line(aes(y=log_mse_beta1_tau_0.05,color="0.05 Quantile")) +
  geom_line(aes(y=log_mse_beta1_tau_0.25,color="0.25 Quantile")) +
  geom_line(aes(y=log_mse_beta1_tau_0.5,color="0.5 Quantile")) +
  geom_line(aes(y=log_mse_beta1_tau_0.75,color="0.75 Quantile")) +
  geom_line(aes(y=log_mse_beta1_tau_0.95,color="0.95 Quantile")) +
  geom_line(aes(y=log_mse_beta1_ols,color="ols")) +
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))

plot_mse_beta1 #plots the results for beta1

# Beta2:
plot_mse_beta2 <- ggplot(data=mse_beta2_df,aes(x=sample_size))+
  ggtitle("Mean-Squared Error of beta2 for increasing sample size") +
  ylab("log of MSE") +
  xlab("Sample Size") +
  geom_line(aes(y=log_mse_beta2_tau_0.05,color="0.05 Quantile")) +
  geom_line(aes(y=log_mse_beta2_tau_0.25,color="0.25 Quantile")) +
  geom_line(aes(y=log_mse_beta2_tau_0.5,color="0.5 Quantile")) +
  geom_line(aes(y=log_mse_beta2_tau_0.75,color="0.75 Quantile")) +
  geom_line(aes(y=log_mse_beta2_tau_0.95,color="0.95 Quantile")) +
  geom_line(aes(y=log_mse_beta2_ols,color="ols")) +
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))

plot_mse_beta2 #plots the results for beta2

# (2) CALCULATING THE CONSISTENCY OF OUR ESTIMATES OVER DIFFERENT SAMPLE SIZES #
#Reminder:
#A parameter is consistent, if the probability that it will diverge (in absolute terms) from its true value 
# approaches zero with increasing sample size (information) n.
#Note, consistency: Prob(|Theta^{hat}_{n} - Theta| > epsilon) -> 0, as n goes to infinity (epsilon is a real number)
# or, equivalently: Prob(|Theta^{hat}_{n} - Theta| <= epsilon) -> 1, as n goes to infinity

#Defining a function that calculates the consistency probability (convergence in 
#probability) for different thresholds and sample sizes
func_consistency_prob_calc_new <- function(all_eps,n,S){
  n_eps <- length(all_eps)
  beta_hats_0_05_quantile <- matrix(NA,nrow=3,ncol=S) #empty matrix for the coefficient estimates for the 0.05 quantile
  beta_hats_0_25_quantile <- matrix(NA,nrow=3,ncol=S) #empty matrix for the coefficient estimates for the 0.25 quantile
  beta_hats_0_5_quantile <- matrix(NA,nrow=3,ncol=S) #empty matrix for the coefficient estimates for the 0.5 quantile
  beta_hats_0_75_quantile <- matrix(NA,nrow=3,ncol=S) #empty matrix for the coefficient estimates for the 0.75 quantile
  beta_hats_0_95_quantile <- matrix(NA,nrow=3,ncol=S) #empty matrix for the coefficient estimates for the 0.95 quantile
  beta_hats_ols <- matrix(NA,nrow=3,ncol=S) # empty matrix for the coefficient estimates of an ols model
  for(i in 1:S){
    df <- dgp_heteroscedastic(n=n, 
                              beta0=5000,
                              beta1=175,
                              beta2=-0.75) # every iteration we draw a new data set
    beta_hats_ols[,i] <- lm(Carbon_Emissions~poly(Income,2,raw=TRUE),
                            data=df)$coefficients
    beta_hats_0_05_quantile[,i] <- rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                                      data=df,
                                      tau=qs[1])$coefficients # stores the coefficient estimates for the 0.05-quantile
    beta_hats_0_25_quantile[,i] <- rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                                      data=df,
                                      tau=qs[2])$coefficients # stores the coefficient estimates for the 0.25-quantile
    beta_hats_0_5_quantile[,i] <- rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                                     data=df,
                                     tau=qs[3])$coefficients # stores the coefficient estimates for the 0.5-quantile
    beta_hats_0_75_quantile[,i] <- rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                                      data=df,
                                      tau=qs[4])$coefficients # stores the coefficient estimates for the 0.75-quantile
    beta_hats_0_95_quantile[,i] <- rq(Carbon_Emissions ~ poly(Income,2, raw=TRUE),
                                      data=df,
                                      tau=qs[5])$coefficients # stores the coefficient estimates for the 0.95-quantile
    
  }
  probs_tau_0.05 <-matrix(NA,nrow=3,ncol=n_eps)
  probs_tau_0.25 <-matrix(NA,nrow=3,ncol=n_eps)
  probs_tau_0.5 <-matrix(NA,nrow=3,ncol=n_eps)
  probs_tau_0.75 <-matrix(NA,nrow=3,ncol=n_eps)
  probs_tau_0.95 <-matrix(NA,nrow=3,ncol=n_eps)
  probs_ols <- matrix(NA,nrow=3,ncol=n_eps)
  for(z in 1:n_eps){
    for(u in 1:3){
      probs_tau_0.05[u,z] <- mean(abs(beta_hats_0_05_quantile[u,]-true_betas[u+1,1])<=abs(all_eps[z]*true_betas[u+1,1]))
      probs_tau_0.25[u,z] <- mean(abs(beta_hats_0_25_quantile[u,]-true_betas[u+1,2])<=abs(all_eps[z]*true_betas[u+1,2]))
      probs_tau_0.5[u,z] <- mean(abs(beta_hats_0_5_quantile[u,]-true_betas[u+1,3])<=abs(all_eps[z]*true_betas[u+1,3]))
      probs_tau_0.75[u,z] <- mean(abs(beta_hats_0_75_quantile[u,]-true_betas[u+1,4])<=abs(all_eps[z]*true_betas[u+1,4]))
      probs_tau_0.95[u,z] <- mean(abs(beta_hats_0_95_quantile[u,]-true_betas[u+1,5])<=abs(all_eps[z]*true_betas[u+1,5]))
      probs_ols[u,z] <- mean(abs(beta_hats_ols[u,]-model_coefficients[u])<=abs(all_eps[z]*model_coefficients[u]))
    }
  }
  output <- data.frame("Epsilon" = all_eps, #the epsilon values ("small" real number)
                       t(probs_tau_0.05), # the resp. probabilities of the resp. betas to differ more than epsilon from the true beta values for this quantile
                       t(probs_tau_0.25), 
                       t(probs_tau_0.5),
                       t(probs_tau_0.75),
                       t(probs_tau_0.95),
                       t(probs_ols)
  )
  output2 <- beta_hats_0_05_quantile
  return(output)
}
all_eps <- seq(0.01, 1, by = 0.01)# A sequence of epsilon (increasing)

################## CONSISTENCY FOR A LARGE SAMPLE SIZE #########################
# Calculating the consistency probability for a large sample: (takes approx.17 seconds)
all_eps <- seq(0.01, 1, by = 0.01)
consistency_prob <- func_consistency_prob_calc_new(all_eps = all_eps,
                                                     n = 2000,
                                                     S = 5000)
consistency_prob

#PLOTTING THE RESULTS
###### ALTERNATIVE PLOTTING #########
plot_consist_beta0 <- ggplot(data=consistency_prob,aes(x=Epsilon))+
  ggtitle("Consistency probability for the beta0 estimates for an increasing epsilon") +
  ylab("Probability") +
  xlab("Epsilon")+
  geom_line(aes(y=X1,color="0.05 Quantile"))+
  geom_line(aes(y=X1.1,color="0.25 Quantile"))+
  geom_line(aes(y=X1.2,color="0.5 Quantile"))+
  geom_line(aes(y=X1.3,color="0.75 Quantile"))+
  geom_line(aes(y=X1.4,color="0.95 Quantile"))+
  geom_line(aes(y=X1.5,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_consist_beta0
plot_consist_beta1 <- ggplot(data=consistency_prob,aes(x=Epsilon))+
  ggtitle("Consistency probability for the beta1 estimates for an increasing epsilon") +
  ylab("Probability") +
  xlab("Epsilon")+
  geom_line(aes(y=X2,color="0.05 Quantile"))+
  geom_line(aes(y=X2.1,color="0.25 Quantile"))+
  geom_line(aes(y=X2.2,color="0.5 Quantile"))+
  geom_line(aes(y=X2.3,color="0.75 Quantile"))+
  geom_line(aes(y=X2.4,color="0.95 Quantile"))+
  geom_line(aes(y=X2.5,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_consist_beta1
plot_consist_beta2 <- ggplot(data=consistency_prob,aes(x=Epsilon))+
  ggtitle("Consistency probability for the beta2 estimates for an increasing epsilon") +
  ylab("Probability") +
  xlab("Epsilon")+
  geom_line(aes(y=X3,color="0.05 Quantile"))+
  geom_line(aes(y=X3.1,color="0.25 Quantile"))+
  geom_line(aes(y=X3.2,color="0.5 Quantile"))+
  geom_line(aes(y=X3.3,color="0.75 Quantile"))+
  geom_line(aes(y=X3.4,color="0.95 Quantile"))+
  geom_line(aes(y=X3.5,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_consist_beta2


##### (3) CALCULATING THE COVERAGE PROBABILITY OVER DIFFERENT SAMPLE SIZES #####
#Reminder: Coverage probability is the number of times that the true parameter value
# falls within a confidence interval constructed from a random sample

#Defining a function that calculates the coverage prob. for our coefficient
# estimates depending on the sample size n
S_cov_prob = 2000 #repeated draws 
# This is now the alternative function for calculating the coverage probability,
# where it now also includes ols
func_cov_prob_calc_n_alt <- function(n){
  conf_intervals_tau_0.05 <- matrix(NA,ncol=3*2,nrow=S_cov_prob) #empty matrix for the confidence intervals for the beta estimates for the 0.05 quantile regression line
  conf_intervals_tau_0.25 <- matrix(NA,ncol=3*2,nrow=S_cov_prob) #empty matrix for the confidence intervals for the beta estimates for the 0.25 quantile regression line
  conf_intervals_tau_0.5 <- matrix(NA,ncol=3*2,nrow=S_cov_prob) #empty matrix for the confidence intervals for the beta estimates for the 0.5 quantile regression line
  conf_intervals_tau_0.75 <- matrix(NA,ncol=3*2,nrow=S_cov_prob) #empty matrix for the confidence intervals for the beta estimates for the 0.75 quantile regression line
  conf_intervals_tau_0.95 <- matrix(NA,ncol=3*2,nrow=S_cov_prob) #empty matrix for the confidence intervals for the beta estimates for the 0.95 quantile regression line
  conf_intervals_ols <- matrix(NA,ncol=3*2,nrow=S_cov_prob)
  for(i in 1:S_cov_prob){
    df <- dgp_heteroscedastic(n=n, 
                              beta0=5000,
                              beta1=175,
                              beta2=-0.75) # every iteration we draw a new data set
    for(u in 1:3){ #for every beta (since we have three beta coefficients)
      conf_intervals_ols[i,2*u-1] <-confint(lm(data=df,
                                               formula=Carbon_Emissions~poly(Income,2,raw=TRUE)),u,level=0.95)[1]
      conf_intervals_ols[i,2*u] <-confint(lm(data=df,
                                             formula=Carbon_Emissions~poly(Income,2,raw=TRUE)),u,level=0.95)[2]
      conf_intervals_tau_0.05[i,2*u-1] <- coef(rq(data=df,
                                                  formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                  tau=0.05,
                                                  ci=TRUE,
                                                  alpha=0.05))[u,2] #stores the lower bound for the confidence level
      conf_intervals_tau_0.05[i,2*u] <- coef(rq(data=df,
                                                formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                tau=0.05,
                                                ci=TRUE,
                                                alpha=0.05))[u,3] #stores the upper bound for the confidence level
      conf_intervals_tau_0.25[i,2*u-1] <- coef(rq(data=df,
                                                  formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                  tau=0.25,
                                                  ci=TRUE,
                                                  alpha=0.05))[u,2]
      conf_intervals_tau_0.25[i,2*u] <- coef(rq(data=df,
                                                formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                tau=0.25,
                                                ci=TRUE,
                                                alpha=0.05))[u,3]
      conf_intervals_tau_0.5[i,2*u-1] <- coef(rq(data=df,
                                                 formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                 tau=0.5,
                                                 ci=TRUE,
                                                 alpha=0.05))[u,2]
      conf_intervals_tau_0.5[i,2*u] <- coef(rq(data=df,
                                               formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                               tau=0.5,
                                               ci=TRUE,
                                               alpha=0.05))[u,3]
      conf_intervals_tau_0.75[i,2*u-1] <- coef(rq(data=df,
                                                  formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                  tau=0.75,
                                                  ci=TRUE,
                                                  alpha=0.05))[u,2]
      conf_intervals_tau_0.75[i,2*u] <- coef(rq(data=df,
                                                formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                tau=0.75,
                                                ci=TRUE,
                                                alpha=0.05))[u,3]
      conf_intervals_tau_0.95[i,2*u-1] <- coef(rq(data=df,
                                                  formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                  tau=0.95,
                                                  ci=TRUE,
                                                  alpha=0.05))[u,2]
      conf_intervals_tau_0.95[i,2*u] <- coef(rq(data=df,
                                                formula=Carbon_Emissions~ poly(Income,2, raw=TRUE),
                                                tau=0.95,
                                                ci=TRUE,
                                                alpha=0.05))[u,3]
    }
  }
  inside_ci_tau_0.05 <- matrix(NA,nrow=S_cov_prob,ncol=3) 
  inside_ci_tau_0.25 <- matrix(NA,nrow=S_cov_prob,ncol=3)
  inside_ci_tau_0.5 <- matrix(NA,nrow=S_cov_prob,ncol=3)
  inside_ci_tau_0.75 <- matrix(NA,nrow=S_cov_prob,ncol=3)
  inside_ci_tau_0.95 <- matrix(NA,nrow=S_cov_prob,ncol=3)
  cov_prob_tau_0.05 <- matrix(NA,nrow=3, ncol=1)
  cov_prob_tau_0.25 <- matrix(NA,nrow=3, ncol=1)
  cov_prob_tau_0.5 <- matrix(NA,nrow=3, ncol=1)
  cov_prob_tau_0.75 <- matrix(NA,nrow=3, ncol=1)
  cov_prob_tau_0.95 <- matrix(NA,nrow=3, ncol=1)
  inside_ci_ols <- matrix(NA,nrow=S_cov_prob,ncol=3)
  cov_prob_ols <- matrix(nrow=3,ncol=1)
  for(b in 1:3){ #for each beta
    for(i in 1:S_cov_prob){ 
      inside_ci_tau_0.05[i,b] <- (conf_intervals_tau_0.05[i,2*b-1]<= true_betas[b+1,1] & true_betas[b+1,1]<=conf_intervals_tau_0.05[i,2*b]) #boolean matrix (TRUE and FALSE), TRUE whenever the true beta is within the ci boundaries, FALSE otherwise
      inside_ci_tau_0.25[i,b] <- (conf_intervals_tau_0.25[i,2*b-1]<= true_betas[b+1,2] & true_betas[b+1,2]<=conf_intervals_tau_0.25[i,2*b])
      inside_ci_tau_0.5[i,b] <- (conf_intervals_tau_0.5[i,2*b-1]<= true_betas[b+1,3] & true_betas[b+1,3]<=conf_intervals_tau_0.5[i,2*b])
      inside_ci_tau_0.75[i,b] <- (conf_intervals_tau_0.75[i,2*b-1]<= true_betas[b+1,4] & true_betas[b+1,4]<=conf_intervals_tau_0.75[i,2*b])
      inside_ci_tau_0.95[i,b] <- (conf_intervals_tau_0.95[i,2*b-1]<= true_betas[b+1,5] & true_betas[b+1,5]<=conf_intervals_tau_0.95[i,2*b])
      inside_ci_ols[i,b] <- (conf_intervals_ols[i,2*b-1]<= model_coefficients[b] & model_coefficients[b]<=conf_intervals_ols[i,2*b])
    }
    cov_prob_tau_0.05[b,1] <- round(sum(inside_ci_tau_0.05[,b],na.rm=TRUE)/S_cov_prob,2) #calculates the coverage probability for the 0.05 quantile
    cov_prob_tau_0.25[b,1] <- round(sum(inside_ci_tau_0.25[,b],na.rm=TRUE)/S_cov_prob,2)
    cov_prob_tau_0.5[b,1] <- round(sum(inside_ci_tau_0.5[,b],na.rm=TRUE)/S_cov_prob,2)
    cov_prob_tau_0.75[b,1] <- round(sum(inside_ci_tau_0.75[,b],na.rm=TRUE)/S_cov_prob,2)
    cov_prob_tau_0.95[b,1] <- round(sum(inside_ci_tau_0.95[,b],na.rm=TRUE)/S_cov_prob,2)
    cov_prob_ols[b,1] <- round(sum(inside_ci_ols[,b],na.rm=TRUE)/S_cov_prob,2)
  }
  output <- rbind(cov_prob_tau_0.05,
                  cov_prob_tau_0.25,
                  cov_prob_tau_0.5,
                  cov_prob_tau_0.75,
                  cov_prob_tau_0.95,
                  cov_prob_ols)
  return(output)
}

#################################
# Calculating the coverage probability for the following sequence

n_seq_cov_prob <- c(50,60,80,100,150,200,300,500,800,1000,2000)
Cov_prob_seq <- sapply(X = n_seq_cov_prob, FUN = function(X) func_cov_prob_calc_n_alt(X)) #takes a few hours
Cov_prob_seq

#Plotting the results:
#First, create a data frame:
cov_prob_df <- data.frame(t(Cov_prob_seq),n_seq_cov_prob)
cov_prob_df # taking a look at the data frame

###### ALTERNATIVE PLOTTING #########
plot_cov_prob_beta0 <- ggplot(data=cov_prob_df,aes(x=n_seq_cov_prob))+
  ggtitle("Coverage probability for the beta0-estimates for increasing sample sizes for alpha = 0.05") +
  ylab("Probability") +
  xlab("Sample Size")+
  geom_line(aes(y=X1,color="0.05 Quantile"))+
  geom_line(aes(y=X4,color="0.25 Quantile"))+
  geom_line(aes(y=X7,color="0.5 Quantile"))+
  geom_line(aes(y=X10,color="0.75 Quantile"))+
  geom_line(aes(y=X13,color="0.95 Quantile"))+
  geom_line(aes(y=X16,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_cov_prob_beta0
plot_cov_prob_beta1 <- ggplot(data=cov_prob_df,aes(x=n_seq_cov_prob))+
  ggtitle("Coverage probability for the beta1-estimates for increasing sample sizes for alpha = 0.05") +
  ylab("Probability") +
  xlab("Sample Size")+
  geom_line(aes(y=X2,color="0.05 Quantile"))+
  geom_line(aes(y=X5,color="0.25 Quantile"))+
  geom_line(aes(y=X8,color="0.5 Quantile"))+
  geom_line(aes(y=X11,color="0.75 Quantile"))+
  geom_line(aes(y=X14,color="0.95 Quantile"))+
  geom_line(aes(y=X17,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_cov_prob_beta1
plot_cov_prob_beta2 <- ggplot(data=cov_prob_df,aes(x=n_seq_cov_prob))+
  ggtitle("Coverage probability for the beta2-estimates for increasing sample sizes for alpha = 0.05") +
  ylab("Probability") +
  xlab("Sample Size")+
  geom_line(aes(y=X3,color="0.05 Quantile"))+
  geom_line(aes(y=X6,color="0.25 Quantile"))+
  geom_line(aes(y=X9,color="0.5 Quantile"))+
  geom_line(aes(y=X12,color="0.75 Quantile"))+
  geom_line(aes(y=X15,color="0.95 Quantile"))+
  geom_line(aes(y=X18,color="ols"))+
  scale_color_manual(name="Quantiles",
                     values= c("0.05 Quantile" = "red",
                               "0.25 Quantile" = "green",
                               "0.5 Quantile" = "blue",
                               "0.75 Quantile" = "yellow",
                               "0.95 Quantile" = "black",
                               "ols"= "orange"))
plot_cov_prob_beta2


################################################################################
######################## SIMULATION STUDY PART II ##############################
################################################################################

#Now that we have shown that quantile regression seems to have a coverage probability converging to the 
# confidence level, is consistent - that is, the probability that the estimator deviates in absolute 
# terms more than epsilon converges to 0, resp. the probability, that the absolute value of this deviation
# is smaller than epsilon converges to 1, as well as that the method's mean-squared error decreases with 
# increasing sample size, we compare quantile regression and 
# ordinary least-squares regression (ols) in different settings, where in theory there should be an 
# advantage in using quantile regression

# (1) Homoscedastic data (starting point)
# (2) Heteroscedastic data
# (3) Outliers
# ((4) Quantile Crossing)


# 1) HOMOSCEDASTIC DATA (STARTING POINT)

#Defining a function that generates a homoscedastic data set 
func_data_benchmark <- function(n,beta0,beta1,beta2){
  x <- rnorm(n, mean = 0, sd =4) # independent variable (Income)
  eps <- rnorm(n, mean =0, sd = 20) # homoscedastic error term
  y <- beta0 + beta1 * x + beta2 * x^2 + eps # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}

#plotting the benchmark data
data_benchmark <- func_data_benchmark(30000,1,5,-1)
plot_data_benchmark <- ggplot(data = data_benchmark,
                    mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Homoscedastic Data (starting point)")+
  ylab("Y") +
  xlab("X") +
  geom_smooth(data = data_benchmark,
              aes(color= "OLS Regression"),
              linetype = 1,
              method = "lm",
              formula = y ~ poly(x,2, raw=TRUE)) +
  geom_quantile(data = data_benchmark,
                aes(color = "Median Quantile Regression"),
                linetype = 1,
                quantiles = 0.5, 
                formula = y ~ poly(x,2, raw=TRUE))+
  scale_fill_manual(values = c( "OLS Regression" = "darkgreen",
                                "Median Quantile_Regression" = "darkred"))+
  labs(color = "Regression Lines")

plot_data_benchmark #shows the plot

# Comparing quantile regression with ordinary least square regression

perform_comparison_qr_ols <- function(train_percentage,S){
  mse_ols <- matrix(NA,nrow=S,ncol=1)
  mse_rq <- matrix(NA,nrow=S,ncol=1)
  for(i in 1:S){
    df <- func_data_benchmark(n=10000,
                              beta0=1,
                              beta1=5,
                              beta2=-1)
    sample <- sample.int(n=nrow(df),size=floor(train_percentage*nrow(df)),replace=F)
    train <- df[sample,]
    test <- df[-sample,]
    fit_ols <- lm(data=train,
                  formula=Carbon_Emissions~poly(Income,2,raw=TRUE))
    fit_rq <- rq(data=train,
                 formula=Carbon_Emissions~poly(Income,2,raw=TRUE),
                 tau=0.5,
                 ci=TRUE,
                 alpha=0.05)
    fitted_values_ols <- predict(fit_ols,
                                 newdata=test)
    fitted_values_rq <- predict(fit_rq,
                                newdata=test)
    mse_ols[i,1] <- mean((test$Carbon_Emissions - fitted_values_ols)^2) #mse ols 
    mse_rq[i,1] <- mean((test$Carbon_Emissions - fitted_values_rq[1])^2) #mse rq 
  }
  ols_beats_rq <- sum(mse_ols<=mse_rq)/S #percentage of times ols has a smaller mse than rq
  output <- data.frame(mse_ols,
                       mse_rq)
  return(ols_beats_rq)
}
# Calculating the result (takes approx 30 seconds to calculate)
mse_comparison_1 <- perform_comparison_qr_ols(
  train_percentage=0.8,
  S=1000)
mse_comparison_1
#apparently ols outperforms quantile regression in this setting it 100% of all cases
#this should be consistent with theory, since without any ols assumption violated, 
#the leas squares estimator is BLUE, and thus the corresponding estimated values
#should depict this.

# 2) HETEROSCEDASTIC DATA

#Defining a function that generates heteroscedastic data
func_data_heteroscedastic <- function(n,beta0,beta1,beta2){
  x <- rnorm(n, mean = 20, sd =1) # independent variable (Income)
  #eps <- rnorm(n,mean=0,sd=2*x)
  #eps <- rnorm(n, sd = exp(2*abs(x)^2)*abs(x)) # heteroscedastic error term
  eps <- rnorm(n, sd = exp(0.5*abs(x*0.2)^2)) # heteroscedastic error term
  y <- beta0 + beta1 * x + beta2 * x^2 + eps*x # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}

#Plotting the data for visualization
data_heteroscedastic <- func_data_heteroscedastic(3000,1,5,-1)
plot_data_heteroscedastic <- ggplot(data = data_heteroscedastic,
                              mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Heteroscedastic Data") +
  ylab("Y") +
  xlab("X") +
  geom_smooth(data = data_heteroscedastic,
              aes(color= "OLS Regression"),
              linetype = 1,
              method = "lm",
              formula = y ~ poly(x,2, raw=TRUE)) +
  geom_quantile(data = data_heteroscedastic,
                aes(color = "Median Quantile Regression"),
                linetype = 1,
                quantiles = 0.5, 
                formula = y ~ poly(x,2, raw=TRUE))+
  scale_fill_manual(values = c( "OLS Regression" = "darkgreen",
                                "Median Quantile_Regression" = "darkred"))+
  labs(color = "Regression Lines")

plot_data_heteroscedastic #shows the plot

#Comparing ols and quantile regression once more
perform_comparison_qr_ols_hs <- function(train_percentage,S){
  mse_ols <- matrix(NA,nrow=S,ncol=1)
  mse_rq <- matrix(NA,nrow=S,ncol=1)
  for(i in 1:S){
    df <- func_data_heteroscedastic(n=3000,
                              beta0=1,
                              beta1=5,
                              beta2=-1)
    sample <- sample.int(n=nrow(df),size=floor(train_percentage*nrow(df)),replace=F)
    train <- df[sample,]
    test <- df[-sample,]
    fit_ols <- lm(data=train,
                  formula=Carbon_Emissions~poly(Income,2,raw=TRUE))
    fit_rq <- rq(data=train,
                 formula=Carbon_Emissions~poly(Income,2,raw=TRUE),
                 tau=0.5,
                 ci=TRUE,
                 alpha=0.05)
    fitted_values_ols <- predict(fit_ols,
                                 newdata=test)
    fitted_values_rq <- predict(fit_rq,
                                newdata=test)
    mse_ols[i,1] <- mean((test$Carbon_Emissions - fitted_values_ols)^2) #test-mse ols 
    mse_rq[i,1] <- mean((test$Carbon_Emissions - fitted_values_rq[1])^2) #test-mse rq 
  }
  ols_beats_rq <- sum(mse_ols<=mse_rq)/S #percentage of times ols has a smaller mse than rq
  output <- data.frame(mse_ols,
                       mse_rq)
  return(ols_beats_rq)
  #return(output) #just in case one wants to have a look on the mses over the S repeated draws
}
# Calculating the result (takes approx 30 seconds to calculate)
mse_comparison_2 <- perform_comparison_qr_ols_hs(
  train_percentage=0.8,
  S=2000)
mse_comparison_2
# Now, with this data generating process , quantile regression seems to outperform 
# ordinary least squares
# NOTE: i am not sure, whether this performance improvement is due to outliers,
# which even though not included directly, are apparent...

# 3) OUTLIERS
#Now we will include outliers explicitly into the dgp

func_data_outliers <- function(n,beta0,beta1,beta2){
  n_outliers <- 0.03*n #number of outliers
  x_no <- n-n_outliers #observations not being outliers
  x <- c(rnorm(n-n_outliers, mean = 0, sd =4),rnorm(n_outliers,mean=0,sd=4)) # independent variable (Income)
  eps <- c(rnorm(n-n_outliers, mean =0, sd = 20),runif(n_outliers,min=1000,max= 2000)) # homoscedastic error term with outliers
  y <- beta0 + beta1 * x + beta2 * x^2 + eps # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}
#Plotting the data for visualization
data_outliers <- func_data_outliers(1000,1,5,-1)
plot_data_outliers <- ggplot(data = data_outliers,
                                    mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Homoscedastic Data with Oultiers") +
  ylab("Y") +
  xlab("X") +
  geom_smooth(data = data_outliers,
              aes(color= "OLS Regression"),
              linetype = 1,
              method = "lm",
              formula = y ~ poly(x,2, raw=TRUE)) +
  geom_quantile(data = data_outliers,
                aes(color = "Median Quantile Regression"),
                linetype = 1,
                quantiles = 0.5, 
                formula = y ~ poly(x,2, raw=TRUE))+
  scale_fill_manual(values = c( "OLS Regression" = "darkgreen",
                                "Median Quantile_Regression" = "darkred"))+
  labs(color = "Regression Lines")

plot_data_outliers #shows the plot

# Once again comparing ols and quantile regression
perform_comparison_qr_ols_outliers <- function(train_percentage,S,n){
  mse_ols <- matrix(NA,nrow=S,ncol=1)
  mse_rq <- matrix(NA,nrow=S,ncol=1)
  for(i in 1:S){
    df_train <- func_data_outliers(n=n*(1-train_percentage),
                                    beta0=1,
                                    beta1=5,
                                    beta2=-1)
    df_test <- func_data_benchmark(n=n*train_percentage,
                                   beta0=1,
                                   beta1=5,
                                   beta2=-1)
    fit_ols <- lm(data=df_train,
                  formula=Carbon_Emissions~poly(Income,2,raw=TRUE))
    fit_rq <- rq(data=df_train,
                 formula=Carbon_Emissions~poly(Income,2,raw=TRUE),
                 tau=0.5,
                 ci=TRUE,
                 alpha=0.05)
    fitted_values_ols <- predict(fit_ols,
                                 newdata=df_test)
    fitted_values_rq <- predict(fit_rq,
                                newdata=df_test)
    mse_ols[i,1] <- mean((df_test$Carbon_Emissions - fitted_values_ols)^2) #mse ols 
    mse_rq[i,1] <- mean((df_test$Carbon_Emissions - fitted_values_rq[1])^2) #mse rq 
  }
  ols_beats_rq <- sum(mse_ols<=mse_rq)/S #percentage of times ols has a smaller mse than rq
  output <- data.frame(mse_ols,
                       mse_rq)
  return(ols_beats_rq)
  #return(output) #just in case one wants to have a look on the mses over the S repeated draws
}
# Calculating the result (takes approx 30 seconds to calculate)
mse_comparison_3 <- perform_comparison_qr_ols_outliers(
  train_percentage=0.8,
  S=3000,
  n=2000)
mse_comparison_3

# with outliers (this large) we can clearly see the advantage of quantile regression in terms
# of relative robustness to outliers compared with ordinary least squares regression.
# This can already be observed in the graph.

# 4) Quantile Crossing
# Reminder: Quantile Crossing can occur, because we are using some bounded data
# to estimate the (intercept and) slope parameters for the regression lines for
# different quantiles.
# i) For linear regression models of polynomial degree 2:
# The crossing issue can resp. will occur, when new observations are "too much" 
# outside the bounds of the data used for estimating the model's parameters, as 
# long as the regression coefficients (slope) of the respective quantile regression
# lines differ.
# ii) For a simple regression model:
# Crossing can occur, even within the bounds, because each quantile is discribed 
# only by a straight line, and these lines will have to cross somewhere as long 
# as the slope parameter estimates differ, obviously.


# 4.1 DGP WITH NO CROSSING ISSUE
# Introducing the general issue of crossing, when estimating several quantiles
# regression lines at the same time

# Quantiles of interest:
qs #as before

# Function that creates a homscedastic dgp
data_no_crossing_dgp <- function(n,beta0,beta1){
  x <- rnorm(n, mean = 10, sd =3) # independent variable (Income)
  eps <- rnorm(n, mean =0, sd = 10) # homoscedastic error term
  y <- beta0 + beta1 * x + eps # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}

#plotting
data_no_crossing <- data_no_crossing_dgp(3000,5,1)
plot_data_no_crossing <- ggplot(data = data_no_crossing,
                                mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Quantile Crossing with Homoscedastic Data")+
  ylab("Y") +
  xlab("X") +
  geom_quantile(data = data_no_crossing,
                aes(color = "Quantile Regression"),
                linetype = 1,
                quantiles = qs, 
                formula = y ~ x,
                colour = "red")
plot_data_no_crossing #shows the plot

# Defining a function, that uses Monte-Carlo simulations in order to create S 
# repeated data draws and estimates for the respective quantile regression line
# coefficients and calculating how often we observe the issue of quantile crossing

func_crossing_detection<- function(S,n,train_percentage){
  crossings <- matrix(NA,nrow=S,ncol=length(qs)) 
  for(i in 1:S){
    df <- data_no_crossing_dgp(n,5,1)
    sample <- sample.int(n=n,size=floor(train_percentage*n),replace=F)
    train <- df[sample,]
    test <- df[-sample,]
    fit <- rq(data=train,
              formula = Carbon_Emissions ~ Income,
              tau = qs)
    fitted_values <- predict(fit,newdata=test)
    crossings[i,1] <- any(fitted_values[,1] > fitted_values[,5]|
                            fitted_values[,1] > fitted_values[,4]|
                            fitted_values[,1] > fitted_values[,3]|
                            fitted_values[,1] > fitted_values[,2])
    crossings[i,2] <- any(fitted_values[,2] > fitted_values[,5]|
                            fitted_values[,2] > fitted_values[,4]|
                            fitted_values[,2] > fitted_values[,3])
    crossings[i,3] <- any(fitted_values[,3] > fitted_values[,5]|
                            fitted_values[,3] > fitted_values[,4])
    crossings[i,4] <- any(fitted_values[,4] > fitted_values[,5])
    crossings[i,5] <- any(fitted_values[,5] < fitted_values[,1]|
                            fitted_values[,5] < fitted_values[,2]|
                            fitted_values[,5] < fitted_values[,3]|
                            fitted_values[,5] < fitted_values[,4])
  }
  output <- list(data.frame("Quantile"= c("0.05","0.25","0.5","0.75","0.95"),
                            "Number_of_Crossings"=colSums(crossings)),frequency_of_crossings = max(colSums(crossings))/S)
  return(output)
}
quant_crossing_no_issue <- func_crossing_detection(
  S=1000,
  n=1000,
  train_percentage = 0.8)
quant_crossing_no_issue
# as we can see from the calculated values, we do not have any issues with quantile
# crossing.


# 4.2 DGP WITH CROSSING ISSUE (due to heteroscedastic data and a simple linear regression model)

# Creating a simple (only one slope parameter) heteroscedastic dgp
data_crossing_dgp <- function(n,beta0,beta1){
  x <- rnorm(n, mean = 10, sd =3) # independent variable (Income)
  eps <- rnorm(n, mean =0, sd = exp(0.2*x)) # homoscedastic error term
  y <- beta0 + beta1 * x + eps # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}

#plotting 
data_crossing <- data_crossing_dgp(3000,5,1)
plot_data_crossing <- ggplot(data = data_crossing,
                              mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Quantile Regression with Heteroscedastic Data")+
  ylab("Y") +
  xlab("X") +
  geom_quantile(data = data_crossing,
                aes(color = "Quantile Regression"),
                linetype = 1,
                quantiles = qs, 
                formula = y ~ x,
                colour = "red")
plot_data_crossing #shows the plot
# As we can see in the plot, the issue of crossing is apparent

# Same function as in 4.1, but now with a differing dgp
func_crossing_detection<- function(S,n,train_percentage){
  crossings <- matrix(NA,nrow=S,ncol=length(qs)) 
  for(i in 1:S){
    df <- data_crossing_dgp(n,5,1)
    sample <- sample.int(n=n,size=floor(train_percentage*n),replace=F)
    train <- df[sample,]
    test <- df[-sample,]
    fit <- rq(data=train,
              formula = Carbon_Emissions ~ Income,
              tau = qs)
    fitted_values <- predict(fit,newdata=test)
    crossings[i,1] <- any(fitted_values[,1] > fitted_values[,5]|
                            fitted_values[,1] > fitted_values[,4]|
                            fitted_values[,1] > fitted_values[,3]|
                            fitted_values[,1] > fitted_values[,2])
    crossings[i,2] <- any(fitted_values[,2] > fitted_values[,5]|
                            fitted_values[,2] > fitted_values[,4]|
                            fitted_values[,2] > fitted_values[,3])
    crossings[i,3] <- any(fitted_values[,3] > fitted_values[,5]|
                            fitted_values[,3] > fitted_values[,4])
    crossings[i,4] <- any(fitted_values[,4] > fitted_values[,5])
    crossings[i,5] <- any(fitted_values[,5] < fitted_values[,1]|
                            fitted_values[,5] < fitted_values[,2]|
                            fitted_values[,5] < fitted_values[,3]|
                            fitted_values[,5] < fitted_values[,4])
  }
  output <- list(data.frame("Quantile"= c("0.05","0.25","0.5","0.75","0.95"),
                            "Number_of_Crossings"=colSums(crossings)),frequency_of_crossings = max(colSums(crossings))/S)
  return(output)
}
quant_crossing_is_issue <- func_crossing_detection(
  S=1000,
  n=10000,
  train_percentage = 0.8)
quant_crossing_is_issue

# 4.3 DEALING WITH THE ISSUE OF CROSSING
# i) Log-transforming the explantory variable
data_crossing_dgp_log <- function(n,beta0,beta1){
  x <- log(rtruncnorm(n, a=0,b=Inf,mean = 10, sd =3)) # independent variable (Income)
  eps <- rnorm(n, mean =0, sd = exp(0.2*abs(x))) # homoscedastic error term
  y <- beta0 + beta1 * x + eps # the dgp (the true relationship of the regressor and the regressand; Carbon Emissions)
  data <- data.frame(Income = x, Carbon_Emissions = y) #storing the vector of y and corresponding x values in a data frame
  return(data) # function returns the data frame
}
data_crossing_log <- data_crossing_dgp_log(3000,5,1)
plot_data_crossing_log <- ggplot(data = data_crossing_log,
                             mapping = aes(x=Income, y=Carbon_Emissions)) +
  geom_point(
    colour = "blue",
    size = 0.5
  ) +
  ggtitle("Quantile Regression with Heteroscedastic Data (log-transformed) ")+
  ylab("Y") +
  xlab("X") +
  geom_quantile(data = data_crossing_log,
                aes(color = "Quantile Regression"),
                linetype = 1,
                quantiles = qs, 
                formula = y ~ x,
                colour = "red")
plot_data_crossing_log #shows the plot
# when looking at the plot, it seems that log transformation of the explanatory 
# variable has dealt with the issue of quantile crossing

# Same function as in 4.1, but now with a log-transformed dgp
func_crossing_detection_log<- function(S,n,train_percentage){
  crossings <- matrix(NA,nrow=S,ncol=length(qs)) 
  for(i in 1:S){
    df <- data_crossing_dgp_log(n,5,1)
    sample <- sample.int(n=n,size=floor(train_percentage*n),replace=F)
    train <- df[sample,]
    test <- df[-sample,]
    fit <- rq(data=train,
              formula = Carbon_Emissions ~ Income,
              tau = qs)
    fitted_values <- predict(fit,newdata=test)
    crossings[i,1] <- any(fitted_values[,1] > fitted_values[,5]|
                            fitted_values[,1] > fitted_values[,4]|
                            fitted_values[,1] > fitted_values[,3]|
                            fitted_values[,1] > fitted_values[,2])
    crossings[i,2] <- any(fitted_values[,2] > fitted_values[,5]|
                            fitted_values[,2] > fitted_values[,4]|
                            fitted_values[,2] > fitted_values[,3])
    crossings[i,3] <- any(fitted_values[,3] > fitted_values[,5]|
                            fitted_values[,3] > fitted_values[,4])
    crossings[i,4] <- any(fitted_values[,4] > fitted_values[,5])
    crossings[i,5] <- any(fitted_values[,5] < fitted_values[,1]|
                            fitted_values[,5] < fitted_values[,2]|
                            fitted_values[,5] < fitted_values[,3]|
                            fitted_values[,5] < fitted_values[,4])
  }
  output <- list(data.frame("Quantile"= c("0.05","0.25","0.5","0.75","0.95"),
                       "Number_of_Crossings"=colSums(crossings)),frequency_of_crossings = max(colSums(crossings))/S)
  return(output)
}
# This takes approx. 1.5 minutes to compute
quant_crossing_issue_log <- func_crossing_detection_log(
  S=1000,
  n=10000,
  train_percentage = 0.8)
quant_crossing_issue_log
# As we can see, the frequency of quantile crossings was significantly reduced
# by this transformation, but unfortunately it is not zero

