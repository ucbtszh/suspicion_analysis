# install.packages(c("data.table", "feather", "lme4", "lmerTest", "lm.beta", "MuMIn", "Rmisc", "ggpubr", "cowplot", "corrplot"))

# data formats
library(data.table)
library(feather)
# running linear mixed-effects models
library(lme4)
library(lmerTest)
library(lm.beta)
# additional stats
# install.packages("MuMIn")
# library(MuMIn)
library(Rmisc)
# for plotting
library(ggpubr)
library(cowplot)
library(corrplot)

# double check working directory
setwd("D:/XUK5-2097647 - S Zheng - DATA/Behavioural cybersecurity/projects/cdg/5finalparams_clean")
getwd()

# study 1 data
data_long <- read.csv("../final_analyses/master_fr_aggregate.csv") # adapt path as needed; fr = study 1
data_long <- data_long[,1:12]

data <- read.csv("../fr_topup_clean_long.csv") # adapt path as needed; fr = study 1
exp1_lieprop <- data_long$lie_prop


# study 2 data
data_long <- read.csv("../final_analyses/master_rs_aggregate.csv") # adapt path as needed; fr = study 1
setDT(data_long)
setkey(data_long, "PID")
data_long <- data_long[,2:9]

data <- read.csv("../rs_clean_long.csv") # adapt path as needed; fr = study 1
exp2_lieprop <- data_long$lie_prop
data <- setDT(data) # make sure the data are of data.table type
colnames(data)

## convert "win lose tie" to binary factor (discrete) variable "subject_lost"
data$subject_lost <- ifelse(data$win_lose_tie == "loss", TRUE, FALSE)
data$outcome_blue <- ifelse(data$outcome == 1, TRUE, FALSE)
# data$subject_lost <- ifelse(data$win_lose_tie == "loss", 1, 0)

include <- c("PID", # unique subject ID
             
             "normed_signed_e_v", # signed expectation violation, i.e., likelihood that the other player reported red or blue (continuous scale)
             "normed_unsigned_e_v", # surprise, i.e., unsigned expectation violation (continuous scale)
             
             "subject_lied", # whether the subject lied (TRUE) on the given trial or not (FALSE); Boolean var, will be treated as discrete variable
             "subject_lost", # whether the subject won (-1), tied (0) or lost (1) the trial
             
             "outcome_blue", # whether the other player reported blue

             "suspicion_rating", # human suspicion
             "pp_lied") # perfect lie detector

data <- data[, ..include] # subset data to only necessary columns
setkey(data, "PID") # use subject IDs as keys

## store all unique PIDs
uuids <- unique(data$PID)
length(uuids) # sanity check: study 1 data has 102 subjects, study 2 has 108

################################################
# DOUBLE CHECK HOW MANY PEOPLE NEVER LIED/LOST #
################################################

never_lied <- 0
for(id in uuids) {
  print(id)
  if(sum(data[id]$subject_lied) == 0) {
    never_lied <- never_lied + 1
  } else {
    print("participant lied at least once")
  }
}
never_lied

never_lost <- 0
for(id in uuids) {
  print(id)
  if(sum(data[id]$subject_lost) == 0) {
    never_lost <- never_lost + 1
  } else {
    print("participant lost at least once")
  }
}
never_lost

############################
# BAYESIAN MODEL AVERAGING #
############################

## define model formulas
# define all possible models with five unique factors
# DEFINE MODEL FORMULAS TO FIT
## define list of parameters to include in constructing model formulas
params <- c("normed_signed_e_v", 
            "normed_unsigned_e_v",
            "subject_lied", 
            "subject_lost",
            "outcome_blue")

## prepare binary grid to get all possible combinations of parameters
params_grid <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1) # needs to correspond to number of unique parameters
names(params_grid) <- params # add parameters as column names

## check if every row has at least one parameter
rowSums(params_grid)
params_grid <- params_grid[-1, ] # first one has 0 parameters, so drop it
params_combi <- apply(params_grid, 1, function(i) which(i == 1)) # get list of parameters in each defined model

## write all model fomrulas in lmer format and add to list
models <- c()
for (combi in params_combi) {
  model <- paste("suspicion_rating ~", paste0(names(combi), collapse=' + '), "+ (1 +", paste0(names(combi), collapse=' + '), "| PID)")
  models <- c(models, model)
}
length(models) # check how many models we defined

## fit function
lmer_fit <- function(data, models) {
  # prepare empty data tables to append to
  fit_estimates_std <- data.table()
  fit_result <- data.table()
  
  for (i in 1:length(models)){
    model <- models[i] # get model formula from list of models
    print(paste0(i, " ", model)) # to check how far into the process we are when running the function
    
    tryCatch(
      {
        fit <- lmer(model, data = data, REML = FALSE, control=lmerControl(optimizer="bobyqa")) # fit linear mixed-effects model according to formula
      },
      warning=function(w) {
        print(w)
      }
    )
    
    # compute standardised parameter estimates
    sdy <- sd(attr(fit, "resp")$y)
    sdx <- apply(attr(fit, "pp")$X, 2, sd)
    sc <- fixef(fit)*sdx/sdy
    std_estimates <- data.frame(t(sc))
    fit_estimates_std <- rbindlist(list(fit_estimates_std, std_estimates), fill = TRUE) # add the std. estimates to the data.table
    
    # get overall fit indicators
    r <- data.table(model=model,
                    BIC=BIC(fit),
                    AIC=AIC(fit))
                    # R2=r.squaredGLMM(fit))
    fit_result <- rbindlist(list(fit_result, r))
  }
  return(cbind(fit_result, fit_estimates_std)) # combine all output before returning
}

BIC_weights <- function(BIC_array) {
  tmp_delta_bics <- abs(BIC_array - min(BIC_array))
  return(exp(-.5*tmp_delta_bics)/sum(exp(-.5*tmp_delta_bics)))
}

param_names <- c("normed_signed_e_v", "normed_unsigned_e_v", "subject_liedTRUE", "subject_lostTRUE", "outcome_blueTRUE")

#### TEST WITH non-normalised EV AND UNSIGNED EV VALUES
# fr_trials$unsigned_e_v <- abs(fr_trials$e_v)
# fr_trials$uuid <- fr_trials$PID
# fr_trials$subject_lost <- ifelse(fr_trials$win_lose_tie == 'loss', 1, 0)
# fr_trials$subject_lost <- as.factor(fr_trials$subject_lost)
# fr_trials$subject_lied <- as.factor(fr_trials$subject_lied)

# tmp_fits <- lmer_fit(fr_trials, models)
# tmp_fits$BIC_weights <- BIC_weights(tmp_fits$BIC)
# round(tmp_fits$BIC_weights, digits = 50)
# params <- c("e_v", 
#             "unsigned_e_v",
#             "subject_lied", 
#             "subject_lost")
# 
# param_names <- c("e_v", "unsigned_e_v", "subject_liedTRUE", "subject_lost1")
# tmp_BICw <- tmp_fits[, ..param_names] * tmp_fits$BIC_weights
# tmp_BICw_m <- apply(tmp_BICw, MARGIN = 2, function(x) { sum(x, na.rm = TRUE)})
# tmp_BICw_sd <- apply(tmp_BICw, MARGIN = 2, function(x) { sd(x, na.rm = TRUE)}) 
# 
# tmp_BICw_m + .95 * (tmp_BICw_sd / sqrt(8))
# tmp_BICw_m - .95 * (tmp_BICw_sd / sqrt(8))
# 
# models_perfect_agent <- c( "pp_lied ~ e_v",
#                            "pp_lied ~ unsigned_e_v",                                                                                                    
#                            "pp_lied ~ e_v + unsigned_e_v",                                                            
#                            "pp_lied ~ subject_lied",                                                                                                                  
#                            "pp_lied ~ e_v + subject_lied",                                                                          
#                            "pp_lied ~ unsigned_e_v + subject_lied",                                                                      
#                            "pp_lied ~ e_v + unsigned_e_v + subject_lied",                              
#                            "pp_lied ~ subject_lost",                                                                                                                  
#                            "pp_lied ~ e_v + subject_lost",                                                                          
#                            "pp_lied ~ unsigned_e_v + subject_lost",                                                                      
#                            "pp_lied ~ e_v + unsigned_e_v + subject_lost",                              
#                            "pp_lied ~ subject_lied + subject_lost",                                                                                    
#                            "pp_lied ~ e_v + subject_lied + subject_lost",                                            
#                            "pp_lied ~ unsigned_e_v + subject_lied + subject_lost",                                        
#                            "pp_lied ~ e_v + unsigned_e_v + subject_lied + subject_lost")
# 
# perfect_detector_fit <- lm_fit(fr_trials, models_perfect_agent)
# perfect_detector_fit$BIC_weights <- BIC_weights(perfect_detector_fit$BIC)
# perfect_detector_fit$AIC_weights <- BIC_weights(perfect_detector_fit$AIC)
# 
# params <- c("e_v", "unsigned_e_v", "subject_liedTRUE", "subject_lost1")
# perfect_detector_weighted_factors <- perfect_detector_fit[, ..params] * perfect_detector_fit$BIC_weights
# perfect_detector_weighted_factors_avg <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sum(x, na.rm=T)})
# perfect_detector_weighted_factors_sd <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sd(x, na.rm=T)})
# round(perfect_detector_weighted_factors_avg, digits = 3)
# 
# # 95%-CI perfect detector
# # upper bound
# perfect_ci_ub <- perfect_detector_weighted_factors_avg + .95 * (perfect_detector_weighted_factors_sd * sqrt(8)) # 8 = number of models in which each parameter is included
# perfect_ci_lb <- perfect_detector_weighted_factors_avg - .95 * (perfect_detector_weighted_factors_sd * sqrt(8))
# 


# study 1 model fits
fr_lmer_fits <- lmer_fit(data, models)
fr_lmer_fits$BIC_weights <- BIC_weights(fr_lmer_fits$BIC)
round(fr_lmer_fits$BIC_weights, digits = 50)

# BIC MODEL PROBABILITY WEGIHTS #
fr_tmp_BICweighted_params <- fr_lmer_fits[, ..param_names] * fr_lmer_fits$BIC_weights
fr_weighted_factors <- apply(fr_tmp_BICweighted_params, MARGIN = 2, function(x) { sum(x, na.rm = TRUE)}) 
fr_weighted_factors
fr_weighted_factors_sd <- apply(fr_tmp_BICweighted_params, MARGIN = 2, function(x) { sd(x, na.rm = TRUE)}) 
fr_weighted_factors_sd

# 95%-CIs
fr_ci_ub <- fr_weighted_factors + .95 * (fr_weighted_factors_sd / sqrt(8))
fr_ci_lb <- fr_weighted_factors - .95 * (fr_weighted_factors_sd / sqrt(8))

fr_human_CIs <- data.frame(mean=fr_weighted_factors, upper=fr_ci_ub, lower=fr_ci_lb)
fr_human_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing", "Lying motivation")
fr_human_CIs$param <- factor(fr_human_CIs$param)
str(fr_human_CIs)

study1_CIs <- ggplot(fr_human_CIs, aes(param, mean)) + 
  geom_point() +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving human suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying motivation",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))
study1_CIs

### study 2
rs_lmer_fits <- lmer_fit(data, models)
rs_lmer_fits$BIC_weights <- BIC_weights(rs_lmer_fits$BIC)
round(rs_lmer_fits$BIC_weights, digits = 3)

# BIC MODEL PROBABILITY WEGIHTS #
rs_weighted_params <- rs_lmer_fits[, ..param_names] * rs_lmer_fits$BIC_weights
rs_weighted_factors_avg <- apply(rs_weighted_params, MARGIN = 2, function(x) { sum(x, na.rm = TRUE)}) 
rs_weighted_factors_sd <- apply(rs_weighted_params, MARGIN = 2, function(x) { sd(x, na.rm = TRUE)}) 
rs_weighted_factors_avg
rs_weighted_factors_sd

# 95%-CIs
rs_ci_ub <- rs_weighted_factors_avg + .95 * (rs_weighted_factors_sd / sqrt(8))
rs_ci_lb <- rs_weighted_factors_avg - .95 * (rs_weighted_factors_sd / sqrt(8))

rs_human_CIs <- data.frame(mean=rs_weighted_factors_avg, upper=rs_ci_ub, lower=rs_ci_lb)
rs_human_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing", "Lying motivation")
rs_human_CIs$param <- factor(rs_human_CIs$param)

round(rs_human_CIs[,1:3], digits = 3)

study2_CIs <- ggplot(rs_human_CIs, aes(param, mean)) + 
  geom_point() +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving human suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying motivation",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))
study2_CIs


#################################################

models_perfect_agent <- c( "pp_lied ~ normed_signed_e_v",
                           "pp_lied ~ normed_unsigned_e_v",                                                                                                    
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v",                                                            
                           "pp_lied ~ subject_lied",                                                                                                                  
                           "pp_lied ~ normed_signed_e_v + subject_lied",                                                                          
                           "pp_lied ~ normed_unsigned_e_v + subject_lied",                                                                      
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied",                              
                           "pp_lied ~ subject_lost",                                                                                                                  
                           "pp_lied ~ normed_signed_e_v + subject_lost",                                                                          
                           "pp_lied ~ normed_unsigned_e_v + subject_lost",                                                                      
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lost",                              
                           "pp_lied ~ subject_lied + subject_lost",                                                                                    
                           "pp_lied ~ normed_signed_e_v + subject_lied + subject_lost",                                            
                           "pp_lied ~ normed_unsigned_e_v + subject_lied + subject_lost",                                        
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + subject_lost",
                           
                           "pp_lied ~ outcome_blue",
                           "pp_lied ~ normed_signed_e_v + outcome_blue",
                           "pp_lied ~ normed_unsigned_e_v + outcome_blue",                                                                                                    
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + outcome_blue",                                                            
                           "pp_lied ~ subject_lied + outcome_blue",                                                                                                                  
                           "pp_lied ~ normed_signed_e_v + subject_lied + outcome_blue",                                                                          
                           "pp_lied ~ normed_unsigned_e_v + subject_lied + outcome_blue",                                                                      
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + outcome_blue",                              
                           "pp_lied ~ subject_lost + outcome_blue",                                                                                                                  
                           "pp_lied ~ normed_signed_e_v + subject_lost + outcome_blue",                                                                          
                           "pp_lied ~ normed_unsigned_e_v + subject_lost + outcome_blue",                                                                      
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lost + outcome_blue",                              
                           "pp_lied ~ subject_lied + subject_lost + outcome_blue",                                                                                    
                           "pp_lied ~ normed_signed_e_v + subject_lied + subject_lost + outcome_blue",                                            
                           "pp_lied ~ normed_unsigned_e_v + subject_lied + subject_lost + outcome_blue",                                        
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + subject_lost + outcome_blue")

lm_fit <- function(data, models) {
  # prepare empty data tables to append to
  fit_estimates_std <- data.table()
  fit_ps <- data.table()
  fit_result <- data.table()
  
  for (i in 1:length(models)){
    model <- models[i] # get model formula from list of models
    print(paste0(i, " ", model)) # to check how far into the process we are when running the function
    
    tryCatch(
      {
        fit <- lm(model, data = data)
      }, warning=function(w) {
        print(w)
      }
    )
    # get standardised beta coefficients)
    coefs <- t(data.frame(lm.beta(fit)$standardized.coefficients)) 
    fit_estimates_std <- rbindlist(list(fit_estimates_std, data.table(coefs)), fill = T)
    
    # get significance (p-values) of each parameter in the model
    fit_sum <- data.frame(t(summary(fit)$coefficients))
    p <- fit_sum[nrow(fit_sum),] # p value is given in the last row
    fit_ps <- rbindlist(list(fit_ps, p), fill = TRUE) # add p-values to data.table; NOTE that the column names will be the same as the std. estimates data
    
    # get overall fit indicators
    r <- data.table(model=model,
                    BIC=BIC(fit),
                    AIC=AIC(fit),
                    R2=summary(fit)$r.squared)
    fit_result <- rbindlist(list(fit_result, r))
  }
  return(cbind(fit_result, fit_estimates_std)) # combine all output before returning ; fit_estimates_std ; fit_estimates_std
}

##### BEST MODEL FOR PERFECT LIE DETECTOR #####
perfect_detector_fit <- lm_fit(data, models_perfect_agent)
perfect_detector_fit$BIC_weights <- BIC_weights(perfect_detector_fit$BIC)
# perfect_detector_fit$AIC_weights <- BIC_weights(perfect_detector_fit$AIC)

params <- c("normed_signed_e_v", "normed_unsigned_e_v", "subject_liedTRUE", "subject_lostTRUE", "outcome_blueTRUE")
perfect_detector_weighted_factors <- perfect_detector_fit[, ..params] * perfect_detector_fit$BIC_weights
# perfect_detector_weighted_factors <- perfect_detector_fit[, ..params] * perfect_detector_fit$AIC_weights
perfect_detector_weighted_factors_avg <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sum(x, na.rm=T)})
perfect_detector_weighted_factors_sd <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sd(x, na.rm=T)})
round(perfect_detector_weighted_factors_avg, digits = 3)

# 95%-CI perfect detector
# upper bound
perfect_ci_ub <- perfect_detector_weighted_factors_avg + .95 * (perfect_detector_weighted_factors_sd * sqrt(8)) # 8 = number of models in which each parameter is included
perfect_ci_lb <- perfect_detector_weighted_factors_avg - .95 * (perfect_detector_weighted_factors_sd * sqrt(8))

### study 1 perfect detector
fr_perfect_CIs <- data.frame(mean=perfect_detector_weighted_factors_avg, upper=perfect_ci_ub, lower=perfect_ci_lb)
fr_perfect_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing", "Lying motivation")
fr_perfect_CIs$param <- factor(fr_perfect_CIs$param)
# round(fr_perfect_CIs[,1:3], digits = 3)
# round(fr_human_CIs[,1:3], digits = 3)

study1_CIs_perfect <- ggplot(fr_perfect_CIs, aes(param, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving perfect detector's suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying motivation",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))
study1_CIs_perfect


### study 2 perfect detector
rs_perfect_CIs <- data.frame(mean=perfect_detector_weighted_factors_avg, upper=perfect_ci_ub, lower=perfect_ci_lb)
rs_perfect_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing", "Lying motivation")
rs_perfect_CIs$param <- factor(rs_perfect_CIs$param)

round(rs_perfect_CIs[,1:3], digits = 3)

study2_CIs_perfect <- ggplot(rs_perfect_CIs, aes(param, mean)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving perfect detector's suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying motivation",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 18))
study2_CIs_perfect

#############################
# INDIVIDUAL BETA ESTIMATES #
#############################

fit_lm_indiv_std <- function(data, uuids) {
  winning_formula <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + subject_lost + outcome_blue"
  winning_formula_noloss <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + outcome_blue"
  winning_formula_nolies <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lost + outcome_blue"
  
  result <- data.table()
  
  for (i in 1:length(uuids)){
    uuid <- uuids[i]
    
    print(paste(i, uuid))
    
    data_subj <- data[uuid]
    
    formula <- ifelse(sum(data_subj$subject_lied) == 0, winning_formula_nolies, 
                      ifelse(sum(data_subj$subject_lost) == 0, winning_formula_noloss, winning_formula))
    
    data_subj$subject_lost <- as.factor(data_subj$subject_lost)
    data_subj$subject_lied <- as.factor(data_subj$subject_lied)
    
    fit <- lm(formula, data = data_subj) # fit winning model
    coefs <- t(data.frame(lm.beta(fit)$standardized.coefficients)) # get standardised beta coefficients)
    result <- rbindlist(list(result, data.table(coefs)), fill = TRUE) # add to list of individual fit results
  }
  return(result)
}

##### INDIV COEFS + MERGE WITH AGGREGATE DATA #####
fr_std_coefs <- fit_lm_indiv_std(data, uuids)
fr_std_coefs$uuid <- uuids
setDT(fr_std_coefs)
setkey(fr_std_coefs, 'uuid')

master_fr_std <- merge(data_long, fr_std_coefs, by.x = 'PID', by.y = 'uuid')
colnames(master_fr_std)

##### AGGREGATE LEVEL REGRESSIONS #####
master <- master_fr_std

master$normed_signed_e_v <- ifelse((master$normed_signed_e_v == "NA") | (master$normed_signed_e_v == Inf), 0, master$normed_signed_e_v)
master$normed_unsigned_e_v <- ifelse((master$normed_unsigned_e_v == "NA") | (master$normed_unsigned_e_v == Inf), 0, master$normed_unsigned_e_v)
master$subject_lostTRUE <- ifelse((master$subject_lostTRUE == "NA") | (master$subject_lostTRUE == Inf), 0, master$subject_lostTRUE)
master$outcome_blueTRUE <- ifelse((master$outcome_blueTRUE == "NA") | (master$outcome_blueTRUE == -Inf) | (master$outcome_blueTRUE == Inf), 0, master$outcome_blueTRUE)

summary(lm('normed_signed_e_v ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('normed_unsigned_e_v ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('subject_liedTRUE ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('subject_lostTRUE ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('outcome_blueTRUE ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))

# CONTRAST PARAMS: WEIGHT OBJECTIVE STAT COEFS VS. INTROSPECTION (SUBJECT LIED)
master$lied_min_sev <- master$subject_liedTRUE - master$normed_signed_e_v
master$lied_min_usev <- master$subject_liedTRUE - master$normed_unsigned_e_v

lm.beta(lm('lied_min_sev ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('lied_min_sev ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('lied_min_usev ~ crt + eq + aq + rgpts + # NO SIGNIFICANT PREDS, contrary to signed 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('mean_suspicion_rating ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))

### SIG + CONSISTENTE EFFECTS ###
lm.beta(lm('mean_suspicion_rating ~ crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))

fr_lm_dprime <- summary(lm('lie_detect_dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))
fr_lm_dprime_std <- lm.beta(lm('lie_detect_dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + eq + aq + rgpts + 
           age + gender + ed_lev +
           lie_prop', data = master))


# plot coefficients
fr_lm_dprime_std <- data.frame(fr_lm_dprime_std$standardized.coefficients)
fr_lm_dprime_std <- fr_lm_dprime_std[c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
fr_lm_dprime_std <- data.frame(estim=fr_lm_dprime_std, params=c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"))
tmp <- fr_lm_dprime$coefficients
tmp <- as.data.frame(tmp)
tmp <- tmp[c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
fr_lm_dprime_std$p <- tmp$`Pr(>|t|)`
fr_lm_dprime_std$params <- as.factor(fr_lm_dprime_std$params)

# setorder(fr_lm_dprime_std, cols="p")
tmp <- fr_lm_dprime_std
tmp$params <- factor(tmp$params)

fr_dprime_task_params <- ggplot(fr_lm_dprime_std, aes(x=factor(params, level=c("normed_unsigned_e_v", "normed_signed_e_v", "subject_lostTRUE", "outcome_blueTRUE", "subject_lied", "lie_prop")),
                y=estim)) +
  geom_bar(stat="identity") +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  geom_point(data = fr_lm_dprime_std[fr_lm_dprime_std$p < .05, ], aes(x=params, y=1), shape = "*", size=10, color="black") +
  ylab("Standardised betas \n predicting d'-scores") +
  xlab("") +
  scale_x_discrete(labels=c("Individual betas \n for unsigned \n expectation violation",
  "Individual betas \n for signed \n expectation violation",
  "Individual betas \n for losing",
  "Individual betas \n for inferred lying \n motivation",
  "Individual betas \n for lying oneself",
  "Tendency \n to lie")) +
  theme_bw() + theme_bw() + theme(text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
fr_dprime_task_params

################# STUDY 2 #################

##### INDIV COEFS + MERGE WITH AGGREGATE DATA #####http://127.0.0.1:13207/graphics/4fcb7f91-ec2b-4255-bf7f-349391f7a0db.png
data_long <- data_long[,2:9]
data_long

rs_std_coefs <- fit_lm_indiv_std(data, uuids)
rs_std_coefs$uuid <- uuids
setDT(rs_std_coefs)
setkey(rs_std_coefs, 'uuid')
master_rs_std <- merge(data_long, rs_std_coefs, by.x = 'PID', by.y = 'uuid')

master_rs_std[is.na(master_rs_std)] <- NA
master_rs_std[master_rs_std==-Inf] <- NA
master_rs_std[master_rs_std==Inf] <- NA

##### AGGREGATE LEVEL REGRESSIONS #####
master <- master_rs_std
dim(master)

summary(lm('normed_signed_e_v ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('normed_signed_e_v ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('normed_unsigned_e_v ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('subject_liedTRUE ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('subject_liedTRUE ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('subject_lostTRUE ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('outcome_blueTRUE ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))

# CONTRAST PARAMS: WEIGHT OBJECTIVE STAT COEFS VS. INTROSPECTION (SUBJECT LIED)
master$lied_min_sev <- master$subject_liedTRUE - master$normed_signed_e_v
master$lied_min_usev <- master$subject_liedTRUE - master$normed_unsigned_e_v
colnames(master)

master$mean_suspicion_rating <- 6-master$mean_honesty_rating

summary(lm('lied_min_sev ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('lied_min_sev ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))

summary(lm('lied_min_usev ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('lied_min_usev ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

# SIG & CONSISTENT EFFECTS:
summary(lm('mean_suspicion_rating ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('mean_suspicion_rating ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
rs_lm_dprime <- summary(lm('lie_detect_dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
rs_lm_dprime_std <- lm.beta(lm('lie_detect_dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

rs_lm_dprime_std <- data.frame(rs_lm_dprime_std$standardized.coefficients)
rs_lm_dprime_std <- rs_lm_dprime_std[c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
rs_lm_dprime_std <- data.frame(estim=rs_lm_dprime_std, params=c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"))
rs_lm_dprime_std$params <- as.factor(rs_lm_dprime_std$params)
tmp <- rs_lm_dprime$coefficients
tmp <- as.data.frame(tmp)
tmp <- tmp[c("outcome_blueTRUE", "subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
rs_lm_dprime_std$p <- tmp$`Pr(>|t|)`


rs_dprime_task_params <- ggplot(rs_lm_dprime_std, aes(x=factor(params, level=c("normed_unsigned_e_v", "normed_signed_e_v", "subject_lostTRUE", "lie_prop", "outcome_blueTRUE", "subject_lied")),
                                                      y=estim)) +
  geom_bar(stat="identity") +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  geom_point(data = fr_lm_dprime_std[fr_lm_dprime_std$p < .05, ], aes(x=params, y=1), shape = "*", size=10, color="black") +
  ylab("Standardised betas \n predicting d'-scores") +
  xlab("") +
  scale_x_discrete(labels=c("Individual betas \n for unsigned \n expectation violation", 
                            "Individual betas \n for signed \n expectation violation", 
                            "Individual betas \n for losing", 
                            "Tendency \n to lie",
                            "Individual betas \n for lying oneself")) +
  theme_bw() + theme_bw() + theme(text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


rs_dprime_task_params <- ggbarplot(data=rs_lm_dprime_std, x="params", y="estim") +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  geom_point(data = rs_lm_dprime_std[rs_lm_dprime_std$p < .05, ], aes(x=params, y=1), shape = "*", size=10, color="black") +
  ylab("Standardised betas predicting d'-scores") +
  scale_x_discrete(labels=c("Tendency \n to lie", "Signed expectation \n violation", 
                            "Unsigned \n expectation violation",
                            "Lying \n oneself", "Losing"))
rs_dprime_task_params

############################
# PLOT PARTIAL REGRESSIONS #
############################

plot_partial <- function(lm1, lm2) {
  resid.1 <- resid(lm1)
  resid.2 <- resid(lm2)

  tmp <- data.frame(resid.1, resid.2)
  
  plot <- ggplot(data=tmp, aes(x=resid.2, y=resid.1)) +
    geom_point() + 
    geom_smooth(method=lm, color="darkgoldenrod") +
    theme_bw() + theme(text = element_text(size=16), panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(plot)
}

str(master)

mod1 <- lm('lie_detect_dprime ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('subject_lostTRUE ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

p1 <- plot_partial(mod1, mod2) + xlab("Losing \n (residuals)") + ylab("Lie detection \n d'-scores (residuals)")

mod1 <- lm('lie_detect_dprime ~ subject_lost1 + 
            subject_liedTRUE + 
            normed_unsigned_e_v +
outcome1 +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('normed_signed_e_v ~ subject_lost1 + 
            subject_liedTRUE +
            normed_unsigned_e_v +
outcome1 +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

p2 <- plot_partial(mod1, mod2) + xlab("Signed expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")

mod1 <- lm('lie_detect_dprime ~ subject_lost1 + 
            subject_liedTRUE + 
            normed_signed_e_v +
            outcome1 +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('normed_unsigned_e_v ~ subject_lost1 + 
            subject_liedTRUE +
            normed_signed_e_v +
            outcome1 +
            crt + aq + eq + rgpts +
           age + gender + ed_lev +
           lie_prop', data = master)

p3 <- plot_partial(mod1, mod2) + xlab("Unsigned expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")

pg_fr <- plot_grid(p1, p2, p3, ncol=4)

# study 2 partial regressions
mod1 <- lm('lie_detect_dprime ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            outcome1 +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('subject_lost1 ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            outcome1 +
            crt +
           age + gender + ed_lev +
           lie_prop', data = master)

p1 <- plot_partial(mod1, mod2) + xlab("Losing \n (residuals)") + ylab("Lie detection \n d'-scores (residuals)")
p1

mod1 <- lm('lie_detect_dprime ~ subject_lost1 + 
            subject_liedTRUE + 
            normed_unsigned_e_v +
            outcome1 + 
            crt +
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('normed_signed_e_v ~ subject_lost1 + 
            subject_liedTRUE +
            normed_unsigned_e_v + outcome1 +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master)

p2 <- plot_partial(mod1, mod2) + xlab("Signed expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")

mod1 <- lm('lie_detect_dprime ~ subject_lost1 + 
            subject_liedTRUE + 
            normed_signed_e_v +
            outcome1 +
            crt +
           age + gender + ed_lev +
           lie_prop', data = master)

mod2 <- lm('normed_unsigned_e_v ~ subject_lost1 + 
            subject_liedTRUE +
            normed_signed_e_v +
            outcome1 +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master)

p3 <- plot_partial(mod1, mod2) + xlab("Unsigned expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")

pg_rs <- plot_grid(p1, p2, p3, ncol=4)
pg_rs

##########################################
# CORRELATION MATRIX ACROSS PARTICIPANTS #
##########################################

data$subject_lied <- as.numeric(data$subject_lied)

# plot correlation matrix
data$subject_lost <- as.numeric(data$subject_lost)
tmp <- data[,c("normed_signed_e_v", "normed_unsigned_e_v", "subject_lied", "subject_lost")]

param_corr_combis <- list(c("normed_signed_e_v", "normed_unsigned_e_v"), 
            c("normed_signed_e_v", "subject_lied"), 
            c("normed_signed_e_v", "subject_lost"), 
            c("normed_unsigned_e_v", "subject_lied"), 
            c("normed_unsigned_e_v", "subject_lost"), 
            c("subject_lied", "subject_lost"))

master_corrs <- data.table()
for (uuid in uuids) {
  data_single <- data[uuid]
  data_single$subject_lied <- as.numeric(data_single$subject_lied)
  data_single$subject_lost <- as.numeric(data_single$subject_lost)
  
  corrs <- c()
  for (combi in param_corr_combis) {
    subset <- data_single[, ..combi]
    corr <- cor.test(subset[,1][[1]], subset[,2][[1]])$estimate
    corrs <- c(corrs, corr)
  }
  corrs <- data.frame(t(corrs))
  print(corrs)
  corrs$PID <- uuid
  
  master_corrs <- rbindlist(list(corrs, master_corrs))
}

colnames(master_corrs)[1:6] <- c("signedEV_unsignedEV", "signedEV_lied",
                                 "signedEV_lost", "unsignedEV_lied",
                                 "unsignedEV_lost", "lied_lost")
str(data)

# t.test(data[data$subject_lied == 1, ]$normed_signed_e_v, data[data$subject_lied == 0, ]$normed_signed_e_v)
# t.test(data[data$subject_lost == 1, ]$normed_signed_e_v, data[data$subject_lost == 0, ]$normed_signed_e_v)

fr_corrs_sd <- lapply(master_corrs[,1:6], function(x) {sd(x, na.rm=T)})
fr_corrs_avg <- lapply(master_corrs[,1:6], function(x) {mean(x, na.rm=T)})
fr_corrs_test <- lapply(master_corrs[,1:6], function(x) {t.test(x, na.rm=T)})

rs_corrs_avg <- lapply(master_corrs[,1:6], function(x) {mean(x, na.rm=T)})
rs_corrs_test <- lapply(master_corrs[,1:6], function(x) {t.test(x, na.rm=T)})

setkey(master_corrs, "PID")

p1 <- sapply(param_corr_combis, function(x) {x[1]})
p2 <- sapply(param_corr_combis, function(x) {x[2]})
corrs_avg <- data.frame(p1, p2,
                        fr_avg=t(data.frame(fr_corrs_avg)),
                        rs_avg=t(data.frame(rs_corrs_avg)))

ggplot(corrs_avg, aes(x = p1, y = p2, fill = fr_avg)) + 
  xlab("") +
  ylab("") +
  geom_tile() +
  geom_text(aes(label = round(fr_avg, digits = 3)), color = "white", size = 4) +
  scale_fill_gradient(low = "grey", high = "red", name="Pearson's r") +
  scale_x_discrete(labels=c("Signed \n expectation \n violation", "Unsigned \n expectation \n violation", "Lying \n oneself")) +
  scale_y_discrete(labels=c("Unsigned \n expectation \n violation", "Lying \n oneself", "Losing")) +
  theme_bw() + theme(text = element_text(size=16))

ggplot(corrs_avg, aes(x = p1, y = p2, fill = rs_avg)) + 
  xlab("") +
  ylab("") +
  geom_tile() +
  geom_text(aes(label = round(rs_avg, digits = 3)), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "white", name="Pearson's r") +
  scale_x_discrete(labels=c("Signed \n expectation \n violation", "Unsigned \n expectation \n violation", "Lying \n oneself")) +
  scale_y_discrete(labels=c("Unsigned \n expectation \n violation", "Lying \n oneself", "Losing")) +
  theme_bw() + theme(text = element_text(size=16))
