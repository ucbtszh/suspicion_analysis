# install.packages(c("data.table", "feather", "lme4", "lmerTest", "lm.beta", "MuMIn","Rmisc", "ggpubr", "cowplot", "corrplot"))

# data formats
library(data.table)
library(feather)
# running linear mixed-effects models
library(lme4)
library(lmerTest)
library(lm.beta)
# additional stats
library(Rmisc)
# for plotting
library(ggpubr)
library(cowplot)
library(corrplot)

lmer_std_ests <- function(fit) {
  sdy <- sd(attr(fit, "resp")$y)
  sdx <- apply(attr(fit, "pp")$X, 2, sd)
  sc <- fixef(fit)*sdx/sdy
  std_estimates <- data.frame(t(sc))
  return(std_estimates)
}

## load datasets
data_long <- read.csv("C:/Users/Sarah Zheng/dev/simpy-suspicion-model/exp3_data_trials.csv")

data <- read.csv("C:/Users/Sarah Zheng/dev/simpy-suspicion-model/exp3_meta.csv")
data_meta <- read.csv("C:/Users/Sarah Zheng/dev/simpy-suspicion-model/exp3_meta.csv")

## process trials data
# subset data to only necessary columns
setDT(data_long)
setkey(data_long, "uuid")

data_long <- data_long[, 3:ncol(data_long)]
colnames(data_long)

data_long$subject_lost <- ifelse(data_long$subject_lost == 'True', TRUE, FALSE)
data_long$subject_lied <- ifelse(data_long$subject_lied == 'True', TRUE, FALSE)
data_long$random_pick_red <- ifelse(data_long$random_pick_col == 'red', 1, 0)
data_long$reward_block <- ifelse(data_long$reward_block == 'high', 1, 0)

## process meta data
setDT(data) # make sure the data are of data.table type
setkey(data, "uuid") # use subject IDs as keys
data$gender <- factor(data$gender)

# sum(data$lie_proportion_block_high_reward > data$lie_proportion_block_low_reward)

# check if people lied more on high reward block
t.test(data$lie_proportion_block_high_reward, data$lie_proportion_block_low_reward, paired = T, alternative = "g")
t.test(data$mean_suspicion_rating_block_high_reward, data$mean_suspicion_rating_block_low_reward, paired = T, alternative = "g")

t.test(data$lie_proportion_block_high_reward, data$lie_proportion_block_low_reward, paired = T)
t.test(data$mean_suspicion_rating_block_high_reward, data$mean_suspicion_rating_block_low_reward, paired = T)

# inspect how much more likely participants lied in the high stakes block if the manipulation worked  
tmp2 <- ifelse(data$lie_proportion_block_high_reward > data$lie_proportion_block_low_reward, (data$lie_proportion_block_high_reward+1) / (data$lie_proportion_block_low_reward+1), NA)
mean(na.omit(tmp2))

## store all unique PIDs
uuids <- unique(data_long$uuid)
length(uuids) # sanity check

## get trials subset of participants for whom manipulation worked
data_succ <- read.csv("C:/Users/Sarah Zheng/dev/simpy-suspicion-model/exp3_meta_success.csv")
uuids_succ <- unique(data_succ$uuid)
data_long_succ <- data_long[uuids_succ]

################################################
# DOUBLE CHECK HOW MANY PEOPLE NEVER LIED/LOST #
################################################

never_lied <- 0
for(id in uuids) {
  print(id)
  if(sum(data_long[id]$subject_lied) == 0) {
    never_lied <- never_lied + 1
  } else {
    print("participant lied at least once")
  }
}
never_lied

never_lost <- 0
for(id in uuids) {
  print(id)
  if(sum(data_long[id]$subject_lost) == 0) {
    never_lost <- never_lost + 1
  } else {
    print("participant lost at least once")
  }
}
never_lost

#################### 
####### PLOT ####### 
#################### 

ggplot(master_data) +
  geom_segment(aes(x=1, xend=2, y=mean_suspicion_rating_block_high_reward, yend=mean_suspicion_rating_block_low_reward, col=uuid), size=.75, show.legend=F) + xlim(.5, 2.5) 

data$diff_avg_suspicion <- data$mean_suspicion_rating_block_high_reward - data$mean_suspicion_rating_block_low_reward
data$diff_lie_prop <- data$lie_proportion_block_high_reward - data$lie_proportion_block_low_reward

ggplot(data, aes(x=diff_lie_prop, y=diff_avg_suspicion)) + 
  geom_point() +  # aes(col=uuid), show.legend=F
  geom_smooth(method="lm", color="darkgoldenrod") +
  xlab("\U0394 Tendency to lie \n(high - low reward block)") +
  ylab("\U0394 Mean suspicion rating \n(high - low reward block)") +
  theme_bw() + theme(text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# plot suspicion predictor estimates
sr_plotdata <- data.frame(model3_std_est)
sr_plotdata <- data.frame(t(sr_plotdata))
sr_plotdata$p <- c(data.frame(summary(fit3)$coefficients)$'Pr...t..')
sr_plotdata <- sr_plotdata[-1, ]
sr_plotdata$params <- c("normed_signed_e_v", "normed_unsigned_e_v", "subject_lied", "subject_lost", "reward_block", "condition")
sr_plotdata$params <- factor(sr_plotdata$params)
colnames(sr_plotdata) <- c("estimate", "p", "params")

ggplot(sr_plotdata, aes(x=factor(params, levels = c("normed_signed_e_v", "normed_unsigned_e_v", "subject_lied", "subject_lost", "reward_block", "condition")), y=estimate)) +
  geom_bar(stat="identity") +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  geom_point(data = sr_plotdata[sr_plotdata$p < .05, ], aes(x=params, y=0.35), shape = "*", size=10, color="black") +
  ylab("Standardised betas \npredicting suspicion") +
  xlab("") +
  scale_x_discrete(labels=c("Signed \n expectation \n violation", 
                            "Unsigned \n expectation \n violation", 
                            "Lynig oneself",
                            "Losing", 
                            "High stakes \n block", 
                            "Block order")) +
  theme_bw() + theme_bw() + theme(text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# install.packages("beeswarm")
library(beeswarm)

tmp <- melt(data, measure.vars = c("lie_proportion_block_high_reward", "lie_proportion_block_low_reward"))
beeswarm(value ~ variable, data=tmp, method="swarm", xlab = '', ylab = 'Tendency to lie', labels = c('High stakes block', 'Low stakes block'))
boxplot(value ~ variable, data = tmp, add = T, names = c("",""), col="#DAA52080")

tmp <- melt(data, measure.vars = c("mean_suspicion_rating_block_high_reward", "mean_suspicion_rating_block_low_reward"))
beeswarm(value ~ variable, data=tmp, method="swarm", xlab = '', ylab = 'Mean suspicion', labels = c('High stakes block', 'Low stakes block'))
boxplot(value ~ variable, data = tmp, add = T, names = c("",""), col="#DAA52080")

########### REPLICATING RESULTS #############

data <- data_long
data$pp_lied <- ifelse(data$pp_lied == "True", TRUE, FALSE)
str(data)

include <- c("uuid", # unique subject ID
             
             "normed_signed_e_v", # signed expectation violation, i.e., likelihood that the other player reported red or blue (continuous scale)
             "normed_unsigned_e_v", # surprise, i.e., unsigned expectation violation (continuous scale)
             
             "subject_lied", # whether the subject lied (TRUE) on the given trial or not (FALSE); Boolean var, will be treated as discrete variable
             "subject_lost", # whether the subject won (-1), tied (0) or lost (1) the trial
             
             "suspicion_rating", # human suspicion
             "pp_lied") # perfect lie detector

data <- data[, ..include] # subset data to only necessary columns
setkey(data, "uuid") # use subject IDs as keys

params <- c("normed_signed_e_v", 
            "normed_unsigned_e_v",
            "subject_lied", 
            "subject_lost")

## prepare binary grid to get all possible combinations of parameters
params_grid <- expand.grid(0:1, 0:1, 0:1, 0:1) # needs to correspond to number of unique parameters
names(params_grid) <- params # add parameters as column names

## check if every row has at least one parameter
rowSums(params_grid)
params_grid <- params_grid[-1, ] # first one has 0 parameters, so drop it
params_combi <- apply(params_grid, 1, function(i) which(i == 1)) # get list of parameters in each defined model

## write all model formulas in lmer format and add to list
models <- c()
for (combi in params_combi) {
  model <- paste("suspicion_rating ~", paste0(names(combi), collapse=' + '), "+ (1 +", paste0(names(combi), collapse=' + '), "| uuid)")
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
                    # R2=r.squaredGLMM(fit)) # requires MuMIn, only works in R >= 4.2
    fit_result <- rbindlist(list(fit_result, r))
  }
  return(cbind(fit_result, fit_estimates_std)) # combine all output before returning
}

BIC_weights <- function(BIC_array) {
  tmp_delta_bics <- abs(BIC_array - min(BIC_array))
  return(exp(-.5*tmp_delta_bics)/sum(exp(-.5*tmp_delta_bics)))
}

param_names <- c("normed_signed_e_v", "normed_unsigned_e_v", "subject_liedTRUE", "subject_lostTRUE")

# model fits
exp3_lmer_fits <- lmer_fit(data, models)
exp3_lmer_fits$BIC_weights <- BIC_weights(exp3_lmer_fits$BIC)
round(exp3_lmer_fits$BIC_weights, digits = 50)

# BIC MODEL PROBABILITY WEGIHTS #
exp3_tmp_BICweighted_params <- exp3_lmer_fits[, ..param_names] * exp3_lmer_fits$BIC_weights
exp3_weighted_factors <- apply(exp3_tmp_BICweighted_params, MARGIN = 2, function(x) { sum(x, na.rm = TRUE)})
exp3_weighted_factors_sd <- apply(exp3_tmp_BICweighted_params, MARGIN = 2, function(x) { sd(x, na.rm = TRUE)}) 

# 95%-CIs
exp3_ci_ub <- exp3_weighted_factors + .95 * (exp3_weighted_factors_sd / sqrt(8))
exp3_ci_lb <- exp3_weighted_factors - .95 * (exp3_weighted_factors_sd / sqrt(8))

exp3_human_CIs <- data.frame(mean=exp3_weighted_factors, upper=exp3_ci_ub, lower=exp3_ci_lb)
exp3_human_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing")
exp3_human_CIs$param <- factor(exp3_human_CIs$param)
str(exp3_human_CIs)

exp3_CIs <- ggplot(exp3_human_CIs, aes(param, mean)) + 
  geom_point() +  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving human suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))
exp3_CIs

#####
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
                           "pp_lied ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + subject_lost")

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
perfect_detector_fit$AIC_weights <- BIC_weights(perfect_detector_fit$AIC)

params <- c("normed_signed_e_v", "normed_unsigned_e_v", "subject_liedTRUE", "subject_lostTRUE")
perfect_detector_weighted_factors <- perfect_detector_fit[, ..params] * perfect_detector_fit$BIC_weights
perfect_detector_weighted_factors_avg <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sum(x, na.rm=T)})
perfect_detector_weighted_factors_sd <- apply(perfect_detector_weighted_factors, MARGIN = 2, function(x) {sd(x, na.rm=T)})
round(perfect_detector_weighted_factors_avg, digits = 3)

# 95%-CI perfect detector
# upper bound
perfect_ci_ub <- perfect_detector_weighted_factors_avg + .95 * (perfect_detector_weighted_factors_sd * sqrt(8)) # 8 = number of models in which each parameter is included
perfect_ci_lb <- perfect_detector_weighted_factors_avg - .95 * (perfect_detector_weighted_factors_sd * sqrt(8))

### study 1 perfect detector
exp3_perfect_CIs <- data.frame(mean=perfect_detector_weighted_factors_avg, upper=perfect_ci_ub, lower=perfect_ci_lb)
exp3_perfect_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing")
exp3_perfect_CIs$param <- factor(exp3_perfect_CIs$param)

exp3_CIs_perfect <- ggplot(exp3_perfect_CIs, aes(param, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving perfect detector's suspicion (95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))
exp3_CIs_perfect

##########
fit_lm_indiv_std <- function(data, uuids) {
  winning_formula <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied + subject_lost"
  winning_formula_noloss <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lied"
  winning_formula_nolies <- "suspicion_rating ~ normed_signed_e_v + normed_unsigned_e_v + subject_lost"
  
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
exp3_std_coefs <- fit_lm_indiv_std(data, uuids)
exp3_std_coefs$uuid <- uuids
setDT(exp3_std_coefs)
setkey(exp3_std_coefs, 'uuid')

master <- merge(data_meta, exp3_std_coefs, by.x = 'uuid', by.y = 'uuid')
colnames(master)

# predict
master[is.na(master)] <- 0
master[master==-Inf] <- 0
master[master==Inf] <- 0

summary(lm('normed_signed_e_v ~ crt + 
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))
summary(lm('normed_unsigned_e_v ~ crt + 
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))
summary(lm('subject_liedTRUE ~ crt +
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))
summary(lm('subject_lostTRUE ~ crt +
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))

# CONTRAST PARAMS: WEIGHT OBJECTIVE STAT COEFS VS. INTROSPECTION (SUBJECT LIED)
master$lied_min_sev <- master$subject_liedTRUE - master$normed_signed_e_v
master$lied_min_usev <- master$subject_liedTRUE - master$normed_unsigned_e_v
colnames(master)

summary(lm('lied_min_sev ~ crt +
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))
summary(lm('lied_min_usev ~ crt + 
           age + gender + edlev +
           lie_proportion_block_high_reward + lie_proportion_block_low_reward', data = master))

# SIG & CONSISTENT EFFECTS:
data_meta$lie_prop <- (data_meta$lie_proportion_block_high_reward + data_meta$lie_proportion_block_low_reward)/2
data_meta$mean_suspicion_rating <- (data_meta$mean_suspicion_rating_block_high_reward + data_meta$mean_suspicion_rating_block_low_reward)/2
data_meta$dprime <- (data_meta$dprime_block_high_reward + data_meta$dprime_block_low_reward)/2

summary(lm('mean_suspicion_rating ~ crt + 
           age + gender + edlev +
           lie_prop', data = data_meta))
lm.beta(lm('mean_suspicion_rating ~ crt + 
           age + gender + edlev +
           lie_prop', data = data_meta))

exp3_lm_dprime <- summary(lm('dprime ~ subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + edlev +
           lie_prop', data = master))
exp3_lm_dprime_std <- lm.beta(lm('dprime ~ subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + edlev +
           lie_prop', data = master))
exp3_lm_dprime_std <- data.frame(exp3_lm_dprime_std$standardized.coefficients)
exp3_lm_dprime_std <- exp3_lm_dprime_std[c("subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
exp3_lm_dprime_std <- data.frame(estim=exp3_lm_dprime_std, params=c("subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"))
exp3_lm_dprime_std$params <- as.factor(exp3_lm_dprime_std$params)

tmp <- exp3_lm_dprime$coefficients
tmp <- as.data.frame(tmp)
tmp <- tmp[c("subject_liedTRUE", "subject_lostTRUE", "normed_signed_e_v", "normed_unsigned_e_v", "lie_prop"), ]
exp3_lm_dprime_std$p <- tmp$`Pr(>|t|)`

exp3_dprime_task_params <- ggplot(exp3_lm_dprime_std, aes(x=factor(params, level=c("normed_unsigned_e_v", "normed_signed_e_v", "subject_lostTRUE", "lie_prop", "subject_lied")),
                                                      y=estim)) +
  geom_bar(stat="identity") +
  theme(legend.position="none", axis.title.x=element_blank()) + 
  geom_point(data = exp3_lm_dprime_std[exp3_lm_dprime_std$p < .05, ], aes(x=params, y=1), shape = "*", size=10, color="black") +
  ylab("Standardised betas \n predicting d'-scores") +
  xlab("") +
  scale_x_discrete(labels=c("Individual betas \n for unsigned \n expectation violation", 
                            "Individual betas \n for signed \n expectation violation", 
                            "Individual betas \n for losing", 
                            "Tendency \n to lie",
                            "Individual betas \n for lying oneself")) +
  theme_bw() + theme_bw() + theme(text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

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

# Partial regressions
mod1 <- lm('dprime ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            crt + 
           age + gender + edlev +
           lie_prop', data = master)
mod2 <- lm('subject_lostTRUE ~ subject_liedTRUE + 
            normed_signed_e_v + 
            normed_unsigned_e_v +
            crt +
           age + gender + edlev +
           lie_prop', data = master)

p1 <- plot_partial(mod1, mod2) + xlab("Losing \n (residuals)") + ylab("Lie detection \n d'-scores (residuals)")
p1

mod1 <- lm('dprime ~ subject_lostTRUE + 
            subject_liedTRUE + 
            normed_unsigned_e_v +
            crt +
           age + gender + edlev +
           lie_prop', data = master)
mod2 <- lm('normed_signed_e_v ~ subject_lostTRUE + 
            subject_liedTRUE +
            normed_unsigned_e_v + 
            crt + 
           age + gender + edlev +
           lie_prop', data = master)

p2 <- plot_partial(mod1, mod2) + xlab("Signed expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")
p2

mod1 <- lm('dprime ~ subject_lostTRUE + 
            subject_liedTRUE + 
            normed_signed_e_v +
            crt +
           age + gender + edlev +
           lie_prop', data = master)
mod2 <- lm('normed_unsigned_e_v ~ subject_lostTRUE + 
            subject_liedTRUE +
            normed_signed_e_v +
            crt + 
           age + gender + edlev +
           lie_prop', data = master)
p3 <- plot_partial(mod1, mod2) + xlab("Unsigned expectation \n violations (residuals)") + ylab("Lie detection \n d'-scores (residuals)")
p3

partial <- plot_grid(p3, p2, p1, ncol=3)
partial

##########################################
# CORRELATION MATRIX ACROSS PARTICIPANTS #
##########################################

# plot correlation matrix
data$subject_lied <- as.numeric(data$subject_lied)
data$subject_lost <- as.numeric(data$subject_lost)

param_corr_combis <- list(c("normed_signed_e_v", "normed_unsigned_e_v"), 
                          c("normed_signed_e_v", "subject_lied"), 
                          c("normed_signed_e_v", "subject_lost"), 
                          c("normed_unsigned_e_v", "subject_lied"), 
                          c("normed_unsigned_e_v", "subject_lost"), 
                          c("subject_lied", "subject_lost"))

master_corrs <- data.table()
for (uuid in uuids) {
  data_single <- data[uuid]
  # data_single$subject_lied <- as.numeric(data_single$subject_lied)
  # data_single$subject_lost <- as.numeric(data_single$subject_lost)
  
  corrs <- c()
  for (combi in param_corr_combis) {
    subset <- data_single[, ..combi]
    corr <- cor.test(subset[,1][[1]], subset[,2][[1]])$estimate
    corrs <- c(corrs, corr)
  }
  corrs <- data.frame(t(corrs))
  print(corrs)
  corrs$uuid <- uuid
  
  master_corrs <- rbindlist(list(corrs, master_corrs))
}

colnames(master_corrs)[1:6] <- c("signedEV_unsignedEV", "signedEV_lied",
                                 "signedEV_lost", "unsignedEV_lied",
                                 "unsignedEV_lost", "lied_lost")

corrs_sd <- lapply(master_corrs[,1:6], function(x) {sd(x, na.rm=T)})
corrs_avg <- lapply(master_corrs[,1:6], function(x) {mean(x, na.rm=T)})
corrs_test <- lapply(master_corrs[,1:6], function(x) {t.test(x, na.rm=T)})

setkey(master_corrs, "uuid")

p1 <- sapply(param_corr_combis, function(x) {x[1]})
p2 <- sapply(param_corr_combis, function(x) {x[2]})
corrs_avg <- data.frame(p2, avg=t(data.frame(corrs_avg)))

ggplot(corrs_avg, aes(x = p1, y = p2, fill = avg)) + 
  xlab("") +
  ylab("") +
  geom_tile() +
  geom_text(aes(label = round(avg, digits = 3)), color = "white", size = 4) +
  scale_fill_gradient(low = "grey", high = "red", name="Pearson's r") +
  scale_x_discrete(labels=c("Signed \n expectation \n violation", "Unsigned \n expectation \n violation", "Lying \n oneself")) +
  scale_y_discrete(labels=c("Unsigned \n expectation \n violation", "Lying \n oneself", "Losing")) +
  theme_bw() + theme(text = element_text(size=16))

### HISTOGRAMS ###
hist_mhr <- gghistogram(data_meta, x="mean_suspicion_rating", bins=50) + 
  theme(legend.position="none") +
  xlab("Mean suspicion ratings") +
  geom_vline(xintercept = mean(data_meta$mean_suspicion_rating, na.rm=T), 
             color = "black", size=1)

hist_lies <- gghistogram(data_meta, x="lie_prop", bins=50) + 
  theme(legend.position="none") +
  xlab("Tendency to lie") +
  geom_vline(xintercept = mean(data_meta$lie_prop), 
             color = "black", size=1)

hist_dprime <- gghistogram(data_meta, x="dprime", bins=50) + 
  theme(legend.position="none") +
  xlab("Lie detection accuracy (d'-scores)") +
  geom_vline(xintercept = mean(data_meta$dprime), 
             color = "black", size=1)

hist_sev <- gghistogram(data_long, x="e_v", bins=50) + 
  theme(legend.position="none") +
  xlab("Level of signed expectation violation") +
  geom_vline(xintercept = mean(data_long$e_v), 
             color = "black", size=1)

plot_grid(hist_mhr, hist_sev, hist_lies, hist_dprime, 
          ncol=2, labels = c("A", "B", "C", "D"))




mod1 <- lm('mean_suspicion_rating ~ crt + 
           age + gender + edlev', data = master)
mod2 <- lm('lie_prop ~ crt + 
           age + gender + edlev', data = master)
p4 <- plot_partial(mod1, mod2) + xlab("Tendency to lie (residuals)") + ylab("Mean suspicion (residuals)")
p4
