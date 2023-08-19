# SET WHICH EXPERIMENT'S DATA TO ANALYSE (1,2 or 3)
experiment_n <- 3


# LOAD DATA SETS
## aggregate data, e.g. demographics
# data_meta <- read.csv("../final_analyses/master_fr_aggregate.csv") # adapt path as needed; fr = study 1
# data_meta <- read.csv("../final_analyses/master_rs_aggregate.csv") # adapt path as needed; rs = study 2
# data_meta <- read.csv("exp3_meta.csv") # study 3

data_meta <- setDT(data_meta) # make sure the data are of data.table type
setkey(data_meta, "PID")


## trials data
# data <- read.csv("../fr_topup_clean_long.csv") # adapt path as needed; fr = study 1
# data <- read.csv("../rs_clean_long.csv") # adapt path as needed; rs = study 2
# data <- read.csv("C:/Users/Sarah Zheng/dev/simpy-suspicion-model/exp3_data_trials.csv") # study 3

data <- setDT(data) # make sure the data are of data.table type
setkey(data, "PID")

data$subject_lost <- ifelse(data$win_lose_tie == "loss", TRUE, FALSE) # convert "win lose tie" to boolean factor variable "subject_lost"
data$outcome_blue <- ifelse(data$outcome == 1, TRUE, FALSE) # convert "outcome" to boolean factor variable indicating whether other participant reported blue, i.e., lying motivation

data <- data[, ..INCLUDE] # subset data to only necessary columns for efficiency

unique(data$PID) == unique(data_meta$PID) # SANITY CHECK
uuids <- data_meta$PID


# DEFINE MODELS TO FIT
params <- c("normed_signed_e_v", 
            "normed_unsigned_e_v",
            "subject_lied", 
            "subject_lost",
            "outcome_blue") # lying motivation, gave rise to signed vs. unsigned version of expectation violation
models <- define_models(params)
length(models)


# BAYESIAN MODEL AVERAGING PROCEDURE
cv_human_CIs <- cv_bayesian_avg_lmer(data, uuids, models, 16)
cv_perfect_CIs <- cv_bayesian_avg_perfect(data, uuids, MODELS_PERFECT, 16)


# INDIVIDUAL PARAMETER ESTIMATES WITH WINNING MODEL PARAMETERS
human_param_betas <- lm_fit_indiv_std(data, uuids)
human_param_betas$`(Intercept)` <- NULL


# COMPUTE SDT METRICS
human_sdt <- compute_sdt_metrics(data, uuids)


# COMPUTE PERFECT DETECTOR'S SDT METRICS
## the following perfect detector estimates are used to construct the perfect_model function
## ONLY add the significant parameters to the model!
mean(cv_perfect_CIs[cv_perfect_CIs$param == 'Signed expectation violation', ]$mean)
mean(cv_perfect_CIs[cv_perfect_CIs$param == 'Unsigned expectation violation', ]$mean)
mean(cv_perfect_CIs[cv_perfect_CIs$param == 'Lying oneself', ]$mean)
mean(cv_perfect_CIs[cv_perfect_CIs$param == 'Losing', ]$mean)
mean(cv_perfect_CIs[cv_perfect_CIs$param == 'Lying motivation', ]$mean)

perfect_sdt <- compute_perfect_sdt_metrics(data, uuids, experiment_n)


## compare human performance with that of perfect detector
wilcox.test(perfect_sdt$dprime, human_sdt$dprime)
wilcox.test(perfect_sdt$beta, human_sdt$beta)

dprime_long <- rbind(data.frame(type='perfect', dprime=perfect_sdt$dprime), data.frame(type='human', dprime=human_sdt$dprime))

ggplot(dprime_long, aes(x=type, y=dprime)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Lie detection accuracy human vs. perfect detector") +
  xlab("") + 
  scale_x_discrete(labels = c("human", "perfect detector"))
  

# participants for whom each predictor is significant
human_ps <- lm_fit_indiv_ps(data, uuids)
human_ps$X.Intercept. <- NULL

print('Proportion of participants for whom each model predictor is significant:')
sapply(human_ps, function(x) {sum(x < .05, na.rm = T)})/length(uuids)


# rerun model fits after excluding participants who thought there was no other player?
# TODO??


# SAVE PROCESSED DATA 
data_meta <- merge(data_meta, human_param_betas, by = 'PID')
data_meta <- merge(data_meta, human_sdt, by.x = 'PID', by.y = 'uuid')
# write.csv(cv_human_CIs, "exp3_cv_human_cis.csv")
# write.csv(cv_perfect_CIs, "exp3_cv_perfect_cis.csv")
# write.csv(data_meta, 'exp3_meta_complete.csv')
# write.csv(human_ps, 'exp3_human_ps.csv')
# write.csv(perfect_sdt, 'exp3_perfect_sdt.csv')


data_meta[is.na(data_meta)] <- 0
data_meta[data_meta==-Inf] <- 0
data_meta[data_meta==Inf] <- 0

data_meta$lie_prop <- (data_meta$lie_proportion_block_high_reward + data_meta$lie_proportion_block_low_reward)/2
data_meta$mean_suspicion_rating <- (data_meta$mean_suspicion_rating_block_high_reward + data_meta$mean_suspicion_rating_block_low_reward)/2

summary(lm('normed_signed_e_v ~ crt + 
           age + gender + edlev +
           lie_prop', data = data_meta))
summary(lm('normed_unsigned_e_v ~ crt + 
           age + gender + edlev +
           lie_prop', data = data_meta))
summary(lm('subject_liedTRUE ~ crt +
           age + gender + edlev +
           lie_prop', data = data_meta))
summary(lm('subject_lostTRUE ~ crt +
           age + gender + edlev +
           lie_prop', data = data_meta))
summary(lm('outcome_blueTRUE ~ crt +
           age + gender + edlev +
           lie_prop', data = data_meta))

summary(lm('mean_suspicion_rating ~ crt + 
           age + gender + edlev +
           lie_prop', data = data_meta))

summary(lm('dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt +
           age + gender + edlev +
           lie_prop', data = data_meta))



