# CROSS-VALIDATION BAYESIAN AVERAGING MODELLING
cv_bay_model_ci <- data.table()

for (i in 1:length(uuids)){
  print(i)
  uuid <- uuids[i]
  tmp <- data[data$PID != uuid, ]
  
  lmer_fits <- lmer_fit(tmp, models)
  lmer_fits$BIC_weights <- BIC_weights(lmer_fits$BIC)
  
  # BIC MODEL PROBABILITY WEGIHTS #
  weighted_params <- lmer_fits[, ..param_names] * lmer_fits$BIC_weights
  weighted_factors_avg <- apply(weighted_params, MARGIN = 2, function(x) { sum(x, na.rm = TRUE)}) 
  weighted_factors_sd <- apply(weighted_params, MARGIN = 2, function(x) { sd(x, na.rm = TRUE)}) 
  
  # 95%-CIs
  # TODO: check if need to replace .95 by z-value 1.96?
  ci_ub <- weighted_factors_avg + .95 * (weighted_factors_sd / sqrt(16)) # 16 = number of models in which each parameter occurs
  ci_lb <- weighted_factors_avg - .95 * (weighted_factors_sd / sqrt(16))
  
  human_CIs <- data.frame(excluded=uuid, mean=weighted_factors_avg, upper=ci_ub, lower=ci_lb)
  human_CIs$param <- c("Signed expectation violation", "Unsigned expectation violation", "Lying oneself", "Losing", "Lying motivation")
  human_CIs$param <- factor(human_CIs$param)
  
  cv_bay_model_ci <- rbindlist(list(human_CIs, cv_bay_model_ci))
  
  print("========================================================================")
}

cv_bay_model_ci[cv_bay_model_ci$upper <= 0, ]

ggplot(cv_bay_model_ci, aes(param, mean)) + 
  geom_violin() +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("Weighted standardised beta") +
  ggtitle("Factors driving human suspicion (cross-validated 95% CIs)") +
  scale_x_discrete(labels= c("Losing",
                             "Lying motivation",
                             "Lying oneself",
                             "Signed \n expectation violation",
                             "Unsigned expectation \n violation")) +
  theme_classic() +
  theme(text = element_text(size = 20))


# COMPUTE BETA CRITERION (BIAS) + PREDICT FROM PERSONL FACTORS AND WEIGHTS
df_sdt <- data.table()

for (uuid in uuids){
  tmp <- data[uuid]
  
  n_hit <- sum((tmp$suspicion_rating > .5) & (tmp$pp_lied == TRUE))
  n_targets <- sum(tmp$pp_lied == TRUE)
  hit_rate <- n_hit / n_targets
  
  n_fa <- sum((tmp$suspicion_rating > .5) & (tmp$pp_lied == FALSE))
  n_nolie <- sum(tmp$pp_lied == FALSE)
  fa_rate <- n_fa / n_nolie

  # dprime
  dprime <- qnorm(hit_rate) - qnorm(fa_rate)

  # beta
  zhr <- qnorm(hit_rate)
  zfar <- qnorm(fa_rate)
  beta <- exp(-zhr * zhr / 2 + zfar * zfar / 2)
  
  n_lied_self <- sum(tmp$col_reported != tmp$col_picked)

  sdt <- data.frame(uuid=uuid, dprime=dprime, beta=beta, hitrate=hit_rate, farate=fa_rate,
                    n_lied_self=n_lied_self)
  
  df_sdt <- rbindlist(list(sdt, df_sdt))
}

setkey(df_sdt, 'uuid')

data_long <- merge(data_long, df_sdt, by.x = 'PID', by.y = 'uuid')

master <- merge(data_long, rs_std_coefs, by.x = 'PID', by.y = 'uuid')
master$mean_suspicion_rating <- 6-master$mean_honesty_rating

master[is.na(master)] <- NA
master[master==-Inf] <- NA
master[master==Inf] <- NA
master$`(Intercept)` <- NULL

##### AGGREGATE LEVEL REGRESSIONS #####
dim(master)
head(master)

summary(lm('normed_signed_e_v ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('normed_unsigned_e_v ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('subject_liedTRUE ~ crt +
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

summary(lm('lied_min_sev ~ crt +
           age + gender + ed_lev +
           lie_prop', data = master))
summary(lm('lied_min_usev ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

# SIG & CONSISTENT EFFECTS:
summary(lm('mean_suspicion_rating ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
lm.beta(lm('mean_suspicion_rating ~ crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

rs_lm_dprime <- summary(lm('dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
rs_lm_dprime_std <- lm.beta(lm('dprime ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

rs_lm_dprime <- summary(lm('beta ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))
rs_lm_dprime_std <- lm.beta(lm('beta ~ outcome_blueTRUE + subject_liedTRUE + subject_lostTRUE + normed_signed_e_v + normed_unsigned_e_v +
            crt + 
           age + gender + ed_lev +
           lie_prop', data = master))

# write.csv(master, file="rs_master_meta.csv")


# RERUN MODEL AFTER EXCLUDING PARTICIPANTS WHO THOUGHT <3 OTHER PARTICIPANTS

