library(tidyverse)
library(haven)

data <- read_dta("broockman_kalla_replication_data.dta")
data <- data %>%
  filter(!is.na(treat_ind),
         !is.na(vf_age))

## Data transformation

all.dv.names.t1 <- c('miami_trans_law_t1', 'miami_trans_law2_t1', 'therm_trans_t1',
                     'gender_norm_sexchange_t1', 'gender_norm_moral_t1',
                     'gender_norm_abnormal_t1', 'gender_norm_trans_moral_wrong_t1')

all.dv.names.t2 <- c('miami_trans_law_t2', 'miami_trans_law2_t2', 'therm_trans_t2',
                     'gender_norm_sexchange_t2', 'gender_norm_moral_t2',
                     'gender_norm_abnormal_t2', 'gender_norm_trans_moral_wrong_t2')

all.dv.names.t3 <- c('miami_trans_law_withdef_t3', 'miami_trans_law2_withdef_t3', 
                     'therm_trans_t3', 'gender_norm_sexchange_t3', 'gender_norm_moral_t3',
                     'gender_norm_abnormal_t3','gender_norm_trans_moral_wrong_t3')

all.dv.names.t4 <- c('miami_trans_law_withdef_t4', 'miami_trans_law2_withdef_t4', 
                     'therm_trans_t4', 'gender_norm_sexchange_t4', 'gender_norm_moral_t4',
                     'gender_norm_abnormal_t4', 'gender_norm_trans_moral_wrong_t4')

## More stuff

trans.tolerance.dvs.t0 <- c('therm_trans_t0', 'gender_norms_sexchange_t0',
                            'gender_norms_moral_t0', 'gender_norms_abnormal_t0')

trans.tolerance.dvs.t1 <- c('therm_trans_t1', 'gender_norm_sexchange_t1',
                            'gender_norm_moral_t1', 'gender_norm_abnormal_t1',
                            'gender_norm_trans_moral_wrong_t1')

trans.tolerance.dvs.t2 <- c('therm_trans_t2', 'gender_norm_sexchange_t2',
                            'gender_norm_moral_t2', 'gender_norm_abnormal_t2',
                            'gender_norm_trans_moral_wrong_t2',
                            'trans_teacher_t2', 'trans_bathroom_t2')

trans.tolerance.dvs.t3 <- c('therm_trans_t3', 'gender_norm_sexchange_t3',
                            'gender_norm_moral_t3', 'gender_norm_abnormal_t3',
                            'gender_norm_trans_moral_wrong_t3',
                            'trans_teacher_t3', 'trans_bathroom_t3')

trans.tolerance.dvs.t4 <- c('therm_trans_t4', 'gender_norm_sexchange_t4',
                            'gender_norm_moral_t4', 'gender_norm_abnormal_t4',
                            'gender_norm_trans_moral_wrong_t4',
                            'trans_teacher_t4', 'trans_bathroom_t4')

## More stuff

trans.law.dvs.t0 <- c('miami_trans_law_t0', 'miami_trans_law2_t0')
trans.law.dvs.t1 <- c('miami_trans_law_t1', 'miami_trans_law2_t1')
trans.law.dvs.t2 <- c('miami_trans_law_t2', 'miami_trans_law2_t2')
# Note: Beginning with t3, the definition was added. 
trans.law.dvs.t3 <- c('miami_trans_law_withdef_t3', 'miami_trans_law2_withdef_t3')
trans.law.dvs.t4 <- c('miami_trans_law_withdef_t4', 'miami_trans_law2_withdef_t4')

##

gender.nonconformity.t0 <- c('gender_norm_looks_t0', 'gender_norm_rights_t0')
gender.nonconformity.t1 <- c('gender_norm_looks_t1', 'gender_norm_rights_t1')
# Note: Beginning with t2, an additional item was added to the measure. 
gender.nonconformity.t2 <- c('gender_norm_looks_t2', 'gender_norm_rights_t2', 
                             'gender_norm_dress_t2')
gender.nonconformity.t3 <- c('gender_norm_looks_t3', 'gender_norm_rights_t3', 
                             'gender_norm_dress_t3')
gender.nonconformity.t4 <- c('gender_norm_looks_t4', 'gender_norm_rights_t4', 
                             'gender_norm_dress_t4')

## Reverse coding 

reverse.coded.items <- c('gender_norms_moral_t0', 'gender_norm_moral_t1',
                         'gender_norm_moral_t2', 'gender_norm_moral_t3',
                         'gender_norm_moral_t4', 'gender_norms_abnormal_t0',
                         'gender_norm_abnormal_t1', 'gender_norm_abnormal_t2',
                         'gender_norm_abnormal_t3','gender_norm_abnormal_t4',
                         'gender_norm_trans_moral_wrong_t1', 
                         'gender_norm_trans_moral_wrong_t2',
                         'gender_norm_trans_moral_wrong_t3',
                         'gender_norm_trans_moral_wrong_t4',
                         'trans_bathroom_t2', 'trans_bathroom_t3',
                         'trans_bathroom_t4', 'gender_norm_looks_t0',
                         'gender_norm_looks_t1', 'gender_norm_looks_t2',
                         'gender_norm_looks_t3', 'gender_norm_looks_t4',
                         'gender_norm_rights_t0', 'gender_norm_rights_t1',
                         'gender_norm_rights_t2', 'gender_norm_rights_t3',
                         'gender_norm_rights_t4', 'gender_norm_dress_t2',
                         'gender_norm_dress_t3', 'gender_norm_dress_t4')
for(item in reverse.coded.items) data[,item] <- -1 * data[,item]

## Factor analysis

# Compute factor analysis outcome
compute.factor.dv <- function(dv.names, respondent.booleans, print.loadings = TRUE){
  responders <- data[respondent.booleans,]
  
  # Factor analysis
  factor.obj <- princomp(responders[, dv.names], cor=TRUE)
  if(print.loadings) print(loadings(factor.obj))
  dv <- as.vector(factor.obj$scores[,1])
  
  # More positive values on the factor should indicate more tolerance; reverse otherwise.
  if(cor(dv, responders$miami_trans_law_t0, use="complete.obs") < 0) dv <- -1 * dv
  
  # Put in the order of the main data frame
  dv.in.order <- dv[match(data$id, responders$id)]
  
  # Rescale to mean 0 sd 1 in placebo group; treatment effects can then be interpreted
  # as the effect in standard deviations the treatment would have among an untreated
  # population.
  dv.in.order <- (dv.in.order - mean(dv.in.order[!data$treat_ind], na.rm=TRUE)) /
    sd(dv.in.order[!data$treat_ind], na.rm=TRUE)
  
  return(as.vector(dv.in.order))
}

# In this code section we implement the procedures describe previously.

# First, misc. housekeeping.
# Recode age for small number of observations where it is missing.
data$vf_age[which(is.na(data$vf_age))] <- mean(data$vf_age, na.rm=TRUE)

# Language of interview
data$survey_language_es[is.na(data$survey_language_es)] <-
  data$survey_language_t0[is.na(data$survey_language_es)] == "ES"
data$survey_language_es[is.na(data$survey_language_es)] <- mean(data$survey_language_es, na.rm = TRUE)

# We subset to only those who came to door. contacted = came to door.
full.data <- data
data <- subset(data, contacted == 1)

# Compute the DVs in line with the above procedures.

# Omnibus DV of all primary outcomes.
data$all.dvs.t1 <- compute.factor.dv(all.dv.names.t1, 
                                     data$respondent_t1==1 & 
                                       !is.na(data$respondent_t1))
data$all.dvs.t2 <- compute.factor.dv(all.dv.names.t2, 
                                     data$respondent_t2==1 & 
                                       !is.na(data$respondent_t2))
data$all.dvs.t3 <- compute.factor.dv(all.dv.names.t3, 
                                     data$respondent_t3==1 & 
                                       !is.na(data$respondent_t3))
data$all.dvs.t4 <- compute.factor.dv(all.dv.names.t4, 
                                     data$respondent_t4==1 & 
                                       !is.na(data$respondent_t4))

# Trans tolerance DV.
data$trans.tolerance.dv.t0 <- compute.factor.dv(trans.tolerance.dvs.t0, 
                                                data$respondent_t0==1 & 
                                                  !is.na(data$respondent_t0))
data$trans.tolerance.dv.t1 <- compute.factor.dv(trans.tolerance.dvs.t1, 
                                                data$respondent_t1==1 & 
                                                  !is.na(data$respondent_t1))
data$trans.tolerance.dv.t2 <- compute.factor.dv(trans.tolerance.dvs.t2, 
                                                data$respondent_t2==1 & 
                                                  !is.na(data$respondent_t2))
data$trans.tolerance.dv.t3 <- compute.factor.dv(trans.tolerance.dvs.t3, 
                                                data$respondent_t3==1 & 
                                                  !is.na(data$respondent_t3))
data$trans.tolerance.dv.t4 <- compute.factor.dv(trans.tolerance.dvs.t4, 
                                                data$respondent_t4==1 & 
                                                  !is.na(data$respondent_t4))

# Law DV.
# Create outcome scale by averaging over the two questions.
data$miami_trans_law_t0_avg <- (data$miami_trans_law_t0 + 
                                  data$miami_trans_law2_t0)/2
data$miami_trans_law_t1_avg <- (data$miami_trans_law_t1 + 
                                  data$miami_trans_law2_t1)/2
data$miami_trans_law_t2_avg <- (data$miami_trans_law_t2 + 
                                  data$miami_trans_law2_t2)/2
# Note: Beginning with t3, the definition was added. 
data$miami_trans_law_t3_avg <- (data$miami_trans_law_withdef_t3 + 
                                  data$miami_trans_law2_withdef_t3)/2
# Note: Only one question was asked in t3 after the ad was shown, so no averaging is required.
data$miami_trans_law_t4_avg <- (data$miami_trans_law_withdef_t4 + 
                                  data$miami_trans_law2_withdef_t4)/2

# Gender Non-Conformity DV
data$gender_nonconformity_t0 <- compute.factor.dv(gender.nonconformity.t0, 
                                                  data$respondent_t0==1 & 
                                                    !is.na(data$respondent_t0))
data$gender_nonconformity_t1 <- compute.factor.dv(gender.nonconformity.t1, 
                                                  data$respondent_t1==1 & 
                                                    !is.na(data$respondent_t1))
data$gender_nonconformity_t2 <- compute.factor.dv(gender.nonconformity.t2, 
                                                  data$respondent_t2==1 & 
                                                    !is.na(data$respondent_t2))
data$gender_nonconformity_t3 <- compute.factor.dv(gender.nonconformity.t3, 
                                                  data$respondent_t3==1 & 
                                                    !is.na(data$respondent_t3))
data$gender_nonconformity_t4 <- compute.factor.dv(gender.nonconformity.t4, 
                                                  data$respondent_t4==1 & 
                                                    !is.na(data$respondent_t4))

# This dummy records whether the intervention was actually delivered vs. was not for any reason.
# Note that we do not use this variable to conduct comparisons only to measure successful contact rates.
data$treatment.delivered <- data$exp_actual_convo == "Trans-Equality" & !is.na(data$canvass_trans_ratingstart)

## Trans law outcomes

data <- data %>% 
  mutate(miami_trans_law_t0_avg = (miami_trans_law_t0 + 
                                     miami_trans_law2_t0)/2,
         miami_trans_law_t1_avg = (miami_trans_law_t1 + 
                                     miami_trans_law2_t1)/2,
         miami_trans_law_t2_avg = (miami_trans_law_t2 + 
                                     miami_trans_law2_t2)/2,
         miami_trans_law_t3_avg = (miami_trans_law_withdef_t3 + 
                                     miami_trans_law2_withdef_t3)/2,
         miami_trans_law_t4_avg = (miami_trans_law_withdef_t4 + 
                                     miami_trans_law2_withdef_t4)/2)



keep_vars <- c(
  'treat_ind',
  'miami_trans_law_t0_avg',
  'miami_trans_law_t3_avg',
  'therm_trans_t0',
  'therm_trans_t2',
  'gender_norm_moral_t1',
  "therm_obama_t1",
  "therm_obama_t0",
  'gender_norms_moral_t0',
  'vf_democrat',
  'ideology_t0', 
  'religious_t0',
  'exposure_trans_t0', 'pid_t0', 'vf_female', 'vf_hispanic',
  'vf_black', 'vf_age')



make_3cats <- function(v) {
  
  case_when(v < 50 ~ 0,
            v == 50 ~ 1,
            v > 50 ~ 2)
  


}

transphobia <- data |>
  select(all_of(keep_vars)) |>
  rename(
    nondiscrim_law_t0 = miami_trans_law_t0_avg,
    nondiscrim_law_t3 = miami_trans_law_t3_avg,
    gender_norm_moral_t0 = gender_norms_moral_t0, 
    treated = treat_ind
    ) |>
  mutate(
    nondiscrim_law_diff = nondiscrim_law_t3 - nondiscrim_law_t0,
    therm_trans_t2 = make_3cats(therm_trans_t2),
    therm_trans_t0 = make_3cats(therm_trans_t0)
  ) |>
  select(
    treated, nondiscrim_law_t3, therm_trans_t2,
     therm_obama_t1, gender_norm_moral_t1,
    ends_with("t0"), starts_with("vf"), everything() 
  )

save(transphobia, file = "../data/transphobia.rda")
write_csv(transphobia, file = "transphobia.csv")
