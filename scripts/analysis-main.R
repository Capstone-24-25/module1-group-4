library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load('data/biomarker-clean.RData')

biomarker_clean

#----------------------------problem 3----------------------
set.seed(101422)
biomarker_split <- biomarker_clean %>% 
  initial_split(prop = 0.8)
training_data <- training(biomarker_split)
testing_data <- testing(biomarker_split)
#-----------t test ----------
# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- training_data %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

##-----------random Forest ----------
# store predictors and response separately
predictors_train <- training_data %>%
  select(-c(group, ados))

response_train<- training_data %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors_train, 
                       y = response_train, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

##-----------LOGISTIC REGRESSION ----------
# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- training_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)


# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')
#---------------------------------------------------------------
ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 16) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 16) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')
# 
# After changing the number of top predictive proteins to 16, 
# We noticed the accuracy increased from 0.64 to 0.81,
# This indicates that the model became better at correctly classifying instances of both ASD and typically developing controls as more predictive proteins were included in the analysis.
# and the area under the curve also increased from 0.795 to 0.875.
# Sensitivity increased from 0.462 to 0.812. And specificity decreased from 0.833 to 0.8.

# fuzzy intersection ---------------------------------

# Conduct t-tests on the dataset
ttests_out <- biomarker_clean %>%
  # Drop ADOS score
  select(-ados) %>%
  # Arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # Nest by protein
  nest(data = c(level, group)) %>% 
  # Compute t-tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # Sort by p-value
  arrange(p_value) %>%
  # Multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank_ttest = row_number(),
         p.adj = m * hm * p_value / rank_ttest)

# Rank proteins from t-tests
ranked_ttest <- ttests_out %>%
  select(protein, rank_ttest)

## RANDOM FOREST
##################

# Store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% 
  pull(group) %>% 
  factor()

# Fit Random Forest model
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = TRUE)

# Compute Random Forest importance scores and ranks
ranked_rf <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(rank_rf = row_number()) %>%
  select(protein, rank_rf)

# Merge the rankings from both methods
combined_ranks <- ranked_ttest %>%
  full_join(ranked_rf, by = "protein") %>%
  # Handle missing ranks by assigning a high rank value (if protein not found in one of the methods)
  mutate(rank_ttest = ifelse(is.na(rank_ttest), max(rank_ttest, na.rm = TRUE) + 1, rank_ttest),
         rank_rf = ifelse(is.na(rank_rf), max(rank_rf, na.rm = TRUE) + 1, rank_rf)) %>%
  # Calculate a combined rank (average rank)
  mutate(composite_rank = (rank_ttest + rank_rf) / 2) %>%
  # Sort by the combined rank
  arrange(composite_rank)

# Select top proteins based on the fuzzy intersection 
proteins_sstar <- combined_ranks %>%
  slice_min(composite_rank, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# Prepare dataset for logistic regression
biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# Partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)


# Fit logistic regression model to the training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# Evaluate errors on the test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy, 
                            roc_auc)

# Make predictions on the test set and calculate metrics

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')
