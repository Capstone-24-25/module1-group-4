library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(glmnet)
library(yardstick)

load('./data/biomarker-clean.RData')

biomarker_clean
################# Question 3 #####################

set.seed(1078422)
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
















































































































































































































































































################# Question 4 #####################
set.seed(1078422)

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

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
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## LASSO REGRESSION
##################

predictors <- biomarker_clean %>% 
  select(-c(group, ados)) %>% 
  as.matrix()

response <- biomarker_clean %>% 
  pull(group) %>% as.factor()

data_split <- initial_split(biomarker_clean, prop = 0.8)
train_data <- training(data_split)  
test_data <- testing(data_split)  

train_predictors <- train_data %>% select(-c(group, ados)) %>% as.matrix()
train_response <- train_data %>% pull(group) %>% as.factor()

lasso_model <- cv.glmnet(
  x = train_predictors,
  y = train_response,
  family = "binomial",
  alpha = 1,         
  nfolds = 10        
)

best_lambda <- lasso_model$lambda.min

lasso_coeffs <- coef(lasso_model, s = best_lambda)

lasso_df <- as.data.frame(as.matrix(lasso_coeffs)) %>%
  rownames_to_column("protein")

str(lasso_df)

lasso_df <- lasso_df %>%
  filter(protein != "(Intercept)") %>%
  rename(coefficient = `s1`) %>%  
  mutate(coefficient = as.numeric(coefficient))  

top_10_proteins <- lasso_df %>%
  mutate(abs_coef = abs(coefficient)) %>%
  arrange(desc(abs_coef)) %>%
  slice_head(n = 10) %>%
  pull(protein)

top_10_proteins

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, top_10_proteins)

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

print(proteins_sstar)
