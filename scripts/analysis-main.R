library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

load('data/biomarker-clean.RData')

biomarker_clean

raw_data <- read.csv("data/biomarker-raw.csv")

#----------------------------Question 1----------------------
set.seed(1234)

#find 3 random proteins from the raw data
num_columns <- length(raw_data)
random_columns <- sample(2:num_columns, 3)
sample_proteins <- raw_data[, random_columns]
sample_proteins <- sample_proteins[-1,]
sample_proteins <- sample_proteins %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))


#look at the distributions of those proteins
for (i in 1:ncol(sample_proteins)) {
  protein_values <- as.numeric(as.character(sample_proteins[, i]))
  
  plot <- ggplot(data.frame(value=sample_proteins[, i]), aes(x=value)) + 
    geom_density(fill='blue', alpha=0.6) +
    labs(title=colnames(sample_proteins)[i], x='Protein Levels') +
    scale_x_continuous(breaks = pretty(range(protein_values, na.rm=T), n = 5)) +
    theme_minimal()
  print(plot)
}





#----------------------------Question 2----------------------

#If running from this file, the data is now back to the original where outliers were trimmed. Data / graphs will look different than in results.
#create function to identify outliers in the data. We use the mean + 3 standard deviations as the threshold
identify_outliers <- function(x) {
  abs_x <- abs(x)
  mean_x <- mean(abs_x, na.rm=T)
  sd_x <- sd(abs_x, na.rm=T)
  outlier_threshold <- mean_x + 3 * sd_x
  return(abs_x > outlier_threshold)
}


#sums the number of outliers per 
outlier_info <- biomarker_clean %>%
  select(-group, -ados) %>%
  mutate(across(everything(), identify_outliers)) %>%
  rowwise() %>%
  summarize(outliers_count = sum(c_across(everything()), na.rm=T),
            .groups = 'drop') %>%
  mutate(subject = row_number())


biomarker_with_outliers <- biomarker_clean %>%
  bind_cols(outlier_info)

outlier_summary <- biomarker_with_outliers %>% 
  group_by(group) %>%
  summarize(total_outliers = sum(outliers_count),
            .groups = 'drop')

outlier_summary2 <- biomarker_with_outliers %>%
  group_by(group) %>%
  summarize(mean_outliers = mean(outliers_count), .groups = 'drop')


# Plot aggregated data
ggplot(outlier_summary2, aes(x = group, y = mean_outliers, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Number of Outliers by Group", 
       x = "Group", 
       y = "Mean Number of Outliers") +
  theme_minimal()




#----------------------------Question 3----------------------
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

#-------------- larger number (more than ten) of top predictive proteins-----------------------------------------
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

# intersection
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


# fuzzy intersection ---------------------------------

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
  # Multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank_ttest = row_number(),
         p.adj = m * hm * p_value / rank_ttest) %>%
  # arrange by adjusted p-value
  arrange(p.adj)%>%
  mutate(rank_ttest = row_number())
  
  
  
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

# Fit model
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

# select top proteins 
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



#----------------------------Question 4----------------------
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

