#### Load in the data #### 

# Load necessary libraries
library(vroom)        # For fast reading of data
library(tidyverse)    # For data manipulation
library(rstan)        # For Bayesian analysis
library(progress) 

install.packages("remotes")
remotes::install_github("dbdahl/cucumber/cucumber") # For progress bar
library(cucumber)  
library(reshape2) # For sampling functions
source('project_functions.R')  # Source custom functions

# Load in the data
data <- vroom("geocaches.csv")
data2 <- vroom("cali_geocaches.csv")
data <- rbind(data, data2)

# Convert certain columns to factors
data$terrain <- as.factor(data$terrain)
data$difficulty <- as.factor(data$difficulty)

# Create a new column to categorize terrain of geocaches
data <- data %>%
  filter(!is.na(favorites)) %>%
  group_by(terrain) %>%
  mutate(group_id = group_indices())

# Define the number of groups
n <- length(unique(data$group_id))

# Set fixed prior parameters
c <- 5
d <- 1
e <- 1
f <- 1

# Set initial values for unknown parameters
B <- 10000  # Number of iterations
alpha_samples <- matrix(rep(rep(0, B), n), ncol = n, nrow = B)
alpha_samples[1,] <- 0.327
beta_samples <- matrix(rep(rep(0, B), n), ncol = n, nrow = B)
beta_samples[1,] <- 2
lambda_samples <- matrix(rep(rep(0, B), n), ncol = n, nrow = B)
lambda_samples[1,] <- 1
alpha <- alpha_samples[1, 1]
lambda <- lambda_samples[1, 1]
beta <- beta_samples[1, 1]

# Create a progress bar
pb <- progress_bar$new(total = n * B)

# Extract favorites data for each group
y_list <- lapply(1:n, function(i) {
  data %>%
    filter(group_id == i) %>%
    pull(favorites) %>%
    as.numeric()
})

# MCMC sampling
for (i in 1:n) {
  y <- y_list[[i]]
  for (j in 2:B) {
    pb$tick()  # Increment progress bar
    beta <- log_full_conditional_beta(alpha = alpha_samples[j - 1, i], lambda = lambda_samples[j - 1, i])
    beta_samples[j, i] <- beta
    alpha <- slice_sampler_stepping_out(alpha_samples[j - 1, i], log_full_conditional_alpha, w = 5, log = TRUE)$x
    alpha_samples[j, i] <- alpha
    lambda <- slice_sampler_stepping_out(lambda_samples[j - 1, i], log_full_conditional_lambda, w = 5, log = TRUE)$x
    lambda_samples[j, i] <- lambda
  }
}

# Burn-in period
burn_in <- 1000
lambda_sample_posterior <- lambda_samples[-c(1:burn_in),]

### Posterior Summary
round(apply(lambda_sample_posterior, 2, quantile, probs = c(0.025, 0.975)), 4)



num_rows <- 9000
num_cols <- nrow(data)

# #plot the samples

plot.ts(lambda_sample_posterior, main = "Posterior Samples of Lambda: Terrain", xlab = "Iteration", ylab = "Lambda")
coda::effectiveSize(lambda_sample_posterior)
apply(lambda_sample_posterior,2,mean)


# Create a vector of colors
# Create a vector of colors
colors <- rainbow(9)

# Assuming lambda_sample_posterior is a data frame with lambda samples as columns
# Let's convert it to long format
lambda_sample_posterior <- as.data.frame(lambda_sample_posterior) # If not already a data frame
lambda_sample_posterior_long <- melt(lambda_sample_posterior)

# Rename columns for clarity
colnames(lambda_sample_posterior_long) <- c("Lambda", "Sample")

# Define labels for legend
legend_labels <- as.character(seq(1, 5, by = 0.5))

# Plot using ggplot2
ggplot(lambda_sample_posterior_long, aes(x = Sample, color = factor(Lambda))) +
  geom_density() +
  scale_color_manual(name = "Difficulty", values = colors, labels = legend_labels) +
  labs(title = "Densities of Difficulty Type", x = "Lambda", y = "Density") +
  theme_minimal() +
  # extend the x axis limits
  coord_cartesian(xlim = c(2.5, 11)) +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(shape = 15))) 

# Change shape of legend keys
# Change shape of legend keys
# Two parameters, intercept (alpha) and slope (beta).
k = 27
initial_parameters <- c(3, 6, 4, 3, 4, 2, 5, 1, 3)

big_like <- sapply(1:n, function(i) {
  param <- initial_parameters[i]
  mle <-  optim(par = param, fn = log_likelihood, y = y_list[[i]], control = list(fnscale = -1), method = "L-BFGS-B", lower = 0.01)$par
  log_likelihood(mle, y_list[[i]])
})


aic1 <- -2 * sum(big_like) + 2 * k

# DIC
sam <- sapply(1:n, function(i) {
  bayes_estimate <- mean(lambda_sample_posterior[,i])
  y = y_list[[i]]
  p_dic <- 2 * (log_likelihood(bayes_estimate,y) - mean(log_likelihood(lambda_sample_posterior[,i],y)))
  log <- log_likelihood(bayes_estimate, y)
  return(c(p_dic, log))
})

test1 <- function(observed, expected) { # Test general goodness-of-fit in 'residuals'
  sum((observed - expected)^2 / expected)
}

test2 <- function(observed, expected) { # Test for no monotonic trend in 'residuals'
  # Spearman correlation is a nonparametric alternative to the regular (Pearson) correlation
  cor(seq_along(observed), observed - expected, method = "spearman")
}



dic4 <- -2 * sum(sam[2,]) + 2 * sum(sam[1,])


lp <- sapply(1:n, function(i) {
  y <- y_list[[i]]
  samples <- lambda_sample_posterior[,i]
  sapply(1:length(y), \(i) 1/mean(1/dpois(y[i], samples)))
})

 lp <- sapply(lp,log)
 lp <- sapply(lp,sum)

#Table with AIC DIC and WAIC
df4 <- data.frame(AIC = aic, DIC = dic, WAIC = waic, LPML = lp)


### Frequentist ### 
# Initial guess for lambdas for all groups
initial_parameters <- c(3, 6, 4, 3, 4, 2, 5, 1, 3)

# Call optim for the first group
optim_results <- sapply(1:n, function(i) {
  optim(par = initial_parameters[i], fn = log_likelihood, y = y_list[[i]], control = list(fnscale = -1), method = "L-BFGS-B", lower = 0.01)$par
})

















