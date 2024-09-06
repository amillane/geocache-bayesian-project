#### PPL ####
#### Load in the data #### 
library(vroom)
library(tidyverse)
library(rstan)

# Load in the data
data1 <- vroom("geocaches.csv")
data2 <- vroom("cali_geocaches.csv")
data <- rbind(data1, data2)

data$terrain <- as.factor(data$terrain)

#making a new column that categorized the terrain of the geocaches
data <- data %>%
  filter(!is.na(favorites)) %>%
  group_by(terrain) %>%
  mutate(group_id = group_indices())

# data <- data %>%
#   filter(!is.na(favorites)) %>%
#   group_by(difficulty) %>%
#   mutate(group_id = group_indices())


# Define the df 
df <- list(
  N = nrow(data),
  id = data$group_id,
  J = 9,
  y = data$favorites,
  a = 5,
  b = 1,
  c = 1,
  d = 1
)

n_cores <- 5
options(mc.cores = n_cores)

rstan_options(auto_write = TRUE)

stan_model <- stan_model(file = "Geocache.stan")

system.time(stan_fit <- stan(model_code = readLines("Geocache.stan"), 
                             data = df, 
                             chains = n_cores, 
                             iter = 5000, 
                             warmup = 1000, 
                             thin = 2))



samples <- rstan::extract(stan_fit)



for(i in 1:9){
  plot(samples$lambda[,i], type = "l")
}

#graph the density of each lambda column in samples
for(i in 1:9){
  density(samples$lambda[,i]) %>%
    plot(main = paste("Density of lambda", i))
}

apply(samples$lambda, 2, mean)
