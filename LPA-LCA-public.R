

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   A gentle introduction to latent class and latent profile analysis
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#
# Johannes Bauer
#

# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(cowplot)
library(mclust)
library(tidyLPA)


# ------------------------------------------------------------------------------
# Example 1: Univariate LPA of student achievement data
# ------------------------------------------------------------------------------


# --- Create a simulated data set of continuious data with two known classes A
# (30%) and B (70%) ---

set.seed(123456)
N1 = 300
N2 = 700
kca <- rnorm(N1, m = 15, sd = 3)
kcb <- rnorm(N2, m = 30, sd = 5)
mix <- c(kca, kcb)
dat1 <-
  data.frame(true_class = rep(c("A", "B"), times = c(300, 700)),
             value = mix)


# --- Plot empirical and population distributions (Fig. 1) ---

# Empirical
a <- ggplot(dat1, aes(x = value, y = ..density..)) +
  geom_histogram(fill = "lightgrey", binwidth = 2) +
  geom_density() +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density") +
  ylim(0, 0.06) +
  xlim(0, 50)

# Population
fmix <- function(x) .3 * dnorm(x, m = 15, sd = 3) + .7 * dnorm(x, m = 30, sd = 5)
fa   <- function(x) .3 * dnorm(x, m = 15, sd = 3)
fb   <- function(x) .7 * dnorm(x, m = 30, sd = 5)

b <- ggplot(data.frame(x = c(0:50)), aes(x = x)) +
  stat_function(fun = fmix) +
  stat_function(fun = fa, geom = "area", fill = 2, alpha = .2) +
  stat_function(fun = fb, geom = "area", fill = 4, alpha = .2) +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density") +
  ylim(0, 0.06) +
  xlim(0, 50)

plot_grid(a, b, labels = c("a", "b"))


# --- LPA using mclust 5 ---

# - Deciding number of classes -

# Model estimation
m1 <- Mclust(data = mix)
plot(m1, what = "BIC", modelName = "V")
LRT <-
  mclustBootstrapLRT(mix, modelName = "V", maxG = 3)  # Takes long!
LRT

# Table of fit indices
ll <- mclustLoglik(mclustBIC(mix))[, 2][1:4]
npara <- NULL
for (i in 1:4) {
  npara[i] <- nMclustParams(modelName = "V", d = 1, G = i)
}
npara
AIC <- 2 * npara - 2 * ll[1:4]
fit1 <- data.frame(
  nClasses = 1:4,
  Log.Likelihood = ll[1:4],
  AIC = AIC,
  BIC = m1$BIC[1:4, 2],
  ICL = mclustICL(mix)[1:4, 2],
  BLRT = c(NA, LRT$obs),
  p = c(NA, LRT$p.value)
)
round(fit1, 2)

# Get more fit indices from tidyLPA
estimate_profiles(mix, n_profiles = 1:4, variances = "varying") %>%
  get_fit()


# - Results for the two-class model -

# Mclust
summary(m1, parameters = TRUE)

# tidyLPA
estimate_profiles(mix, n_profiles = 2, variances = "varying") %>%
  get_estimates()

# Average latent class probabilities for most likely latent class membership
round(aggregate(
  x = 1 - m1$uncertainty,
  by = list(m1$classification),
  FUN = "mean"
) , 2)


# - Misclassification: True class vs. most likely class membership -
mlcm <- m1$classification
addmargins(table(dat1$true_class, mlcm))
classError(dat1$true_class, mlcm)



# ------------------------------------------------------------------------------
# Example 2: Multivariate LPA of work recovery
# ------------------------------------------------------------------------------

# We analyze simulated data for a LPA with three indicators of work recovery to
# mimic results given in Gillet et al. (2020), Study 1 (Tab S11). Variables
# are: overcommitment (overcom), psychological detachment (pdetach) and
# rumination (rumin). The last column contains class membership (lc). N = 500.

# For reading in the data, put the data file in the same folder as the R-script
# file.

wr <- read.table("wr_data.txt", header = TRUE)
head(wr)
boxplot(wr)



# --- LPA using mclust and tidyLPA---

# - Select best model and number of classes -

# Model fit from mclust
m1 <- Mclust(wr, G = 1:8)
summary(m1)
plot(m1, what = "BIC")  # Model fit for various parameterizations
plot(m1, what = "BIC", model = "VVI")  # Model fit for VVI-Model (= variances varying across classes, covariances restricted to zero).


# BLRT (we estimate only up to 5 classes estimation  takes very long!)
set.seed(123456)
LRT <- mclustBootstrapLRT(wr, modelName = "VVI", maxG = 5)  
LRT

# More fit indiced from tidyLPA
set.seed(123456)
m1t <- wr %>%
  estimate_profiles(n_profiles = 1:8,
                    variances = "varying",
                    covariances = "zero")
compare_solutions(m1t)  # What ist the best fitting model?
fitm1t <- get_fit(m1t)
fitm1t
fitm1t %>%
  select(Classes, AIC, CAIC, BIC, SABIC) %>%
  pivot_longer(cols = AIC:SABIC,
               names_to = "Index",
               values_to = "value") %>%
  ggplot(aes(x = Classes, y = value, shape = Index)) +
  geom_line() +
  geom_point(size = 2) +
  theme_few() +
  labs(shape = " Index")



# - Refit the identified best model with four latent classes -

# Get model parameters
m1.4 <- Mclust(wr, G = 4)
summary(m1.4, parameters = TRUE)

# Plot mean profiles of the classes
m1$parameters$mean %>%
  data.frame() %>%
  rename(
    Class1 = X1,
    Class2 = X2,
    Class3 = X3,
    Class4 = X4
  ) %>%
  mutate(Var = row.names(m1$parameters$mean)) %>%
  pivot_longer(!Var,
               names_to = "Class",
               values_to = "Mean") %>%
  ggplot(aes(
    x = Var,
    y = Mean,
    col = Class,
    group = Class
  )) +
  geom_hline(yintercept = 0, col = "grey60") +
  geom_line() +
  geom_point(aes(shape = Class), size = 2) +
  ylim(-3, 3) +
  scale_x_discrete(labels = c("Overcommitment",
                              "Psy. Detachment",
                              "Rumination")) +
  theme_few() +
  scale_color_few(palette = "Medium") +
  labs(x = "", col = "", shape = "")


# Average Latent Class Probabilities for Most Likely Latent Class Membership
round(aggregate(
  x = 1 - m1$uncertainty,
  by = list(m1$classification),
  FUN = "mean"
), 2)


# Plot mixture densities through dimensionality reduction
clust <- MclustDR(m1)

plot(clust, what = "boundaries", ngrid = 200)  # shows cluster boundaries:
# clusters are well separated
plot(clust, what = "density", dimens = 1)
