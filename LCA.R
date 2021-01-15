
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
library(tidyLPA)
library(ggthemes)
library(cowplot)
library(psych)
library(mclust)
# library(poLCA)

# ------------------------------------------------------------------------------
# Part 1: UNIVARIATE FINITE NORMAL MIXTURE MODEL
# ------------------------------------------------------------------------------

# cf. example in Sterba (2013)

# Create a data set with two known classes A and B
set.seed(123456)   
N1 = 300
N2 = 700
kca <- rnorm(N1, m = 15, sd = 3)
kcb <- rnorm(N2, m = 30, sd = 5)
joint <- c(kca, kcb)
dat1 <- data.frame(true_class = rep(c("A", "B"), times = c(300, 700)), 
                   value = joint)

describeBy(dat1$value, group = dat1$true_class)
boxplot(dat1$value ~ dat1$true_class)

datt1 <- tibble(dat1)

# Plot empirical distributions
plot(density(joint), col = 1, lty = 1,# within-class densities weighted by their respective class probabilities (Sterba, 2013, S. 780)
    type = "l", xlim = c(0, 50), ylim = c(0, .07), lwd = 2,
    main = "Illustrative univariate finite normal mixture density:\ntwo within class distributions and combined distibution",
    xlab = "Achievement", ylab = "Probability Density") 
lines(density(kca, weights = rep(1/length(joint), N1)), col = 4, lty = 2, lwd = 2)
lines(density(kcb, weights = rep(1/length(joint), N2)), col = 2, lty = 2, lwd = 2)

# Plot theoretical distributions
curve(.3 * dnorm(x, m = 15, sd = 3) + .7 *dnorm(x, m = 30, sd = 5), xlim = c(0, 50), col = "grey60", lty = 1, lwd = 2,
      main = "Illustrative univariate finite normal mixture density:\ntwo within class distributions and combined distibution",
      xlab = "Achievement", ylab = "Probability Density")
curve(.3 * dnorm(x, m = 15, sd = 3), col = 4, lty = 2, lwd = 2, add = TRUE)
curve(.7 * dnorm(x, m = 30, sd = 5), col = 2, lty = 2, lwd = 2, add = TRUE)

# Using ggplot: Empirical

# ggplot(dat1, aes(x = value )) +
#   geom_histogram(aes(y = ..density.., fill = true_class), position = "identity", alpha = .3) + 
#   geom_density(aes(x = value)) +
#   theme_few()
# 
# ggplot(dat1, aes(x = value )) +
#   geom_density(aes(y = .3 * ..density.., fill = true_class), position = "identity", alpha = .3) + 
#   geom_density(aes(x = value)) +
#   theme_few()

ggplot() +
  geom_density(dat = dat1, aes(x = value, y = ..density..), size = 1) + 
  geom_density(dat = dat1[  1:300, ], aes(x = value, y = .3 * ..density..), fill = 2, alpha = .3, col = 0) +
  geom_density(dat = dat1[301:700, ], aes(x = value, y = .7 * ..density..), fill = 4, alpha = .3, col = 0) +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density", title = "Illustrative univariate finite normal mixture density: within class and joint distibutions")


fjoint <- function(x) .3 * dnorm(x, m = 15, sd = 3) + .7 *dnorm(x, m = 30, sd = 5)
fa     <- function(x) .3 * dnorm(x, m = 15, sd = 3)
fb     <- function(x) .7 *dnorm(x, m = 30, sd = 5)

# ggplot(data.frame(x = c(0:50)), aes(x = x)) +
#   stat_function(fun = fjoint, aes(col = "Joint"), size = 1, lty = 1) +
#   stat_function(fun = fa, aes(col = "Class A"), size = 1, lty = 1) +
#   stat_function(fun = fb, aes(col = "Class B"), size = 1, lty = 1) +
#   theme_few() +
#   labs(colour = "Distribution")


# FIGURE 1

a <- ggplot(dat1, aes(x = value, y = ..density..)) +
  geom_histogram(fill = "lightgrey", binwidth = 2) +
  geom_density() +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density") +
  ylim(0, 0.06) +
  xlim(0, 50)
a

b <- ggplot(data.frame(x = c(0:50)), aes(x = x)) +
  stat_function(fun = fjoint) +
  stat_function(fun = fa, geom = "area", fill = 2, alpha = .2) +
  stat_function(fun = fb, geom = "area", fill = 4, alpha = .2) +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density") +
  ylim(0, 0.06) +
  xlim(0, 50)
b

plot_grid(a, b, labels = c("a", "b"))
ggsave(filename = "../Abbildungen/fig1_univ_mixt.png", height = 10, width = 22, units = "cm", dpi = 300, type = "cairo")



ggplot(data.frame(x = c(0:50)), aes(x = x)) +
  stat_function(fun = fjoint) +
  stat_function(fun = fa, geom = "area", fill = 2, alpha = .2) +
  stat_function(fun = fb, geom = "area", fill = 4, alpha = .2) +
  theme_few() +
  labs(x = "Achievement", y = "Probability Density") +
  ylim(0, 0.06) +
  xlim(0, 50) +
  geom_vline(xintercept = 21, col = "grey60", lty = 2) +
  geom_point(x = 21, y = .3 * dnorm(21, 15, 3), col = 2) +
  geom_point(x = 21, y = .7 * dnorm(21, 30, 5), col = 4) +
  geom_point(x = 21, y = (.3 * dnorm(21, 15, 3) + .7 * dnorm(21, 30, 5)), col = 1) 





# Univariate density estimation with mclust (https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html)

help(package="mclust")
?mclustModelNames

clust <- densityMclust(joint)
summary(clust, parameters = TRUE)
str(clust)
table(clust$classification)
true <- rep(c("kca", "kcb"), times = c(N1, N2))
classified <- clust$classification
table(true, classified)

classError(classified, true)

plot(clust, what = "BIC") # E = equal variances across classes, V = varying variances
# plot(clust, what = "diagnostic")
plot(clust, what = "density", data = joint, breaks = 25)
plot(clust, what = "diagnostic", type = "cdf")
plot(clust, what = "diagnostic", type = "qq")


# --- Clustering

# Using BIC
BIC <- mclustBIC(joint)
plot(BIC)
summary(BIC)
mclustLoglik(BIC)
ll <- mclustLoglik(BIC)[,2]

# Using ICL
ICL <- mclustICL(joint)
ICL
summary(ICL)
plot(ICL)
icl <- ICL[, 2]

# General clustering
m1 <- Mclust(joint, x = BIC)
summary(m1, parameters = TRUE)
str(m1)

# BLRT (takes long!)
set.seed(123456)
LRT <- mclustBootstrapLRT(joint, modelName = "V", maxG = 3)
LRT
str(LRT)
lrts <- c(NA, LRT$obs)
plrt <- c(NA, LRT$p.value)

# Compute AIC manually
npara <- NULL
for(i in 1:4) {npara[i] <- nMclustParams(modelName = "V", d = 1, G = i)}
npara

AIC <- 2*npara - 2*ll[1:4]

fit1 <- data.frame(nClasses = 1:4, 
                   Log.Likelihood = ll[1:4],
                   AIC = AIC,
                   BIC = m1$BIC[1:4, 2],
                   ICL = icl[1:4],
                   BLRT = lrts,
                   p = plrt)
round(fit1, 2)





plot(m1, what = "BIC") # E = equal variances across classes, V = varying variances
plot(m1, what = "classification")
plot(m1, what = "uncertainty")




table(m1$classification)
true <- rep(c("Class A", "Class B"), times = c(N1, N2))
classified <- m1$classification
addmargins(table(true, classified))
addmargins(prop.table(table(true, classified)))*100

err <- classError(classified, true)
err
length(err$misclassified)


# Average Latent Class Probabilities for Most Likely Latent Class Membership
round(aggregate(x = 1 - m1$uncertainty, by = list(m1$classification), FUN = "mean"), 2)
round(aggregate(x = m1$z, by = list(m1$classification), FUN = "mean"), 2)  # m1$z enthält die posterior probabilities

pposterior_cA <- m1$z[m1$classification == 1, 1]
pposterior_cB <- m1$z[m1$classification == 2, 2]
round(mean(pposterior_cA), 2)
round(mean(pposterior_cB), 2)


# Probabilities of class membership with score of 21 (Oberski, p. 279)
ptot <- .3 * dnorm(21, 15, 3) + .7 * dnorm(21, 30, 5)
pca <- .3 * dnorm(21, 15, 3) / ptot 
pcb <- .7 * dnorm(21, 30, 5) / ptot 
c(pca, pcb)


# using tidyLPA
m2 <- dat1$value %>% 
  estimate_profiles(1:4, variances = "varying") 
m2
get_fit(m2)

# using tidyLPA with Mplus and compare

# Takes very long!
# m3 <- datt1 %>% 
#   dplyr::select(value) %>% 
#   estimate_profiles(1:4, variances = "varying", package = "MplusAutomation",ANALYSIS = "starts = 2000, 500;")
# m3

# m3$model_2_class_4$warnings

# get_fit(m3)
# 
# plot(get_fit(m2)$BIC, type = "b", ylim = c(6800, 7050))
# lines(get_fit(m3)$BIC, type = "b", col = 2, lty = 2)
# lines(get_fit(m2)$AIC, type = "b", col = 3)
# lines(get_fit(m3)$AIC, type = "b", col = 4, lty = 2)


rm(list = ls())


# ------------------------------------------------------------------------------
# Example 2: Multivariate LPA of work recovery
# ------------------------------------------------------------------------------

# These are simulated data for a LPA with three indicators of work recovery to 
# mimic results given in Gillet et al. (2020), Study 1 (Tab S11). Variables
# are: overcommitment (overcom), psychological detachment (pdetach) and
# rumination (rumin). The last column contains class membership (lc). N = 500.

wr <- read.table("multi_lpa.dat", col.names = c("overcom",
                                                "pdetach",
                                                "rumin",
                                                "lc"))
wr <- wr[, -4]
head(wr)



# --- Clustering using mclust - 

# - Select best model
m1 <- Mclust(wr, G = 1:8)
summary(m1)
plot(m1, what = "BIC")

m1_icl <- mclustICL(wr, G = 1:8)
plot(m1_icl)


set.seed(123456)
LRT <- mclustBootstrapLRT(wr, modelName = "VVI", maxG = 5)  # Only for 5 classes sinc it takes long!
LRT


# Refit identified best model with four latent classes
m1 <- Mclust(wr, G = 4)
summary(m1, parameters = TRUE)

m1$parameters$mean %>% 
  data.frame() %>% 
  rename(Class1 = X1, 
         Class2 = X2, 
         Class3 = X3, 
         Class4 = X4) %>% 
  mutate(Var = row.names(m1$parameters$mean)) %>% 
  pivot_longer(!Var,
               names_to = "Class",
               values_to = "Mean") %>% 
  ggplot(aes(x = Var, y = Mean, col = Class, group = Class)) +
  geom_hline(yintercept = 0, col = "grey60") +
  geom_line() +
  geom_point(aes(shape = Class)) +
  ylim(-3, 3) +
  scale_x_discrete(labels = c("Overcommitment",
                            "Psy. Detachment",
                            "Rumination")) +
  theme_few() +
  scale_color_few(palette = "Medium") +
  labs(x = "", col = "", shape = "")

plot(m1, what = "classification")  
plot(m1, what = "uncertainty")


# Average Latent Class Probabilities for Most Likely Latent Class Membership
round(aggregate(x = 1 - m1$uncertainty, 
                by = list(m1$classification), 
                FUN = "mean"), 2)


# Plot mixture densities through dimensionality reduction
clust <- MclustDR(m1)
plot(clust, what = "boundaries", ngrid = 200)  # shows cluster boundaries: clusters 
                                               # are clearly separated
plot(clust, what = "density", dimens = 1)


# - Using tidyLPA -

m4 <- wr %>% 
  estimate_profiles(n_profiles = 1:8, 
                    variances = "varying", 
                    covariances = "zero")
m4
compare_solutions(m4)
fitm4 <- get_fit(m4)
fitm4

# plot(fitm4$AIC, type = "b")
# lines(fitm4$CAIC, type = "b", col = 2)
# lines(fitm4$BIC, type = "b", col = 3)
# lines(fitm4$SABIC, type = "b", col = 4)
# legend(legend = c("AIC", "CAIC", "BIC", "SABIC"), 
#        col = 1:4, lty = 1, x = 7, y = 3800)

fitm4 %>% 
  select(Classes, AIC, CAIC, BIC, SABIC) %>% 
  pivot_longer(col = !Classes,
               names_to = "Index",
               values_to = "value") %>% 
  ggplot(aes(x = Classes, y = value, shape = Index)) +
  geom_line() +
  geom_point(size = 2) +
  theme_few() +
  labs(shape = "Fit Index")


m4.4 <- wr %>% 
  estimate_profiles(n_profiles = 4, 
                  variances = "varying",
                  covariances = "zero") 
m4.4

estm4.4 <- get_estimates(m4.4)
estm4.4
get_data(m4.4)

plot_profiles(m4.4, add_line = TRUE) 

# Latent class proportions
m4.4 %>% 
  get_data() %>% 
  select(CPROB1:CPROB4) %>% 
  colMeans()

# ------------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~          END           ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------------------------------------------------------


sessionInfo()
