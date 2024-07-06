source("ZI.R")
library(gamlss.dist)
library(tidyverse)
library(RColorBrewer)

## reading data
trajan <- as.data.frame(readRDS("trajan_recoded.rds"))

# Convenience functios
logit <- function(p) log(p / (1 - p))
expit <- function(x) exp(x) / (1 + exp(x))

# Functions to find pit0 for each type
dTypeA0 <- function(gamma, mu) rep(expit(gamma), length(mu)) # gamma <= 0
dTypeB0 <- function(gamma, mu) exp(-exp(gamma) * mu) # gamma >= 0
dTypeC0 <- function(gamma, mu) (1 - exp(gamma)) + exp(gamma) * exp(-mu) # gamma <= 0
dTypeD0 <- function(gamma, mu) expit(gamma + logit(exp(-mu))) # gamma <-0

rho.P <- function(pit0, mu) {
  (1 - pit0) / (1 - exp(-mu))
}
disp.P <- function(pit0, mu) {
  rho <- (1 - pit0) / (1 - exp(-mu)) # wrt Poisson
  (1 - rho) * mu
} ### in terms of Poisson par and not actual ZI mean

# functions to use in finding parameters given 'target' points specified as pois.p, alt.p
# eg for Poisson P(0) = 0.1, find parameters to deliver altered prob  of 0.2
fNBI <- function(phi) dNBI(0, mu = -log(pois.p), phi) - alt.p
fNBII <- function(phi) dNBII(0, mu = -log(pois.p), phi) - alt.p
tiny <- 1e-6

###### Points at which the graph is constrained to go through
pois.p <- 0.2
alt.p <- 0.4

# now find pars
phi.NBI <- as.numeric(uniroot(fNBI, c(tiny, 10))[1])
# check dNBI(0,mu=-log(pois.p),sigma=phi.NBI)
phi.NBII <- as.numeric(uniroot(fNBII, c(tiny, 50))[1])
# check dNBII(0,mu=-log(pois.p),sigma=phi.NBII)
gammaA <- logit(alt.p) ### changed to logit link.
# check dTypeA0(gammaA,mu=-log(pois.p))
gammaB <- log(-log(alt.p)) - log(-log(pois.p))
# check dTypeB0(gammaB,mu=-log(pois.p))
gammaC <- log(1 - alt.p) - log(1 - pois.p)
# check dTypeC0(gammaC,mu=-log(pois.p))
gammaD <- logit(alt.p) - logit(pois.p)
# check dTypeD0(gammaD,mu=-log(pois.p))

# Set up a range of means
m <- exp(seq(-5, 7, .1))
pois0 <- dPO(0, m)

plot_data1 <- tibble(
  Base_x = pois0,
  `Base (Poisson)` = pois0,
  `A` = dTypeA0(gammaA, m),
  `B` = dTypeB0(gammaB, m),
  `C` = dTypeC0(gammaC, m),
  `D` = dTypeD0(gammaD, m),
  `NB-quad` = dNBI(0, m, phi.NBI),
  `NB-lin` = dNBII(0, m, phi.NBII)
) %>%
  pivot_longer(
    names_to = "Type", values_to = "Prob",
    -Base_x
  ) %>%
  mutate(Type = factor(Type,
    levels = c("Base (Poisson)", "A", "B", "C", "D", "NB-lin", "NB-quad")
  ))

cols <- brewer.pal(7, "Set1")
cols[6] <- "#000000"
linetypes <- c(1, 1, 1, 1, 1, 3, 2)
png("figure3.png", units = "in", res = 800, w = 14, h = 5)
ggplot(
  plot_data1,
  aes(
    x = Base_x,
    y = Prob,
    colour = Type,
    linetype = Type
  )
) +
  scale_colour_manual(
    name = "Type",
    labels = c("Base (Poisson)", "A", "B", "C", "D", "NB-lin", "NB-quad"),
    values = cols
  ) +
  scale_linetype_manual(
    name = "Type",
    labels = c("Base (Poisson)", "A", "B", "C", "D", "NB-lin", "NB-quad"),
    values = linetypes
  ) +
  geom_line(size = 1.0) +
  xlab(expression(pi[0]^P)) +
  ylab(expression(tilde(pi)[0])) +
  # labs(x = "Base probability of 0 (Poisson)",
  #      y = "Altered\nprobability\nof 0") +
  theme_bw() +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 10)) +
  theme(
    axis.title.y =
      element_text(
        angle = 0,
        vjust = 1,
        hjust = 0
      )
  ) +
  theme(legend.title = element_blank()) +
  coord_fixed()
dev.off()
