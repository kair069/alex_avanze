source("ZI.R")

## reading data
trajan <- as.data.frame(readRDS('trajan_recoded.rds'))

trajan_nblin <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBII)
trajan_nbquad <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBI)
trajan_po_typeA <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "A", family = "poisson")
trajan_po_typeB <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "B", family = "poisson")
trajan_po_typeC <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "C", family = "poisson")
trajan_po_typeD <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "D", family = "poisson")

grid_size <- 1000
mu_grid <- seq(0.01, 100, length = grid_size)

pi0 <- dPO(0, mu = mu_grid)
pit_nblin <- dNBII(0, mu = mu_grid, sigma = exp(trajan_nblin$sigma.coefficients))
pit_nbquad <- dNBI(0, mu = mu_grid, sigma = exp(trajan_nbquad$sigma.coefficients))
pit0_A <- dZI(0, mu = mu_grid, q = plogis(trajan_po_typeA$coef[9]), ZI_type = "A", family = "poisson")
pit0_B <- dZI(0, mu = mu_grid, q = exp(trajan_po_typeB$coef[9]), ZI_type = "B", family = "poisson")
pit0_C <- dZI(0, mu = mu_grid, q = plogis(trajan_po_typeC$coef[9]), ZI_type = "C", family = "poisson")
pit0_D <- dZI(0, mu = mu_grid, q = exp(trajan_po_typeD$coef[9]), ZI_type = "D", family = "poisson")

all_pi <- data.frame(pi0 = pi0,
                     pit0 = c(pit_nblin,pit_nbquad,pit0_A,pit0_B,pit0_C,pit0_D),
                     Model = rep(c("NB-lin","NB-quad","ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic"), each = grid_size))

trajan_p0 <- trajan %>%
  group_by(hormone3, period3) %>%
  summarise(p0 = sum(nshoots == 0)/n(),
            mu = mean(nshoots))

trajan_p0_mu <- trajan_p0
trajan_p0_mu$mu_nblin <- exp(coef(trajan_nblin)) %>% as.numeric
trajan_p0_mu$mu_nbquad <- exp(coef(trajan_nbquad)) %>% as.numeric
trajan_p0_mu$mu_A <- exp(trajan_po_typeA$coef[-9])
trajan_p0_mu$mu_B <- exp(trajan_po_typeB$coef[-9])
trajan_p0_mu$mu_C <- exp(trajan_po_typeC$coef[-9])
trajan_p0_mu$mu_D <- exp(trajan_po_typeD$coef[-9])

obs_pi_pi0 <- tibble(trajan_p0_mu,
                     `NB-lin` = dPO(0, mu_nblin),
                     `NB-quad` = dPO(0, mu_nbquad),
                     `ZI type A: Hurdle` = dPO(0, mu_A),
                     `ZI type B: Poisson-hurdle` = dPO(0, mu_B),
                     `ZI type C: Mixture` = dPO(0, mu_C),
                     `ZI type D: Logistic` = dPO(0, mu_D)) %>%
  dplyr::select(-c(4:10)) %>%
  pivot_longer(cols = 4:9,
               names_to = "Model",
               values_to = "pi0")

obs_pi_pit0 <- tibble(trajan_p0_mu,
                      `NB-lin` = dNBII(0, mu, sigma = exp(trajan_nblin$sigma.coefficients)),
                      `NB-quad` = dNBII(0, mu, sigma = exp(trajan_nbquad$sigma.coefficients)),
                      `ZI type A: Hurdle` = dZI(0, mu, q = plogis(trajan_po_typeA$coef[3]), ZI_type = "A", family = "poisson"),
                      `ZI type B: Poisson-hurdle` = dZI(0, mu, q = exp(trajan_po_typeB$coef[3]), ZI_type = "B", family = "poisson"),
                      `ZI type C: Mixture` = dZI(0, mu, q = plogis(trajan_po_typeC$coef[3]), ZI_type = "C", family = "poisson"),
                      `ZI type D: Logistic` = dZI(0, mu, q = exp(trajan_po_typeD$coef[3]), ZI_type = "D", family = "poisson")) %>%
  dplyr::select(-c(4:10)) %>%
  pivot_longer(cols = 4:9,
               names_to = "Model",
               values_to = "pit0")

obs_pi <- left_join(obs_pi_pi0, obs_pi_pit0)

png("figure2.png", w = 8, h = 5, res = 800, units = "in")
all_pi %>%
  ggplot(aes(x = pi0, y = pit0, col = Model)) +
  theme_bw() +
  geom_line() +
  geom_point(data = obs_pi, aes(y = p0), alpha = .7) +
  ylim(0, 1) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  coord_fixed() +
  xlab(expression(pi[0])) +
  ylab(expression(tilde(pi)[0])) +
  facet_wrap(~ Model) +
  theme(legend.position = "none")
dev.off()
