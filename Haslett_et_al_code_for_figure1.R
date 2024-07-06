source("ZI.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

library(gridExtra)
library(cowplot)
library(ggpubr)

## reading data
trajan <- as.data.frame(readRDS('trajan_recoded.rds'))

trajan_po <- glm(nshoots ~ 0 + period3 : hormone3, data = trajan, family = poisson)
trajan_nblin <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBII)
trajan_nbquad <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBI)
trajan_po_typeA <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "A", family = "poisson")
trajan_po_typeB <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "B", family = "poisson")
trajan_po_typeC <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "C", family = "poisson")
trajan_po_typeD <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "D", family = "poisson")

trajan_summary <- trajan %>%
  group_by(hormone3, period3) %>%
  summarise(p0_obs = sum(nshoots == 0)/n(),
            mean_obs = mean(nshoots)) %>%
  ungroup()

trajan_estimates <- trajan_summary %>%
  mutate(mean_A_int = predict(trajan_po_typeA,
                              newdata = trajan_summary[,1:2],
                              what = "mu", type = "response") %>%
           get_mean(mu = .,
                    q = coef(trajan_po_typeA, what = "zero") %>% plogis,
                    ZI_type = "A"),
         mean_B_int = predict(trajan_po_typeB,
                              newdata = trajan_summary[,1:2],
                              what = "mu", type = "response") %>%
           get_mean(mu = .,
                    q = coef(trajan_po_typeB, what = "zero") %>% exp,
                    ZI_type = "B"),
         mean_C_int = predict(trajan_po_typeC,
                              newdata = trajan_summary[,1:2],
                              what = "mu", type = "response") %>%
           get_mean(mu = .,
                    q = coef(trajan_po_typeC, what = "zero") %>% plogis,
                    ZI_type = "C"),
         mean_D_int = predict(trajan_po_typeD,
                              newdata = trajan_summary[,1:2],
                              what = "mu", type = "response") %>%
           get_mean(mu = .,
                    q = coef(trajan_po_typeD, what = "zero") %>% exp,
                    ZI_type = "D"),
         p0_A_int = coef(trajan_po_typeA, what = "zero") %>% plogis %>% rep(8) %>%
           dZI(x = 0, mu = predict(trajan_po_typeA,
                                   newdata = trajan_summary[,1:2],
                                   what = "mu", type = "response"),
               q = ., ZI_type = "A", family = "poisson"),
         p0_B_int = coef(trajan_po_typeB, what = "zero") %>% exp %>% rep(8) %>%
           dZI(x = 0, mu = predict(trajan_po_typeB,
                                   newdata = trajan_summary[,1:2],
                                   what = "mu", type = "response"),
               q = ., ZI_type = "B", family = "poisson"),
         p0_C_int = coef(trajan_po_typeC, what = "zero") %>% plogis %>% rep(8) %>%
           dZI(x = 0, mu = predict(trajan_po_typeC,
                                   newdata = trajan_summary[,1:2],
                                   what = "mu", type = "response"),
               q = ., ZI_type = "C", family = "poisson"),
         p0_D_int = coef(trajan_po_typeD, what = "zero") %>% exp %>% rep(8) %>%
           dZI(x = 0, mu = predict(trajan_po_typeD,
                                   newdata = trajan_summary[,1:2],
                                   what = "mu", type = "response"),
               q = ., ZI_type = "D", family = "poisson"))

trajan_estimates <- trajan_estimates %>%
  mutate(mean_po_int = predict(trajan_po,
                               newdata = trajan_summary[,1:2],
                               type = "response"),
         mean_nblin_int = predict(trajan_nblin,
                                  newdata = trajan_summary[,1:2],
                                  type = "response"),
         mean_nbquad_int = predict(trajan_nbquad,
                                  newdata = trajan_summary[,1:2],
                                  type = "response"),
         p0_po_int = dPO(0, mu = mean_po_int),
         p0_nblin_int = dNBII(0, mu = mean_nblin_int,
                              sigma = predict(trajan_nblin,
                                              what = "sigma",
                                              type = "response")[1]),
         p0_nbquad_int = dNBI(0, mu = mean_nbquad_int,
                             sigma = predict(trajan_nbquad,
                                             what = "sigma",
                                             type = "response")[1]))

trajan_estimates$var_obs <- trajan %>%
  group_by(hormone3, period3) %>%
  summarise(var_obs = var(nshoots)) %>%
  pull(var_obs) %>%
  as.numeric

get_var <- Vectorize(function(mu, q, ZI_type) {
  ex <- get_mean(mu = mu, q = q, ZI_type = ZI_type)
  ex2 <- sum((0:200)^2 * dZI(0:200, mu = mu, q = q, ZI_type = ZI_type, family = "poisson"))
  return(ex2 - ex^2)
}, "mu")

trajan_estimates <- trajan_estimates %>%
  mutate(var_po_int = mean_po_int,
         var_nblin_int = mean_nblin_int + mean_nblin_int^2 * exp(predict(trajan_nblin, "sigma")[1]),
         var_nbquad_int = mean_nbquad_int + mean_nbquad_int^2 * exp(predict(trajan_nbquad, "sigma")[1]),
         var_A_int = get_var(mu = mean_A_int, q = trajan_po_typeA$coef[9] %>% plogis, ZI_type = "A"),
         var_B_int = get_var(mu = mean_B_int, q = trajan_po_typeB$coef[9] %>% exp, ZI_type = "B"),
         var_C_int = get_var(mu = mean_C_int, q = trajan_po_typeC$coef[9] %>% plogis, ZI_type = "C"),#(1 - p0_C_int) * mean_C_int * (1 + mean_C_int * p0_C_int),
         var_D_int = get_var(mu = mean_D_int, q = trajan_po_typeD$coef[9] %>% exp, ZI_type = "D"))

trajan_zero <- trajan_estimates[,grep("p0", names(trajan_estimates))]
trajan_mean <- trajan_estimates[,grep("mean", names(trajan_estimates))]
trajan_var <- trajan_estimates[,grep("var", names(trajan_estimates))]

trajan_zero <- trajan_zero %>%
  pivot_longer(2:8,
               names_to = "model",
               values_to = "pit0")
trajan_zero$model <- as.factor(trajan_zero$model)
levels(trajan_zero$model) <- c("ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic","NB-lin","NB-quad","Poisson")
trajan_zero$model <- factor(trajan_zero$model, levels = c("Poisson","NB-lin","NB-quad","ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic"))

trajan_mean <- trajan_mean %>%
  pivot_longer(2:8,
               names_to = "model",
               values_to = "fitted_mean")
trajan_mean$model <- as.factor(trajan_mean$model)
levels(trajan_mean$model) <- c("ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic","NB-lin","NB-quad","Poisson")
trajan_mean$model <- factor(trajan_mean$model, levels = c("Poisson","NB-lin","NB-quad","ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic"))

trajan_var <- trajan_var %>%
  pivot_longer(2:8,
               names_to = "model",
               values_to = "fitted_var")
trajan_var$model <- as.factor(trajan_var$model)
levels(trajan_var$model) <- c("ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic","NB-lin","NB-quad","Poisson")
trajan_var$model <- factor(trajan_var$model, levels = c("Poisson","NB-lin","NB-quad","ZI type A: Hurdle","ZI type B: Poisson-hurdle","ZI type C: Mixture","ZI type D: Logistic"))

part_a1 <- trajan_zero %>%
  filter(model %in% c("Poisson","NB-lin","NB-quad")) %>%
  ggplot(aes(x = p0_obs, y = pit0)) +
  theme_bw() +
  geom_point(aes(colour = model), alpha = .7) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  ylab("Fitted probability of zero") +
  xlab("Observed proportion of zeros") +
  coord_fixed() + 
  facet_wrap(~ model, nrow = 1) +
  theme(panel.spacing.x = unit(4, "mm")) +
  ggtitle("(a)") +
  xlim(0,1) + ylim(0,1) +
  scale_color_manual(values = c("black",gg_color_hue(6)[1:2])) +
  labs(color = "Model") +
  theme(legend.position = "none")

part_a2 <- trajan_mean %>%
  filter(model %in% c("Poisson","NB-lin","NB-quad")) %>%
  ggplot(aes(x = mean_obs, y = fitted_mean)) +
  theme_bw() +
  geom_point(aes(colour = model), alpha = .7) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  ylab("Fitted means") +
  xlab("Observed means") +
  coord_fixed() + 
  facet_wrap(~ model, nrow = 1) +
  theme(panel.spacing.x = unit(4, "mm")) +
  ggtitle("") +
  xlim(0,8.5) + ylim(0,8.5) +
  scale_color_manual(values = c("black",gg_color_hue(6)[1:2])) +
  labs(color = "Model") +
  theme(legend.position = "none")

part_b1 <- trajan_zero %>%
  filter(!(model %in% c("Poisson","NB-lin","NB-quad"))) %>%
  ggplot(aes(x = p0_obs, y = pit0)) +
  theme_bw() +
  geom_point(aes(colour = model), alpha = .7) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  ylab("Fitted probability of zero") +
  xlab("Observed proportion of zeros") +
  coord_fixed() + 
  facet_wrap(~ model, nrow = 1) +
  theme(panel.spacing.x = unit(4, "mm")) +
  ggtitle("(b)") +
  xlim(0,1) + ylim(0,1) +
  scale_color_manual(values = gg_color_hue(6)[3:6]) +
  labs(color = "Model") +
  theme(legend.position = "none")

part_b2 <- trajan_mean %>%
  filter(!(model %in% c("Poisson","NB-lin","NB-quad"))) %>%
  ggplot(aes(x = mean_obs, y = fitted_mean)) +
  theme_bw() +
  geom_point(aes(colour = model), alpha = .7) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  ylab("Fitted means") +
  xlab("Observed means") +
  coord_fixed() + 
  facet_wrap(~ model, nrow = 1) +
  theme(panel.spacing.x = unit(4, "mm")) +
  ggtitle("") +
  xlim(0,8.5) + ylim(0,8.5) +
  scale_color_manual(values = gg_color_hue(6)[3:6]) +
  labs(color = "Model") +
  theme(legend.position = "none")

p1 <- ggarrange(part_a1, part_a2, nrow = 2, heights = c(4,4), align = "v")
p2 <- ggarrange(part_b1, part_b2, nrow = 2, heights = c(4,4), align = "v")

png("figure1.png", units = "in", res = 800, w = 12, h = 5)
ggarrange(p1, p2, ncol = 2, widths = c(6,8))
dev.off()