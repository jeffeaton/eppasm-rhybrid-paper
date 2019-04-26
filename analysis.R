#' # 0. Load libraries
library(tidyverse)
library(gridExtra)
library(grid)
library(here)

#' # 1. Load inputs, data, and model outputs

#' inputs and data
inputs <- readRDS(here("inputs", "eppasm-inputs.rds"))

hhsprev <- readRDS(here("inputs", "hhsprev.rds")) %>%
  mutate(unregion = unregion %>% fct_expand("All") %>%
           fct_relevel("All", "Eastern", "Southern", "Middle", "Western"))

ancsite <- readRDS(here("inputs", "ancsite.rds")) %>%
  mutate(unregion = unregion %>% fct_expand("All") %>%
           fct_relevel("All", "Eastern", "Southern", "Middle", "Western"))

ancrtcens <- readRDS(here("inputs", "ancrtcens.rds")) %>%
  mutate(unregion = unregion %>% fct_expand("All") %>%
           fct_relevel("All", "Eastern", "Southern", "Middle", "Western"))


#' outputs
out <- readRDS(here("fits", "outputs.rds")) %>%
  mutate(unregion = unregion %>% fct_expand("All") %>%
           fct_relevel("All", "Eastern", "Southern", "Middle", "Western"))

param <- readRDS(here("fits", "parameters.rds")) %>%
  ungroup %>%
  mutate(unregion = unregion %>% fct_expand("All") %>%
           fct_relevel("All", "Eastern", "Southern", "Middle", "Western"))

unregion <- sapply(inputs, attr, "unregion")
country <- sapply(inputs, attr, "country")
eppregion <- sapply(inputs, attr, "region")


#' # Table 1: summary of data inputs

table1 <- hhsprev %>%
  group_by(unregion, country) %>%
  summarise(eppregions = n_distinct(eppregion),
            surveys = n_distinct(year)) %>%
  left_join(
    ancsite %>%
    group_by(unregion, country) %>%
    summarise(ancss_sites = n_distinct(site),
              ancss_obs = sum(type == "ancss"),
              ancrt_obs = sum(type == "ancrt"))
  ) %>%
  left_join(
    ancrtcens %>%
    group_by(unregion, country) %>%
    summarise(ancrtcens_years = n_distinct(year))
  ) %>%
  as.data.frame %>%
  bind_rows(
  (.) %>% select(-country) %>% group_by(unregion) %>% summarise_all(sum, na.rm=TRUE),
  (.) %>% select(-unregion, -country) %>% summarise_all(sum, na.rm=TRUE)
  )

table1 %>% write.csv(here("figures", "table1.csv"), row.names = FALSE, na = "")

#' Table 1 brackets
hhsprev %>%
  group_by(unregion, country, eppregion) %>%
  filter(n() > 1) %>%
  group_by(unregion) %>%
  summarise(nreg = n_distinct(paste(country, eppregion)),
            nsurv = n_distinct(paste(country, year)))


#' # Adjust ANC site prevalence data

ancsite <- ancsite %>%
  mutate(pstar = (prev * n + 0.5)/(n + 1),
         W = qnorm(pstar),
         v = 2 * pi * exp(W^2) * pstar * (1 - pstar)/n) %>%
  left_join(param %>%
            filter(modlab == "rhybrid", param == "ancbias") %>%
            select(country, eppregion, ancbias = mean)) %>%
  mutate(prev_adj = pnorm(W - ancbias))



#' # Figure 1: region fits

#' Calculate equilibrium prior

get_max_data_year <- function(obj){
  eppd <- attr(obj, "eppd")
  max(
    eppd$ancsitedat %>% filter(used) %>% .$year,
    eppd$ancrtcens %>% .$year,
    eppd$hhs %>% filter(used) %>% .$year
  )
}

lastyear <- data.frame(country = country,
                       eppregion = eppregion,
                       year = inputs %>% vapply(get_max_data_year, numeric(1)))

equil_rt <- out %>%
  filter(indicator == "prev", modlab == "rspline") %>%
  inner_join(lastyear) %>%
  transmute(unregion,
            country,
            eppregion,
            req = 1 / (11.5 * (1 - mean)),
            log_req = log(req))


create_fig <- function(icountry, ieppreg, nrow = 2, ncol = 2, base_size = 11, fatten = 3){

  df <- filter(out, country == icountry, eppregion == ieppreg, year %in% 1986:2020, modlab %in% c("rspline", "rhybrid"))
  p1dat <- hhsprev %>%
    filter(country == icountry, eppregion == ieppreg) %>%
    mutate(W = qnorm(prev),
           se_W = sqrt(2 * pi * exp(W^2) * se^2),
           cil = pnorm(W - qnorm(0.975)*se_W),
           ciu = pnorm(W + qnorm(0.975)*se_W))

  p1artprop <- df %>%
    filter(indicator == "artprop15to49", modlab == "rhybrid") %>%
    filter(mean > 0, year <= 2017)

  p2ancsite <- ancsite %>% filter(country == icountry, eppregion == ieppreg, used)
  p2anccens <- ancrtcens %>% filter(country == icountry, eppregion == ieppreg)

  th_fig1 <- list(
    scale_fill_brewer(palette = "Set1"),
    scale_color_brewer(palette = "Set1"),
    theme_light(base_size),
    scale_x_continuous(element_blank()),
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(t=4, r=4))
  )

  p1 <- df %>%
    filter(indicator == "prev") %>%
    ggplot(aes(year, mean, ymin = lower, ymax = upper, fill = modlab, color = modlab)) +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_ribbon(data = p1artprop, alpha = 0.3, color = NA) +
    geom_line(data = p1artprop, linetype = "dashed") +
    geom_line() +
    geom_pointrange(aes(year, prev, ymin = cil, ymax = ciu), p1dat, fatten = fatten, inherit.aes = FALSE) +
    ## scale_y_continuous("Prevalence", labels = scales::percent_format(1)) +
    scale_y_continuous(element_blank(), labels = scales::percent_format(1)) +
    th_fig1

  p2 <- df %>%
    filter(indicator == "pregprev") %>%
    ggplot(aes(year, mean, ymin = lower, ymax = upper, fill = modlab, color = modlab)) +
    geom_point(aes(year, prev_adj, group = site, shape = type), p2ancsite, color = "grey", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(aes(year, prev_adj, group = site), p2ancsite, color = "grey", alpha = 0.5, inherit.aes = FALSE) +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_line() +
    geom_point(aes(year, prev), p2anccens, color = "darkgreen", shape = "square", inherit.aes = FALSE) +
    ## scale_y_continuous("Preg. wom. prev.", labels = scales::percent_format(1)) +
    scale_y_continuous(element_blank(), labels = scales::percent_format(1)) +
    th_fig1

  p3 <- df %>%
    filter(indicator == "incid") %>%
    ggplot(aes(year, mean, ymin = lower, ymax = upper, fill = modlab, color = modlab)) +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_line() +
    ## scale_y_continuous("Incidence per 1000", labels = scales::number_format(scale = 1e3)) +
    scale_y_continuous(element_blank(), labels = scales::number_format(scale = 1e3)) +
    th_fig1

  p4 <- df %>%
    filter(indicator == "log r(t)") %>%
    ggplot(aes(year, mean, ymin = lower, ymax = upper, fill = modlab, color = modlab)) +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_line() +
    geom_hline(aes(yintercept = log_req), equil_rt %>% filter(country == icountry, eppregion == ieppreg), linetype = "dashed") +
    ## scale_y_continuous("log r(t)") +
    scale_y_continuous(element_blank()) +
    th_fig1

  tt <- textGrob(paste0(icountry, "\n", sub("SOUTH\\_", "", ieppreg)), rot = 90, gp = gpar(fontsize = base_size+2, font = 2))
  ## grid.arrange(p1, p2, p3, p4, ncol=ncol, nrow=nrow, top = textGrob(paste(icountry, "-", ieppreg), x = 0.02, hjust = 0, gp = gpar(font = 2)))
  grid.arrange(tt, p1, p2, p3, p4, nrow = 1, widths = c(1, 4, 4, 4, 4))
}

## pdf(here("figures", "summary-outputs.pdf"), h=6, w=6.5)
## Map(create_fig, country, eppregion)
## dev.off()


fig1a <- create_fig("Kenya", "Eastern", 1, 4, base_size = 8, fatten = 1)
fig1b <- create_fig("Malawi", "Central Region", 1, 4, base_size = 8, fatten = 1)
fig1c <- create_fig("Ethiopia", "Amhara Urban", 1, 4, base_size = 8, fatten = 1)
fig1d <- create_fig("Mozambique", "SOUTH_Maputo Provincia", 1, 4, base_size = 8, fatten = 1)


tgp <- gpar(fontsize = 10, font = 2)
tx <- 0.58
ty <- 0.07
tj <- "bottom"

tt <- grid.arrange(grob(),
                   textGrob("Prevalence (15-49y)", gp = tgp, x = tx, y = ty, just = tj),
                   textGrob("Pregant women\nprevalence", gp = tgp, x = tx, y = ty, just = tj),
                   textGrob("Incidence per 1000\n(15-49y)", gp = tgp, x = tx, y = ty, just = tj),
                   textGrob("log r(t)", gp = tgp, x = tx, y = ty, just = tj),
                   nrow = 1, widths = c(1, 4, 4, 4, 4))

leg <- df %>%
  ggplot(aes(year, mean, ymin = lower, ymax = upper, fill = modlab, color = modlab)) +
  geom_line() +
  geom_ribbon(alpha = 0.3, color = NA) +
  scale_fill_brewer("model", palette = "Set1",
                    labels = c("rhybrid" = "r-hybrid", "rspline" = "r-spline")) +
  scale_color_brewer("model", palette = "Set1",
                     labels = c("rhybrid" = "r-hybrid", "rspline" = "r-spline")) +
  theme_light(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))



fig1 <- arrangeGrob(tt, fig1a, fig1b, fig1c, fig1d, cowplot::get_legend(leg), nrow = 6, heights = c(1, 4, 4, 4, 4, 1))
ggsave(here("figures", "figure1.png"), fig1, width = 6.8, height = 6.4)


#' For text

out_fig1 <- out %>%
  filter(country == "Kenya" & eppregion == "Eastern" | 
         country == "Malawi" & eppregion == "Central Region" |
         country == "Ethiopia" & eppregion == "Amhara Urban" |
         country == "Mozambique" & eppregion == "SOUTH_Maputo Provincia")

                    
#' # Figure 2: posterior parameters for r-hybrid model

rlog_prior <- data.frame(param = c("rlog.log.r0", "rlog.log.r1", "rlog.log.alpha", "rlog.tmid"),
                         parlabel = factor(1:4, 1:4, c("r[0]", "r[infinity]", "log(alpha)", "t[mid]")),
                         pr_mean = c(log(0.35), log(0.09), log(0.2), 1993),
                         pr_sd = c(0.5, 0.3, 0.5, 5),
                         scale = c(4.3, 5, 3.3, 3.3)) %>%
  mutate(pr_min = pr_mean - scale * pr_sd,
         pr_max = pr_mean + scale * pr_sd,
         pr_q1 = pr_mean - qnorm(0.975) * pr_sd,
         pr_q2 = pr_mean + qnorm(0.975) * pr_sd)

fig2dat <- param %>%
  ungroup %>%
  filter(modlab == "rhybrid", grepl("^rlog", param)) %>%
  bind_rows((.) %>% mutate(unregion = factor("All", levels(.$unregion)))) %>%
  left_join(rlog_prior)
  
fig2 <- fig2dat %>%
  ggplot(aes(unregion, mean, color = unregion, fill = unregion)) +
  geom_hline(aes(yintercept = pr_mean), rlog_prior, linetype = "dashed") +
  geom_blank(aes(y = pr_max), rlog_prior, inherit.aes = FALSE) +
  geom_blank(aes(y = pr_min), rlog_prior, inherit.aes = FALSE) +
  geom_hline(aes(yintercept = pr_q1), rlog_prior, linetype = "dotted") +
  geom_hline(aes(yintercept = pr_q2), rlog_prior, linetype = "dotted") +
  geom_violin(trim = FALSE, color = NA, fill = "grey80") +
  geom_dotplot(binaxis = "y", stackdir ="center", alpha = 0.8, dotsize = 1.0, binpositions = "all") +
  stat_summary(fun.y=mean, aes(ymin=..y.., ymax=..y..), color = "black", geom="errorbar", width = 0.5, size=0.75) +
  facet_wrap(~parlabel, scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(element_blank(), expand = expand_scale()) +
  scale_x_discrete(element_blank()) +
  scale_color_discrete(l = 55) +
  scale_fill_discrete(l = 55) +
  theme_light(10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        strip.text = element_text(color = "black", face = "bold", size = 10, margin = margin(t=0.15, b=0.15, unit = "lines")),
        strip.background = element_rect(fill = "grey90"))

ggsave(here("figures", "figure2.png"), fig2, width = 4, height = 4, units = "in")


#' For text

fig2dat %>%
  group_by(param) %>%o
  summarise(mean(mean) %>% round(2),
            quantile(mean, 0.25) %>% round(2),
            quantile(mean, 0.75) %>% round(2)) %>%
  as.data.frame 

fig2dat %>%
  group_by(unregion, param) %>%
  summarise(mean(mean) %>% round(2),
            quantile(mean, 0.25) %>% round(2),
            quantile(mean, 0.75) %>% round(2)) %>%
  as.data.frame 

#' # Table 2: LOO comparisons

outpred <- readRDS(here("fits", "outpred.rds"))

tab2dat <- outpred %>%
  mutate(unregion = unregion %>% fct_expand("All")) %>%
  bind_rows({.} %>% mutate(unregion = factor("All", levels(.$unregion)))) %>%
  select(unregion, country, eppregion, year, mostrecent, model, crps, elpd_q, qq) %>%
  left_join(
    {.} %>%
    filter(model == "rhybrid") %>%
    select(-model) %>%
    rename(crps_rhyb = crps, elpd_q_rhyb = elpd_q, qq_rhyb = qq)
  )

table2 <- tab2dat %>%  
  group_by(unregion, model) %>%
  summarise(crps_mean = 100*mean(crps),
            crps_diff = 100*mean(crps - crps_rhyb),
            crps_diff_se = 100*sd(crps - crps_rhyb) / sqrt(n()),
            elpd = sum(elpd_q),
            elpd_diff = sum(elpd_q - elpd_q_rhyb),
            elpd_diff_se = sd(elpd_q - elpd_q_rhyb) * sqrt(n()),
            cov80 = mean(qq > 0.1 & qq < 0.9),
            cov95 = mean(qq > 0.025 & qq < 0.975)) %>%
  mutate(crps_mean = sprintf("%.2f", crps_mean),
         crps_diff = sprintf("%.2f (%.2f)", crps_diff, crps_diff_se),
         crps_diff_se = NULL,
         elpd = sprintf("%.2f", elpd),
         elpd_diff = sprintf("%.1f (%.1f)", elpd_diff, elpd_diff_se),
         elpd_diff_se = NULL,
         cov80 = sprintf("%.1f%%", 100*cov80),
         cov95 = sprintf("%.1f%%", 100*cov95))

table3 <- tab2dat %>%
  filter(mostrecent) %>%
  group_by(unregion, model) %>%
  summarise(crps_mean = 100*mean(crps),
            crps_diff = 100*mean(crps - crps_rhyb),
            crps_diff_se = 100*sd(crps - crps_rhyb) / sqrt(n()),
            elpd = sum(elpd_q),
            elpd_diff = sum(elpd_q - elpd_q_rhyb),
            elpd_diff_se = sd(elpd_q - elpd_q_rhyb) * sqrt(n()),
            cov80 = mean(qq > 0.1 & qq < 0.9),
            cov95 = mean(qq > 0.025 & qq < 0.975)) %>%
  mutate(crps_mean = sprintf("%.2f", crps_mean),
         crps_diff = sprintf("%.2f (%.2f)", crps_diff, crps_diff_se),
         crps_diff_se = NULL,
         elpd = sprintf("%.2f", elpd),
         elpd_diff = sprintf("%.1f (%.1f)", elpd_diff, elpd_diff_se),
         elpd_diff_se = NULL,
         cov80 = sprintf("%.1f%%", 100*cov80),
         cov95 = sprintf("%.1f%%", 100*cov95))

table2 %>% write.csv(here("figures", "table2.csv"), row.names = FALSE, na = "")
table3 %>% write.csv(here("figures", "table3.csv"), row.names = FALSE, na = "")

tab2dat %>%
  filter(model == "rhybrid") %>%
  count(unregion)
