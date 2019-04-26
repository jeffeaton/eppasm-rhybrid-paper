setwd("~/Dropbox/Documents/Research/age-specific-incidence/epp-asm-paper/analysis/")

#' # 0. Load packages and model inputs

library(tidyverse)
devtools::load_all("~/Documents/Code/R/eppasm-new-master")

#' # 1. Load model inputs

inputs <- here::here("inputs", "eppasm-inputs.rds") %>% readRDS


#' # 2. Select countries with 2+ surveys and create replicate of
#'      each input with single survey missing.

nsurv <- inputs %>% lapply(attr, "eppd") %>% lapply("[[", "hhs") %>% sapply(nrow)

create_loo <- function(obj){

  eppd <- attr(obj, "eppd")
  
  hhs <- eppd$hhs %>%
    arrange(year) %>%
    mutate(idx = row_number())
    
  val <- hhs %>% select(country, eppregion, year, idx) %>%
    arrange(year, idx) %>%
    mutate(mostrecent = idx == nrow(hhs))

  val$obj <- lapply(val$idx, function(idx_loo) {
    attr(obj, "eppd")$hhs <- hhs[-idx_loo,]
    obj
  })

  val$out <- lapply(val$idx, function(idx_loo){hhs[idx_loo,]})

  val
}

inputs_loo <- lapply(inputs[nsurv > 1], create_loo) %>% bind_rows


#' # 3. Setup DIDE cluster

workdir <- "/Volumes/jeff/eppasm-paper"
dir.create(workdir)

mrc_config <- didehpc::didehpc_config(credentials = "jwe08", workdir=workdir, use_common_lib=FALSE, cluster="mrc", use_workers = TRUE)
pkgsrc <- provisionr::package_sources(local=c("~/Documents/Code/R/eppasm-new-master"),
                                      github=c("jeffeaton/anclik/anclik", "jeffeaton/epp"))
ctx <- context::context_save(file.path(workdir, "context3"), packages=c("eppasm", "scoringRules"), package_sources=pkgsrc)

mrcq <- didehpc::queue_didehpc(ctx, config=mrc_config, initialise=TRUE)

#' # 4. Fit model


#' Full data

mrcq$lapply(inputs, fitmod, name = "fit_rspline", eppmod = "rspline", equil.rprior = TRUE,
            B0=1e5, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rtrend", eppmod = "rtrend",
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid", eppmod = "rhybrid",
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)


#' Full data without pregprev

mrcq$lapply(inputs, fitmod, name = "fit_rspline_nopreg", eppmod = "rspline", equil.rprior = TRUE, pregprev = FALSE,
            B0=1e5, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rtrend_nopreg", eppmod = "rtrend",  pregprev = FALSE,
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid_nopreg", eppmod = "rhybrid", pregprev = FALSE,
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)


#' ## leave-one-out fits

mrcq$lapply(inputs_loo$obj, fitmod, name = "loo_rspline", eppmod="rspline", equil.rprior=TRUE,
            B0=1e5, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs_loo$obj, fitmod, name = "loo_rhybrid", eppmod="rhybrid",
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)
mrcq$lapply(inputs_loo$obj, fitmod, name = "loo_rtrend", eppmod="rtrend",
            B0=1e4, B=1e3, number_k = 500, opt_iter = 1:5*5)


#' # 3. Simulate outputs

rspline <- mrcq$task_bundle_get("fit_rspline")$results()
rtrend <- mrcq$task_bundle_get("fit_rtrend")$results()
rhybrid <- mrcq$task_bundle_get("fit_rhybrid")$results()

rspline_nopreg <- mrcq$task_bundle_get("fit_rspline_nopreg")$results()
rtrend_nopreg <- mrcq$task_bundle_get("fit_rtrend_nopreg")$results()
rhybrid_nopreg <- mrcq$task_bundle_get("fit_rhybrid_nopreg")$results()

#' Extend projection for r-hybrid model
rhybrid <- lapply(rhybrid, extend_projection, 52L)
rhybrid_nopreg <- lapply(rhybrid_nopreg, extend_projection, 52L)

#' Simulate outputs
unregion <- sapply(inputs, attr, "unregion")
country <- sapply(inputs, attr, "country")
eppregion <- sapply(inputs, attr, "region")

mrcq$mapply(tidy_output, rspline, "rspline", country, eppregion, FALSE, name = "out_rspline")
mrcq$mapply(tidy_output, rtrend, "rtrend", country, eppregion, FALSE, name = "out_rtrend")
mrcq$mapply(tidy_output, rhybrid, "rhybrid", country, eppregion, FALSE, name = "out_rhybrid")

mrcq$mapply(tidy_output, rspline_nopreg, "rspline: no preg", country, eppregion, FALSE, name = "out_rspline_nopreg")
mrcq$mapply(tidy_output, rtrend_nopreg, "rtrend: no preg", country, eppregion, FALSE, name = "out_rtrend_nopreg")
mrcq$mapply(tidy_output, rhybrid_nopreg, "rhybrid: no preg", country, eppregion, FALSE, name = "out_rhybrid_nopreg")


mrcq$lapply(rspline, get_param, name = "param_rspline")
mrcq$lapply(rtrend,  get_param, name = "param_rtrend")
mrcq$lapply(rhybrid, get_param, name = "param_rhybrid")

mrcq$lapply(rspline_nopreg, get_param, name = "param_rspline_nopreg")
mrcq$lapply(rtrend_nopreg,  get_param, name = "param_rtrend_nopreg")
mrcq$lapply(rhybrid_nopreg, get_param, name = "param_rhybrid_nopreg")


out <- c(mrcq$task_bundle_get("out_rspline")$results(),
         mrcq$task_bundle_get("out_rtrend")$results(),
         mrcq$task_bundle_get("out_rhybrid")$results(),
         mrcq$task_bundle_get("out_rspline_nopreg")$results(),
         mrcq$task_bundle_get("out_rtrend_nopreg")$results(),
         mrcq$task_bundle_get("out_rhybrid_nopreg")$results()) %>%
  lapply(
    function(x)
      bind_rows(x$core,
                x$pregprev %>%
                filter(agegr == "15-49") %>%
                select(country:year, mean:upper) %>%
                mutate(indicator = "pregprev"))
  ) %>%
  bind_rows %>%
  filter(year %in% 1980:2020) %>%
  left_join(data.frame(unregion, country) %>% unique) %>%
  select(unregion, everything())


fn <- function(set, modlab){
  Map(data.frame,
      unregion = unregion,
      country = country,
      eppregion = eppregion,
      modlab = modlab,
      sampleid = set %>% lapply(nrow) %>% lapply(seq_len),
      set) %>%
    bind_rows %>%
    gather(param, value, -(unregion:sampleid))
}

param <- bind_rows(
  mrcq$task_bundle_get("param_rspline")$results() %>% fn("rspline"),
  mrcq$task_bundle_get("param_rtrend")$results() %>% fn("rtrend"),
  mrcq$task_bundle_get("param_rhybrid")$results() %>% fn("rhybrid"),
  mrcq$task_bundle_get("param_rspline_nopreg")$results() %>% fn("rspline: no preg"),
  mrcq$task_bundle_get("param_rtrend_nopreg")$results() %>% fn("rtrend: no preg"),
  mrcq$task_bundle_get("param_rhybrid_nopreg")$results() %>% fn("rhybrid: no preg")
)

param <- param %>%
  group_by(unregion, country, eppregion, modlab, param) %>%
  summarise(mean = mean(value, na.rm=TRUE),
            se = sd(value, na.rm=TRUE),
            median = median(value, na.rm=TRUE),
            lower = quantile(value, 0.025, na.rm=TRUE),
            upper = quantile(value, 0.975, na.rm=TRUE))

            
saveRDS(out, here::here("fits", "outputs.rds"))
saveRDS(param, here::here("fits", "parameters.rds"))


#' # 4. Simulate out of sample predictions for LOO fits

loo_rspline <- mrcq$task_bundle_get("loo_rspline")$results()
loo_rtrend <- mrcq$task_bundle_get("loo_rtrend")$results()
loo_rhybrid <- mrcq$task_bundle_get("loo_rhybrid")$results()

#' Extend r-hybrid simulation
loo_rhybrid <- lapply(loo_rhybrid, extend_projection, 52L)

#' Simulate likelihood for withheld survey prevalence

inputs_loo <- inputs_loo %>%
  left_join(data.frame(unregion = unregion, country) %>% unique) %>%
  select(unregion, everything()) %>%
  mutate(label = paste(country, "-", eppregion))

mrcq$mapply(outpred_prev, rspline[inputs_loo$label], inputs_loo$out, name = "inpred_rspline")
mrcq$mapply(outpred_prev, rtrend[inputs_loo$label], inputs_loo$out, name = "inpred_rtrend")
mrcq$mapply(outpred_prev, rhybrid[inputs_loo$label], inputs_loo$out, name = "inpred_rhybrid")

mrcq$mapply(outpred_prev, loo_rspline, inputs_loo$out, name="outpred_rspline")
mrcq$mapply(outpred_prev, loo_rtrend, inputs_loo$out, name="outpred_rtrend")
mrcq$mapply(outpred_prev, loo_rhybrid, inputs_loo$out, name="outpred_rhybrid")

outpred <- bind_rows(
  mrcq$task_bundle_get("outpred_rspline")$results() %>% bind_rows() %>% mutate(model = "rspline"),
  mrcq$task_bundle_get("outpred_rtrend")$results() %>% bind_rows() %>% mutate(model = "rtrend"),
  mrcq$task_bundle_get("outpred_rhybrid")$results() %>% bind_rows() %>% mutate(model = "rhybrid")
) %>%
  right_join(inputs_loo %>% select(unregion, country, eppregion, year, label, mostrecent), .)

inpred <- bind_rows(
  mrcq$task_bundle_get("inpred_rspline")$results() %>% bind_rows() %>% mutate(model = "rspline"),
  mrcq$task_bundle_get("inpred_rtrend")$results() %>% bind_rows() %>% mutate(model = "rtrend"),
  mrcq$task_bundle_get("inpred_rhybrid")$results() %>% bind_rows() %>% mutate(model = "rhybrid")
) %>%
  right_join(inputs_loo %>% select(unregion, country, eppregion, year, label, mostrecent), .)

saveRDS(outpred, here::here("fits", "outpred.rds"))
saveRDS(inpred, here::here("fits", "inpred.rds"))
