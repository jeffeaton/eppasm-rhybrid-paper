setwd("~/Dropbox/Documents/Research/age-specific-incidence/epp-asm-paper/analysis-revision1/")

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

#' Full data, various knot spacings
#' 
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid_dk1", eppmod = "rhybrid", rw_dk = 1,
            B0=1e4, B=1e3, number_k = 1000, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid_dk2", eppmod = "rhybrid", rw_dk = 2,
            B0=1e4, B=1e3, number_k = 1000, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid_dk3", eppmod = "rhybrid", rw_dk = 3,
            B0=1e4, B=1e3, number_k = 1000, opt_iter = 1:5*5)
mrcq$lapply(inputs, fitmod, name = "fit_rhybrid_dk4", eppmod = "rhybrid", rw_dk = 4,
            B0=1e4, B=1e3, number_k = 1000, opt_iter = 1:5*5)


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

rhybrid_dk1 <- mrcq$task_bundle_get("fit_rhybrid_dk1")$results()
rhybrid_dk2 <- mrcq$task_bundle_get("fit_rhybrid_dk2")$results()
rhybrid_dk3 <- mrcq$task_bundle_get("fit_rhybrid_dk3")$results()
rhybrid_dk4 <- mrcq$task_bundle_get("fit_rhybrid_dk4")$results()


rspline_nopreg <- mrcq$task_bundle_get("fit_rspline_nopreg")$results()
rtrend_nopreg <- mrcq$task_bundle_get("fit_rtrend_nopreg")$results()
rhybrid_nopreg <- mrcq$task_bundle_get("fit_rhybrid_nopreg")$results()

#' Extend projection for r-hybrid model
rhybrid <- lapply(rhybrid, extend_projection, 52L)
rhybrid_nopreg <- lapply(rhybrid_nopreg, extend_projection, 52L)

rhybrid_dk1 <- lapply(rhybrid_dk1, extend_projection, 52L)
rhybrid_dk2 <- lapply(rhybrid_dk2, extend_projection, 52L)
rhybrid_dk3 <- lapply(rhybrid_dk3, extend_projection, 52L)
rhybrid_dk4 <- lapply(rhybrid_dk4, extend_projection, 52L)

#' Simulate outputs
unregion <- sapply(inputs, attr, "unregion")
country <- sapply(inputs, attr, "country")
eppregion <- sapply(inputs, attr, "region")

mrcq$mapply(tidy_output, rspline, "rspline", country, eppregion, FALSE, name = "out_rspline",)
mrcq$mapply(tidy_output, rtrend, "rtrend", country, eppregion, FALSE, name = "out_rtrend")
mrcq$mapply(tidy_output, rhybrid, "rhybrid", country, eppregion, FALSE, name = "out_rhybrid")

mrcq$mapply(tidy_output, rspline_nopreg, "rspline: no preg", country, eppregion, FALSE, name = "out_rspline_nopreg")
mrcq$mapply(tidy_output, rtrend_nopreg, "rtrend: no preg", country, eppregion, FALSE, name = "out_rtrend_nopreg")
mrcq$mapply(tidy_output, rhybrid_nopreg, "rhybrid: no preg", country, eppregion, FALSE, name = "out_rhybrid_nopreg")

mrcq$mapply(tidy_output, rhybrid_dk1, "rhybrid: dk=1", country, eppregion, FALSE, name = "out_rhybrid_dk1")
mrcq$mapply(tidy_output, rhybrid_dk2, "rhybrid: dk=2", country, eppregion, FALSE, name = "out_rhybrid_dk2")
mrcq$mapply(tidy_output, rhybrid_dk3, "rhybrid: dk=3", country, eppregion, FALSE, name = "out_rhybrid_dk3")
mrcq$mapply(tidy_output, rhybrid_dk4, "rhybrid: dk=4", country, eppregion, FALSE, name = "out_rhybrid_dk4")


mrcq$lapply(rspline, get_param, name = "param_rspline")
mrcq$lapply(rtrend,  get_param, name = "param_rtrend")
mrcq$lapply(rhybrid, get_param, name = "param_rhybrid")

mrcq$lapply(rspline_nopreg, get_param, name = "param_rspline_nopreg")
mrcq$lapply(rtrend_nopreg,  get_param, name = "param_rtrend_nopreg")
mrcq$lapply(rhybrid_nopreg, get_param, name = "param_rhybrid_nopreg")


out <- c(mrcq$task_bundle_get("out_rspline")$results(),
         mrcq$task_bundle_get("out_rtrend")$results(),
         mrcq$task_bundle_get("out_rhybrid")$results(),
         mrcq$task_bundle_get("out_rhybrid_dk4")$results(),
         mrcq$task_bundle_get("out_rhybrid_dk3")$results(),
         mrcq$task_bundle_get("out_rhybrid_dk2")$results(),
         mrcq$task_bundle_get("out_rhybrid_dk1")$results(),
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
  filter(year %in% 1980:2021) %>%
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


#' number of IMIS iterations

imis_iter <- list("rhybrid: dk=5" = rhybrid,
     "rhybrid: dk=4" = rhybrid_dk4,
     "rhybrid: dk=3" = rhybrid_dk3,
     "rhybrid: dk=2" = rhybrid_dk2,
     "rhybrid: dk=1" = rhybrid_dk1,
     "rspline" = rspline,
     "rtrend" = rtrend) %>%
  lapply(lapply, "[[", "stat") %>%
  lapply(sapply, nrow) %>%
  Map(data.frame, modlab = names(.), region = lapply(., names), imis_iter = .) %>%
  bind_rows()
     
     
saveRDS(out, here::here("fits", "outputs.rds"))
saveRDS(param, here::here("fits", "parameters.rds"))
saveRDS(imis_iter, here::here("fits", "imis_iter.rds"))



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


#' # 5. Malawi Central outputs for appendix

fit <- rhybrid[["Malawi - Central Region"]]
fit <- extend_projection(fit, 52L)

ss <- fit$fp$ss

year <- fit$fp$ss$proj_start + 1:fit$fp$ss$PROJ_YEARS - 1L


## simulate model projections
param_list <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

fp_list <- lapply(param_list, function(par) update(fit$fp, list=par))
mod_list <- lapply(fp_list, simmod)

add_mod_names <- function(mod) {
  dimnames(mod) <- list(age = 15:80,
                        sex = c("male", "female"),
                        hivstatus = c("negative", "positive"),
                        year = year)
  
  dimnames(attr(mod, "infections")) <- dimnames(mod)[c(1,2,4)]
  dimnames(attr(mod, "hivdeaths")) <- dimnames(mod)[c(1,2,4)]
  mod
}

merge_mod <- function(mod) {
  as.data.frame.table(apply(mod, c(1:2, 4), sum),
                      responseName = "totpop") %>%
    left_join(as.data.frame.table(mod[,,2,], responseName = "hivpop"),
              by = c("age", "sex", "year")) %>%
    left_join(as.data.frame.table(attr(mod, "infections"), responseName = "infections"),
              by = c("age", "sex", "year")) %>%
    left_join(as.data.frame.table(attr(mod, "hivdeaths"), responseName = "hivdeaths"),
              by = c("age", "sex", "year")) %>%
  type.convert() %>%
    left_join(
    {.} %>%
    transmute(age, sex, year = year + 1, totpop_last = totpop, suscpop_last = totpop - hivpop),
    by = c("age", "sex", "year")
    )
}


calc_agecat_output <- function(mod) {

  mod %>%
    add_mod_names() %>%
    merge_mod() %>%
    mutate(agecat = cut(age, c(15, 25, 35, 50, Inf), c("15-24", "25-34", "35-49", "50+"), TRUE, FALSE)) %>%
    group_by(sex, agecat, year) %>%
    summarise_at(vars(totpop:suscpop_last), sum) %>%
    group_by(sex, agecat, year) %>%
    transmute(prev = hivpop / totpop,
              incid = infections / suscpop_last,
              aidsmx = hivdeaths / suscpop_last)
}

calc_age1_output <- function(mod, year = seq(1985, 2020, 5)) {

  mod %>%
    add_mod_names() %>%
    merge_mod() %>%
    filter(year %in% !!year)%>%
    group_by(sex, age, year) %>%
    summarise_at(vars(totpop:suscpop_last), sum) %>%
    group_by(sex, age, year) %>%
    transmute(prev = hivpop / totpop,
              incid = infections / suscpop_last,
              aidsmx = hivdeaths / suscpop_last)
}


system.time(out_agecat <- parallel::mclapply(mod_list, calc_agecat_output, mc.cores = 10))

mwi_central_agecat_out <- out_agecat %>%
  bind_rows() %>%
  gather(indicator, value, prev, incid, aidsmx) %>%
  group_by(sex, agecat, year, indicator) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            lower = quantile(value, 0.025, na.rm = TRUE),
            upper = quantile(value, 0.975, na.rm = TRUE))


system.time(out_age1 <- parallel::mclapply(mod_list, calc_age1_output, mc.cores = 10))

mwi_central_age1_out <- out_age1 %>%
  bind_rows() %>%
  gather(indicator, value, prev, incid, aidsmx) %>%
  group_by(sex, age, year, indicator) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            lower = quantile(value, 0.025, na.rm = TRUE),
            upper = quantile(value, 0.975, na.rm = TRUE))



saveRDS(mwi_central_agecat_out, here::here("fits", "mwi_central_agecat_out.rds"))
saveRDS(mwi_central_age1_out, here::here("fits", "mwi_central_age1_out.rds"))
