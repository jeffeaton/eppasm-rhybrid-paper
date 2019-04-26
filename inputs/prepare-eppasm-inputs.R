
#' # 0. Load packages

library(tidyverse)
devtools::load_all("~/Documents/Code/R/eppasm-new-master")

options(mc.cores = parallel::detectCores())


#' # 1. Prepare inputs from 2018 Spectrum files
#'
#' * Exclude due to concentrated epidemic structure: Cape Verde, Comores, Madagascar, Mauritania, Mauritius, Niger, Senegal
#' * Exclude due to not using EPP: South Africa
#' * Exclude due to no reliable household sureys: Nigeria (37 states), South Sudan

exclude_str <- c("(^Cabo\\_Verde)", "(^Comores)", "(^Madagascar)", "(^Mauritanie)",
                 "(^Mauritius)", "(^Niger\\_)", "(^Nigeria\\_)", "(Senegal)", "(South Africa)",
                 "(^South\\_Sudan)", "(^Togo\\_2018)") %>%
  paste0(collapse = "|")

pjnz <- list.files("~/Documents/Data/Spectrum files/2018 final/SSA", "PJNZ$", full.names=TRUE, ignore.case = TRUE)
pjnz <- pjnz[grep(exclude_str, basename(pjnz), invert=TRUE)]


#' Read inputs
inputs <- parallel::mclapply(pjnz, prepare_spec_fit, proj.end=2021.5)


## Handle Ethiopia region names + U/R

et_idx <- sapply(inputs, attr, "country") == "Ethiopia"
et_reg <- Map(paste, lapply(inputs[et_idx], attr, "region"),
              lapply(inputs[et_idx], lapply, attr, "region"))
et_reg <- lapply(et_reg, sub, pattern = " total", replacement = "", ignore.case = TRUE)

et_in <- Map(Map, "attr<-", inputs[et_idx], "region", et_reg)
et_in <- Map("attr<-", et_in, "country", "Ethiopia")
et_in <- Map("attr<-", et_in, "region", lapply(inputs[et_idx], attr, "region"))

inputs[et_idx] <- et_in

## Equatorial Guinea is a single national EPP fit: assign region name as country name

attr(inputs[sapply(inputs, attr, "country") == "Equatorial Guinea"][[1]], "region") <- "Equatorial Guinea"


## If a subnational spectrum file with 1 EPP region, update region name to match spectrum file
## Also trim trailing spaces
which_single <- sapply(inputs, length) == 1
inputs[which_single] <- lapply(inputs[which_single], function(x){ attr(x[[1]], "region") <- gsub("\\s+$", "", attr(x, "region")); x})

## Ensure the EPP country name matches spectrum
table(sapply(lapply(inputs, "[[", 1), attr, "country") == sapply(inputs, attr, "country"))
inputs <- lapply(inputs, function(x){lapply(x, "attr<-", "country", attr(x, "country"))})

## Unlist regions
inputs <- unlist(inputs, FALSE)

names(inputs) <- paste0(sapply(inputs, attr, "country"), " - ", sapply(inputs, attr, "region"))


#' # 2. Assign UN regions

unregion  <- bind_rows(
  data.frame(unregion = "Western",
             country = c("Benin", "Burkina Faso", "Cabo Verde", "Côte d'Ivoire",
                         "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia",
                         "Mali", "Mauritania", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo")),
  data.frame(unregion = "Southern",
             country = c("Botswana", "Swaziland", "Lesotho", "Namibia", "South Africa")),
  data.frame(unregion = "Middle",
             country = c("Angola", "Cameroon", "Central African Republic", "Chad", "Democratic Republic of the Congo",
                         "Equatorial Guinea", "Gabon", "Congo", "São Tomé and Príncipe")),
  data.frame(unregion = "Eastern",
             country = c("Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar",
                         "Malawi", "Mauritius", "Mayotte", "Mozambique", "Réunion", "Rwanda", "Seychelles",
                         "Somalia", "South Sudan", "United Republic of Tanzania", "Uganda", "Zambia", "Zimbabwe"))
)

inputs  <- inputs %>% Map("attr<-", ., "unregion", unregion$unregion[match(sapply(., attr, "country"), unregion$country)])


#' # 2. Review prevalence data

unregion <- sapply(inputs, attr, "unregion")
country <- sapply(inputs, attr, "country")
eppregion <- sapply(inputs, attr, "region")

hhsprev <- inputs %>%
  lapply(attr, "eppd") %>%
  lapply("[[", "hhs") %>%
  Map(data.frame, unregion = unregion, country = country, eppregion = eppregion, .) %>%
  bind_rows

#' Check 'not used' surveys
#' * Equatorial Guinea 2008: not sure why not used, use this
#' * Congo Rural: this was an urban only survey -- not sure what this data are (don't use)
#' * Tanzania: update this with re-analysed survey data
#' * Kenya Northeastern: drop these
#' * Rwanda: Not sure why they haven't used, but don't use

hhsprev %>% filter(!used)
hhsprev %>% filter(country == "Equatorial Guinea")

hhsprev <- hhsprev %>%
  mutate(used = if_else(country == "Equatorial Guinea", TRUE, used)) %>%
  filter(used)

#' Update Tanzania and Zambia surveys to re-analysed data using current boundaries
#' For Zambia 2013-14, 

hhsprev <- hhsprev %>%
  filter(!country %in% c("United Republic of Tanzania", "Zambia"),
         !(country == "Uganda" & year == "2011")) %>%
  bind_rows(
    read.csv(here::here("inputs", "tanzania-region-hiv-prevalence.csv")) %>%
    transmute(unregion = "Eastern",
              country = "United Republic of Tanzania",
              eppregion = fct_recode(region, "Dar es Salaam" = "Dar Es Salaam", "Dodmoa" = "Dodoma", "Shinyaga" = "Shinyanga"),
              year,
              sex = "both",
              agegr = "15-49",
              n = counts,
              prev = prev,
              se = se,
              deff = 2,
              deff_approx = 2,
              used = TRUE),
    read.csv(here::here("inputs", "zambia-15to49-prev.csv")) %>%
    transmute(unregion = "Eastern", country, eppregion, year, sex = "both", agegr = "15-49", n, prev, se, deff, deff_approx = 2, used = TRUE)
  )

inputs <- inputs %>%
  lapply(function(x){
    attr(x, "eppd")$hhs <- filter(hhsprev, country == attr(x, "country"), eppregion == attr(x, "region"), used)
    if(nrow(attr(x, "eppd")$hhs) < 2) warning(attr(x, "country"), " - ", attr(x, "region"), " found fewer than 2 HHS prevalence")
    x
  })

#' Remove Namibia ANC-RT census prevalence because it is national data in each region
#' Tanzania - Geita has an errant input (n but no prevalence), remove

idx <- sapply(inputs, attr, "country") %in% c("Namibia", "United Republic of Tanzania")
inputs[idx] <- lapply(inputs[idx], function(x){attr(x, "eppd")$ancrtcens <- NULL; x})

#' Remove pseudo sites
inputs <- inputs %>%
  lapply(function(x){
    attr(x, "eppd")$ancsitedat <- attr(x, "eppd")$ancsitedat %>% filter(!grepl("(pseudo)|(^Garissa 2$)", site, ignore.case = TRUE))
    x
  })
  

#' Extract ANC data

ancrtcens <- inputs %>%
  lapply(attr, "eppd") %>%
  lapply("[[", "ancrtcens")

ancrtcens[sapply(ancrtcens, is.null)] <- replicate(sum(sapply(ancrtcens, is.null)), data.frame())
has_ancrtcens <- sapply(ancrtcens, nrow) > 0

ancrtcens <- Map(data.frame,
                 unregion = unregion[has_ancrtcens],
                 country = country[has_ancrtcens],
                 eppregion = eppregion[has_ancrtcens],
                 ancrtcens[has_ancrtcens]) %>%
     bind_rows()


ancsite <- inputs %>%
  lapply(attr, "eppd") %>%
  lapply("[[", "ancsitedat") %>%
  Map(data.frame, unregion = unregion, country = country, eppregion = eppregion, .) %>%
  bind_rows()

ancrtcens %>% count(year) %>% as.data.frame
ancsite %>% count(year) %>% as.data.frame

#' # 3. Update FRR parameters 

frr_country <- here::here("inputs", "frr_country_sex12m.csv") %>%
  read.csv(stringsAsFactors = FALSE) %>%
  mutate(sex12m_z = (sex12m - 40) / 15)
frr_age <- here::here("inputs", "frr_age.csv") %>% read.csv(stringsAsFactors = FALSE)
frr_cd4cat <- here::here("inputs", "frr_cd4.csv") %>% read.csv(stringsAsFactors = FALSE)
frr_art_age <- here::here("inputs", "frr_art.csv") %>% read.csv(stringsAsFactors = FALSE)
frr_15to19sex12m <- here::here("inputs", "frr_15to19sex12m_z.csv") %>% read.csv(header = FALSE) %>% as.numeric


add_frr_noage_fp <- function(obj){

  country_name <- attr(obj, "country")
  
  region_name <- filter(frr_country, country == country_name)$region
  
  if(length(region_name) != 1)
    stop(paste("Country", country_name, "not found by create_frr_fp()"))

  frr_beta_country <- filter(frr_country, country == country_name)$frr_country
  if(length(frr_beta_country) != 1)
    frr_beta_country <- 1.0

  frr_beta_cd4 <- frr_cd4cat$est
  frr_beta_age <- filter(frr_age, region == region_name)$est * frr_beta_country

  sex12m_z <- filter(frr_country, country == country_name)$sex12m_z
  frr_15to19 <- frr_beta_age[1] * frr_15to19sex12m ^ sex12m_z

  frr_beta_art_age <- frr_art_age$est * frr_beta_country

  ## Construct FRR parameter inputs
  frr_cd4 <- array(NA, c(7, 8, 52))
  frr_cd4[] <- outer(frr_beta_cd4, c(frr_15to19, frr_15to19, frr_beta_age[2:7]), "*")
  
  frr_art <- array(NA, c(3, 7, 8, 52))
  frr_art[1,,,] <- frr_cd4
  frr_art[2:3,,,] <- rep(frr_beta_art_age[c(1, 1:7)], each = 2*7)

  attr(obj, "specfp")$frr_cd4 <- frr_cd4
  attr(obj, "specfp")$frr_art <- frr_art

  obj
}

inputs <- inputs %>% lapply(add_frr_noage_fp)


saveRDS(inputs, here::here("inputs", "eppasm-inputs.rds"))
saveRDS(hhsprev, here::here("inputs", "hhsprev.rds"))
saveRDS(ancsite, here::here("inputs", "ancsite.rds"))
saveRDS(ancrtcens, here::here("inputs", "ancrtcens.rds"))
