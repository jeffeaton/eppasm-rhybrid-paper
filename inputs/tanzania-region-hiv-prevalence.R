library(rdhs)
library(tidyverse)
library(survey)


#' ## Target regions
#'
#' Load shape file for current regions from 2017 MIS

sh17 <- st_read("~/Documents/Data/shape files/DHS/TZ2017MIS/shps/sdr_subnational_boundaries.shp")


#' ## Tanzania 2003 AIS

dhs_datasets(countryIds = "TZ",
             surveyYear = 2003,
             fileFormat = "FL",
             fileType = c("IR", "MR", "AR", "GE"))

ir <- get_datasets("TZIR4AFL.zip")[[1]] %>% readRDS
ar <- get_datasets("tzar4afl.zip")[[1]] %>% readRDS
ge <- get_datasets("TZGE4CFL.zip")[[1]] %>% readRDS

## Load shape files for survey admin areas
sh03 <- st_read("~/Documents/Data/shape files/DHS/TZ2003AIS/shps/sdr_subnational_boundaries.shp")

ggplot() +
  geom_sf(data = mutate(sh03, region2003 = REGNAME),
          aes(fill = region2003), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA)
  
## Assign coordiantes to survey cluster dataset
ge <- st_as_sf(ge, coords = c("LONGNUM", "LATNUM"))


## Merge shape file polygons to survey clusters
ge <- st_join(ge, sh, join = st_intersects)

## Confirm clusters that changed region
ggplot() +
  geom_sf(data = mutate(sh03, region2003 = REGNAME),
          aes(fill = region2003), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf(data = filter(ge, as.character(DHSREGNA) != as.character(REGNAME)),
          size = 0.5)

ggplot() +
  geom_sf(data = mutate(sh03, region2003 = REGNAME),
          aes(fill = region2003), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf_text(data = filter(ge, as.character(DHSREGNA) != as.character(REGNAME)),
               aes(label = DHSCLUST), size = 2.0)

## Three clusters look erroneously reclassified due to geomasking: 112, 40, 345

filter(ge, DHSCLUST %in% c(112, 40, 345))
ge$REGNAME[ge$DHSCLUST == 345] <- "Iringa"
ge$REGNAME[ge$DHSCLUST == 40] <- "Tanga"
ge$REGNAME[ge$DHSCLUST == 112] <- "Dodoma"

set_na <- function(x, na_codes = 9){ x[x %in% na_codes] <- NA; x }

dat <- ar %>%
  left_join(transmute(ge, hivclust = as.integer(DHSCLUST), region = REGNAME)) %>%  
  left_join(
    transmute(ir,
              hivclust = v001,
              hivnumb = v002,
              hivline = v003,
              age = v012,
              sex = as_factor(aidsex))
  ) %>%
  mutate(prev = as.integer(set_na(hiv03, 4:9) > 0),
         hivweight = hiv05 / 1e6,
         region = fct_collapse(region, "Zanzibar" = c("Kaskazini Pemba", "Kaskazini Unguja", "Kusini Pemba", "Kusini Unguja", "Mjini Magharibi")),
         region = factor(region))

filter(dat, is.na(region))

des <- svydesign(ids = ~hivclust, data = filter(dat, !is.na(hivweight)), weights = ~hivweight, nest = TRUE)

## Check that national age 15-49 prevalence matches Stat compiler
svyciprop(~prev, subset(des, age %in% 15:49), na.rm=TRUE)
dhs_data(countryIds = "TZ", indicatorIds = "HA_HIVP_B_HIV")

tz03 <- svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), svymean) %>%
  left_join(svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), unwtd.count) %>% select(-se)) %>%
  data.frame(surveyid = "TZ2003AIS", year = 2003, .)



#' ## Tanzania 2007 AIS

dhs_datasets(countryIds = "TZ",
             surveyYear = 2007,
             fileFormat = "FL",
             fileType = c("IR", "MR", "AR", "GE"))

ir <- get_datasets("TZIR51FL.zip")[[1]] %>% readRDS
ar <- get_datasets("tzar51fl.zip")[[1]] %>% readRDS
ge <- get_datasets("TZGE52FL.zip")[[1]] %>% readRDS

## Load shape files for survey admin areas
sh07 <- st_read("~/Documents/Data/shape files/DHS/TZ2007AIS/shps/sdr_subnational_boundaries.shp")

ggplot() +
  geom_sf(data = mutate(sh07, region2007 = REGNAME),
          aes(fill = region2007), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA)
  
## Assign coordiantes to survey cluster dataset
ge <- st_as_sf(ge, coords = c("LONGNUM", "LATNUM"))

## Merge shape file polygons to survey clusters
ge <- st_join(ge, sh, join = st_intersects)

## Clusters missing geocodes, assign to old region
filter(ge, is.na(REGNAME))

ge <- mutate(ge, REGNAME = ifelse(is.na(REGNAME), as.character(ADM1FIPSNA), as.character(REGNAME)))



## Confirm clusters that changed region
ggplot() +
  geom_sf(data = mutate(sh07, region2007 = REGNAME),
          aes(fill = region2007), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf(data = filter(ge, as.character(ADM1FIPSNA) != as.character(REGNAME)),
          size = 0.5)

ggplot() +
  geom_sf(data = mutate(sh07, region2007 = REGNAME),
          aes(fill = region2007), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf_label(data=sh17, aes(label = REGNAME), color = "red", size = 1.0) +
  geom_sf_text(data = filter(ge, as.character(ADM1FIPSNA) != as.character(REGNAME)),
               aes(label = DHSCLUST), size = 2.0)

## Three clusters look erroneously reclassified due to geomasking: 55, 114, 438

filter(ge, DHSCLUST %in% c(55, 114, 438))
ge$REGNAME[ge$DHSCLUST == 55] <- "Manyara"
ge$REGNAME[ge$DHSCLUST == 114] <- "Pwani"
ge$REGNAME[ge$DHSCLUST == 438] <- "Manyara"

set_na <- function(x, na_codes = 9){ x[x %in% na_codes] <- NA; x }

dat <- ar %>%
  left_join(transmute(ge, hivclust = as.integer(DHSCLUST), region = REGNAME)) %>%  
  left_join(
    transmute(ir,
              hivclust = v001,
              hivnumb = v002,
              hivline = v003,
              age = v012,
              sex = as_factor(aidsex))
  ) %>%
  mutate(prev = as.integer(set_na(hiv03, 4:9) > 0),
         hivweight = hiv05 / 1e6,
         region = fct_collapse(region, "Zanzibar" = c("Kaskazini Pemba", "Kaskazini Unguja", "Kusini Pemba", "Kusini Unguja", "Mjini Magharibi")),
         region = fct_relevel(region, "Zanzibar", after = Inf))


des <- svydesign(ids = ~hivclust, data = filter(dat, !is.na(hivweight)), weights = ~hivweight, nest = TRUE)

## Check that national age 15-49 prevalence matches Stat compiler
svyciprop(~prev, subset(des, age %in% 15:49), na.rm=TRUE)
dhs_data(countryIds = "TZ", indicatorIds = "HA_HIVP_B_HIV")

tz07 <- svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), svymean) %>%
  left_join(svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), unwtd.count) %>% select(-se)) %>%
  data.frame(surveyid = "TZ2007AIS", year = 2007, .)


#' ## Tanzania 2012 AIS

dhs_datasets(countryIds = "TZ",
             surveyYear = 2012,
             fileFormat = "FL",
             fileType = c("IR", "MR", "AR", "GE"))

ir <- get_datasets("TZIR6AFL.zip")[[1]] %>% readRDS
ar <- get_datasets("tzar6Afl.zip")[[1]] %>% readRDS
ge <- get_datasets("TZGE6AFL.zip")[[1]] %>% readRDS

## Load shape files for survey admin areas
sh12 <- st_read("~/Documents/Data/shape files/DHS/TZ2012AIS/shps/sdr_subnational_boundaries.shp")

ggplot() +
  geom_sf(data = mutate(sh12, region2012 = REGNAME),
          aes(fill = region2012), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA)
  
## Assign coordiantes to survey cluster dataset
ge <- st_as_sf(ge, coords = c("LONGNUM", "LATNUM"))

## Merge shape file polygons to survey clusters
ge <- st_join(ge, sh, join = st_intersects)

## Clusters missing geocodes, assign to old region
filter(ge, is.na(REGNAME))

ge <- mutate(ge, REGNAME = ifelse(is.na(REGNAME), as.character(ADM1NAME), as.character(REGNAME)))



## Confirm clusters that changed region
ggplot() +
  geom_sf(data = mutate(sh12, region2012 = REGNAME),
          aes(fill = region2012), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf(data = filter(ge, as.character(ADM1NAME) != as.character(REGNAME)),
          size = 0.5)

ggplot() +
  geom_sf(data = mutate(sh12, region2012 = REGNAME),
          aes(fill = region2012), color = "blue", alpha = 0.3) +
  geom_sf(data=sh17, color = "red", fill = NA) +
  geom_sf_label(data=sh17, aes(label = REGNAME), color = "red", size = 1.0) +
  geom_sf_text(data = filter(ge, as.character(ADM1NAME) != as.character(REGNAME)),
               aes(label = DHSCLUST), size = 2.0)

## All looks good, no clusters look erroneously reclassified due to geomasking

set_na <- function(x, na_codes = 9){ x[x %in% na_codes] <- NA; x }

dat <- ar %>%
  left_join(transmute(ge, hivclust = as.integer(DHSCLUST), region = REGNAME)) %>% 
  left_join(
    transmute(ir,
              hivclust = v001,
              hivnumb = v002,
              hivline = v003,
              age = v012,
              sex = as_factor(aidsex))
  ) %>%
  mutate(prev = as.integer(set_na(hiv03, 4:9) > 0),
         hivweight = hiv05 / 1e6,
         region = fct_collapse(region, "Zanzibar" = c("Kaskazini Pemba", "Kaskazini Unguja", "Kusini Pemba", "Kusini Unguja", "Mjini Magharibi")),
         region = fct_relevel(region, "Zanzibar", after = Inf))

des <- svydesign(ids = ~hivclust, data = filter(dat, !is.na(hivweight)), weights = ~hivweight, nest = TRUE)

## Check that national age 15-49 prevalence matches Stat compiler
svyciprop(~prev, subset(des, age %in% 15:49), na.rm=TRUE)
dhs_data(countryIds = "TZ", indicatorIds = "HA_HIVP_B_HIV")

tz12 <- svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), svymean) %>%
  left_join(svyby(~prev, ~region, subset(des, age %in% 15:49 & !is.na(prev)), unwtd.count) %>% select(-se)) %>%
  data.frame(surveyid = "TZ2012AIS", year = 2012, .)

#' ## Save outputs

rbind(tz03, tz07, tz12) %>%
  arrange(region, year) %>%
  write.csv("tanzania-region-hiv-prevalence.csv", row.names = FALSE)
