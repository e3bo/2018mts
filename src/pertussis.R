#!/usr/bin/env Rscript

library(lubridate)
library(padr) ## for pad()
library(tidyverse)
library(tsibble) ## for slide*_dbl()
library(ggforce) ## for facet_wrap_paginate()
library(gganimate)

## National level, yearly aggregated analysis

pertu <- read.table("usa-national.csv")
pertu[,2] <- as.integer(gsub(",", "", pertu[,2]))
names(pertu) <- c("year", "USA")

pertc <- read.table("candada.csv", header = TRUE)
names(pertc) <- c("year", "Canada")

pertu %>% inner_join(pertc, by="year") -> pertna
pertna %>% gather("Canada", "USA", key = "country", value = "cases") -> perttd

find_peaks <- function(x){
  st <- diff(x) > 0
  n <- length(st)
  peaks <- integer()
  for(i in seq(2, n)){
    if(!st[i] && st[i - 1]){
      peaks <- c(peaks, i)
    }
  }
  peaks
}

ca_peaks <- find_peaks(pertna$Canada) + min(pertna$year) - 1
us_peaks <- find_peaks(pertna$USA) + min(pertna$year) - 1

pdf <- data.frame(country = "Canada", pyear = ca_peaks)
pdf <- rbind(pdf, data.frame(country = "USA", pyear = us_peaks))

pts <- data.frame(year = pertna$year,
                  is_usa_pk = pertna$year %in% us_peaks,
                  is_can_pk = pertna$year %in% ca_peaks)

colwidth <- 84
wind <- 10
yrcut <- 1976
yrcutst <- yrcut - wind + 1
pts %>%
  mutate(USA_MA = slide_dbl(is_usa_pk, mean, .size =  wind)) %>%
  mutate(CA_MA = slide_dbl(is_can_pk, mean, .size = wind)) %>%
    mutate(UCcor = slide2_dbl(is_usa_pk, is_can_pk, cor, .size = wind)) %>%
      mutate(USA_MIEP = 1 / USA_MA, CA_MIEP = 1 / CA_MA) %>%
        mutate(year_first = year - wind + 1) -> pts_stats




xlabstrng <- paste0("First year of ", wind, "-year window")
pts_stats %>% select(CA_MIEP, USA_MIEP, UCcor, year_first) %>%
  gather("CA_MIEP", "USA_MIEP", "UCcor", key = "stat", value = "value") %>%
    mutate(stat_labels = factor(stat,
               levels = c("CA_MIEP", "USA_MIEP", "UCcor"),
               labels = c("Canada, years / peak", "USA, years / peak",
                   "Canada-USA, correlation of peak years"))) %>%
      filter(year_first < yrcutst) %>%
      ggplot(aes(x = year_first, y = value)) + geom_col(width = 0.5) +
        facet_wrap(~ stat_labels, ncol = 1, scales = "free_y") +
          ylab("Sliding window estimate") + xlab(xlabstrng) +
            theme_classic() -> g
ggsave("sliding-window-ests.pdf", plot = g, width = colwidth, units = "mm")


pts %>% gather("is_usa_pk", "is_can_pk", key = variable, value = "is_peak") %>%
  mutate(country = as.character(factor(variable,
             levels = c("is_usa_pk", "is_can_pk"),
             labels = c("USA", "Canada")))) %>%
    select(-"variable") -> ptsl

perttd %>% inner_join(ptsl, by = c("year", "country")) -> pdata

pdata %>% filter(year < yrcut) %>%
  ggplot(aes(x = year, y = cases, fill = is_peak)) +
  geom_col() +
  facet_wrap(~country, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Reported pertussis cases") +
  scale_fill_manual(values=c("#636363", "#bdbdbd"),
                       name="Type of year",
                       breaks=c(TRUE, FALSE),
                       labels=c("Peak", "Non-peak")) +
    theme_bw() -> g



ggsave("national-cases-ts.pdf", plot = g, width = colwidth * 2,
       height = colwidth, units = "mm")


## load in weekly, state-level data and perform consistency checks

pert <- read.csv("US.27836007.csv")

strt_cum <- table(pert$PeriodStartDate, pert$PartOfCumulativeCountSeries)
is_b4_76 <- ymd(rownames(strt_cum)) < ymd("1976-01-01")
stopifnot(all(strt_cum[is_b4_76,"1"] == 0))
## Assert that no cumulative time series in early data

pert %>% filter(is.na(Admin2Name) & is.na(CityName) &
                Admin1Name %in% c(toupper(state.name), "DISTRICT OF COLUMBIA") &
                  PartOfCumulativeCountSeries == 0) %>%
  mutate(start = ymd(PeriodStartDate),
         end = ymd(PeriodEndDate)) -> pbar

stopifnot(all(pbar$end - pbar$start == 6))
## Assert that all intervals are weekly

pbar %>% filter(start < ymd("1960-01-01")) %>%
  filter(start > ymd("1938-06-13")) %>%
  group_by(Admin1Name) %>%
    select(Admin1Name, CountValue, start) %>% pad() -> foobar



### Check the distribution of missingness cf. Magpantay et al.

foobar %>% group_by(Admin1Name) %>%
  summarise(fracNA = mean(is.na(CountValue))) -> fracNas

missingness_threshold <- 0.5
included_states <- fracNas$Admin1Name[fracNas$fracNA < missingness_threshold]
foobar %>% filter(Admin1Name %in% included_states) -> foo2


# Analysis

### Perform imputation via linear interpolation.

foo2 %>% group_by(Admin1Name) %>%
  mutate(ctsi = zoo::na.approx(CountValue), na.rm = FALSE) -> foo3
foo3 %>% group_by(Admin1Name) %>%
  summarise(fracNA = mean(is.na(ctsi))) %>% as.data.frame


nice_names <- c(state.name, "District of Columbia")
key <- match(tolower(foo3$Admin1Name),
             tolower(nice_names))

foo3$State <- nice_names[key]

foo3 %>% ggplot(aes(x = start, y = ctsi, fill = is.na(CountValue))) +
  geom_col() + facet_wrap(~State, scales = "free_y") +
  labs(x = "Start of week", y = "Reported pertussis cases") +
  scale_fill_manual(values=rev(c("#636363", "#bdbdbd")),
                       name="Type of observation",
                       breaks=c(TRUE, FALSE),
                    labels=c("Imputed", "Not imputed")) +
    theme_bw() +
   theme(legend.position = "top") -> g


g + facet_wrap_paginate(~State, nrow = 8, ncol = 3, page = 1,
                        scales = "free_y") -> g1
g + facet_wrap_paginate(~State, nrow = 8, ncol = 3, page = 2,
                        scales = "free_y") -> g2

ggsave("state-time-series-p1.pdf", plot = g1, width = 7.5, height = 10)
ggsave("state-time-series-p2.pdf", plot = g2, width = 7.5, height = 10)

### The interpolation seems to fill in values in a reasonable way.

## Now let's calculate synchrony in moving windows


foo3 %>% select(Admin1Name, start, ctsi) %>%
  spread(key = Admin1Name, value = ctsi) -> foo4

sync_dist <- function(tsmat, nrand = 100){
  ## implements bootstrapping procedure for evaluating confidence
  ## interval of average pairwise correllation in a matrix
  ##
  ## based on Statistical Methods in Ottar N. Bjørnstad, Nils Chr. Stenseth,
  ## Takashi Saitoh, Ecology, 1 March 1999
  n <- ncol(tsmat)
  syncs <- numeric(nrand)
  for(i in seq(0, nrand)){
    if (i != 0) {
      sel <- sample.int(n, replace = TRUE)
    } else {
      sel <- 1:n
    }
    cormat <- cor(tsmat[, sel], use = "pairwise.complete.obs")
    nterms <- cum <- 0
    for (j in 1:n){
      for (k in 1:n){
        if(sel[j] != sel[k] && !is.na(cormat[j,k])){
          nterms <- nterms + 1
          cum <- cum + cormat[j,k]
        }
      }
    }
    if (i == 0) {
      sync0 = cum / nterms
    } else {
      syncs[i] <- cum / nterms
    }
  }
  perc_ci <- quantile(syncs, probs = c(0.025,  0.975))
  list(syncs_bootstrapped = syncs, sync_data = sync0, perc_95ci = perc_ci)
}

get_syncs <- function(tsmat, winsize, nrand, remove_linear = TRUE,
                      pwr_transform = 0.5){
  sync <- list()
  winstarts <- seq_len(nrow(tsmat) - winsize + 1)
  for (i in winstarts){
    inds <- seq(i, i + winsize - 1)
    sel <- tsmat[inds, -1] ^ pwr_transform
    resids <- sel
    for(j in 1:ncol(sel)){
      y <- sel[, j]
      if(!any(is.na(y)) && remove_linear) {
        resids[, j] <- y - predict(lm(y~inds))
      } else {
        resids[, j] <- y
      }
    }
    sync[[i]] <- sync_dist(resids, nrand = nrand)
  }
  sync
}


syncs4 <- get_syncs(foo4, 52 * 8, nrand = 200)
point_ests <- sapply(syncs4, "[[", "sync_data")
ci <- sapply(syncs4, "[[", "perc_95ci")

start_inds <- seq_along(point_ests)
window_starts <- foo4$start[start_inds]
plot(window_starts, point_ests, ylim = c(0, .3))
lines(window_starts, ci[1,])
lines(window_starts, ci[2,])

## There's a slight increase, but it seems significant because the 95%
## percentile intervals are just overlapping.

system.time(syncs5 <- get_syncs(foo4, 52 * 8, pwr_transform = 1, nrand = 200))
point_ests5 <- sapply(syncs5, "[[", "sync_data")
ci5 <- sapply(syncs5, "[[", "perc_95ci")

inpermm <- 0.0394
pdf("synchrony-estimates.pdf", width = colwidth * inpermm,
    height = colwidth * inpermm * 1.2)
par(c(5, 4, 0, 0) + 0.1)
start_inds5 <- seq_along(point_ests5)
window_starts5 <- foo4$start[start_inds5]
plot(window_starts5, point_ests5, ylim = c(0, .3),
     ylab = "Average pairwise correlation", xlab = "First year of 8-year window", type = 'l')
lines(window_starts5, ci5[1,], col = "grey")
lines(window_starts5, ci5[2,], col = "grey")
dev.off()

## A similar result holds when the data are not square-root transformed.


foo3 %>% mutate(year = year(start)) %>%  filter(year > 1938 & year < 1955) %>%
  group_by(year) %>% summarise(annualct = sum(ctsi)) -> natagg

perttd %>% filter(country == "USA") %>% inner_join(natagg, by = "year") -> usdata
usdata %>% ggplot(aes(cases, annualct)) + geom_text(aes(label = year))


## As an additional consistency check we compare national totals of
## cases. Only 1939 looks different. This seems unlikely to affect the
## results of our statistical tests.

## Now check that dominant period also lengthens.

foo3 %>% select(Admin1Name, start, ctsi) %>% group_by(start) %>%
  summarise(usct = sum(ctsi)) -> usweekly

usweekly$lm <- predict(lm(usct ~ start, data = usweekly))
usweekly %>% mutate(resid = usct - lm,
                    start_decimal = decimal_date(start)) -> usweekly2

plot(resid ~ start, data = usweekly2)
bwd <- as.data.frame(usweekly2[, c("start_decimal", "resid")])
twt <- biwavelet::wt(bwd)

pdf("wavelet-spectrum.pdf", width = colwidth * inpermm)
par(mfrow = c(2, 1))
par(mar = c(1, 4.1, 4.1, 2.1))
ylim <- c(0, max(usweekly2$usct) * 1.05)
plot(usct ~ start, data = usweekly2, xaxs = "i", yaxs = "i", xaxt = "n",
     xlab = "", type = "h", ylim = ylim, col = "grey",
     ylab = "Reported pertussis cases")
lines(lm ~ start, data = usweekly2)
par(mar = c(5.1, 4.1, 0, 2.1))
plot(twt)
dev.off()

## Spatial analysis

## Animation of cases on map

library(sf)
library(tigris)

us_geo <- states(class = "sf")
laea <- st_crs("+proj=laea +lat_0=30 +lon_0=-95")
us_geo <- st_transform(us_geo, laea)

plot(us_geo["NAME"], graticule = TRUE, key.pos = NULL, axes = TRUE)

foo3 %>% inner_join(us_geo, by = c("State" = "NAME") ) %>% st_sf -> bar

pmap <- ggplot() + geom_sf(data = bar, aes(fill = ctsi, geometry = geometry)) +
  transition_manual(start) + labs(title = 'Reporting week: {current_frame}') +
    guides(fill = guide_legend(title = "Reported\npertussis\ncases"))

animate(pmap, renderer = ffmpeg_renderer())

anim_save("pertussis-animation")

iqrng <- function(x) unname(diff(quantile(x, probs = c(0.25, 0.75))))

foo3 %>% group_by(Admin1Name) %>%
  mutate(reports_cnts = (ctsi - quantile(ctsi, probs = 0.5)) / diff(range(ctsi))) -> foo7

foo7 %>% inner_join(us_geo, by = c("State" = "NAME") ) %>% st_sf %>%
ggplot(aes(fill = reports_cnts)) + geom_sf() +
  transition_manual(start) + labs(title = 'Reporting week: {current_frame}') +
    guides(fill = guide_legend(title = "Reported\npertussis\ncases\nscaled")) -> pmap2

animate(pmap2, renderer = ffmpeg_renderer())

anim_save("pertussis-animation-scaled.mp4")

## Relationship between strength of correlation with distance and population size

i1 <- 1:(nrow(foo4)/2)
i2 <- (nrow(foo4)/2 + 1):nrow(foo4)

key <- match(colnames(foo4)[-1], toupper(state.name))

longlat <- cbind(state.center$x[key], state.center$y[key])

n <- nrow(latlong)
dmat <- matrix(nrow = n, ncol = n)
for(i in seq_len(n)){
  for(j in seq_len(n)){
    dmat[i, j] <- geosphere::distHaversine(longlat[i,], longlat[j,])
  }
}

par(mfrow = c(2, 1))

cmat1 <- cor(foo4[i1,-1], use = "pairwise", method = "spearman")

p1 <- cbind(as.numeric(dmat[lower.tri(dmat)]), as.numeric(cmat1[lower.tri(cmat1)]))
plot(p1)
lines(lowess(p1, delta = diff(range(p1[,1], na.rm = TRUE)) / 100), col = 2)

cor(as.numeric(cmat1[lower.tri(cmat1)]), as.numeric(dmat[lower.tri(dmat)]), use="complete.obs", method = "kendall")

cmat2 <- cor(foo4[i2,-1], use = "pairwise", method = "spearman")

p2 <- cbind(as.numeric(dmat[lower.tri(dmat)]), as.numeric(cmat2[lower.tri(cmat2)]))
plot(p2)
lines(lowess(p2, delta = diff(range(p2[,1], na.rm = TRUE)) / 100), col = 2)

cor(as.numeric(cmat2[lower.tri(cmat2)]), as.numeric(dmat[lower.tri(dmat)]), use="complete.obs", method = "kendall")


plot(as.numeric(dmat[lower.tri(dmat)]), as.numeric((cmat2 - cmat1)[lower.tri(cmat2)]))

cor(foo4[i1, -1])
