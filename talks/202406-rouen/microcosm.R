library(PLNmodels)
library(tidyverse)
library(viridis)
library(lemon)

microcosm <- readRDS("microcosm_reduced.rds")
microcosm$Abundance <- microcosm$Abundance[, colMeans(microcosm$Abundance > 0) > 0.05]
microcosm$site_time <- droplevels(microcosm$site_time)

PLN           <- PLN(Abundance ~ 1        + offset(log(Offset)), data = microcosm)
PLN_site      <- PLN(Abundance ~ 0 + site + offset(log(Offset)), data = microcosm)
PLN_time      <- PLN(Abundance ~ 0 + time + offset(log(Offset)), data = microcosm)
PLN_site_time <- PLN(Abundance ~ 0 + site_time + offset(log(Offset)), data = microcosm)
ZI            <- ZIPLN(Abundance ~ 1 + offset(log(Offset))           , data = microcosm)
ZI_site       <- ZIPLN(Abundance ~ 1 + offset(log(Offset)) | 0 + site,  data = microcosm)
ZI_time       <- ZIPLN(Abundance ~ 1 + offset(log(Offset)) | 0 + time,  data = microcosm)
ZI_site_time  <- ZIPLN(Abundance ~ 1 + offset(log(Offset)) | 0 + site_time,  data = microcosm)
ZI_site_PLN_site <- ZIPLN(Abundance ~ 0 + site + offset(log(Offset)) | 0 + site,  data = microcosm)
ZI_time_PLN_time <- ZIPLN(Abundance ~ 0 + time + offset(log(Offset)) | 0 + time,  data = microcosm)
ZIPLN_time_PLN_site <- ZIPLN(Abundance ~ 0 + site_time + offset(log(Offset)) | 0 site_time,  data = microcosm)

ZI_site_PLN_time <- ZIPLN(Abundance ~ 0 + time + offset(log(Offset)) | 0 + site,  data = microcosm)
ZI_site_PLN_site_time <- ZIPLN(Abundance ~ 0 + site + time + offset(log(Offset)) | 0 + site + time,  data = microcosm)

rbind(
  PLN$criteria,
  PLN_site$criteria,
  PLN_time$criteria,
  PLN_site_time$criteria,
  ZI$criteria,
  ZI_site$criteria,
  ZI_time$criteria,
  ZI_site_PLN_site$criteria,
  ZI_time_PLN_time$criteria
) %>%
  as.data.frame(row.names = c("PLN", "PLN site", "PLN time", "PLN site and time",
                              "ZI single", "ZI site", "ZI time",
                              "ZI site PLN site", "ZI time PLN time")) %>%
  knitr::kable()

# distribution of zeros
true_zeros <- microcosm$Abundance == 0
n_zeros <- sum(true_zeros)
n_nzeros <- sum(!true_zeros)
distr_zeros <-
  data.frame(
  fitted   = c(as.vector(fitted(PLN)[true_zeros]), as.vector(fitted(ZI)[true_zeros]),
               as.vector(fitted(ZI_time)[true_zeros]), as.vector(fitted(ZI_site)[true_zeros])),
  obs      = rep(1:n_zeros, 4),
  method   = factor(rep(c("PLN", "ZIPLN (single param)", "ZI (covar time)", "ZI (covar site)"), each = n_zeros))
) %>%
  ggplot() +
  geom_rug(aes(x = fitted, color = method)) +
  geom_point(aes(x = obs, y = fitted, color = method), size = 0.3, alpha = 0.1) +
  ggtitle("Distribution of observed zeros as fitted by the models") +
  xlab(label = "observations #") +
  scale_color_viridis_d() + theme_bw()

distr_zeros_logscale <-
  data.frame(
    fitted   = c(as.vector(fitted(PLN_std)[true_zeros]), as.vector(fitted(ZI_single)[true_zeros]),
                 as.vector(fitted(ZI_cols)[true_zeros]), as.vector(fitted(ZI_site)[true_zeros])),
    obs      = rep(1:n_zeros, 4),
    method   = factor(rep(c("PLN", "ZIPLN (single param)", "ZI (1 par per species)", "ZI (1 per group of site)"), each = n_zeros))
  )
) %>%
  ggplot() +
  geom_rug(aes(x = fitted, color = method), show.legend = FALSE) +
  geom_point(aes(x = obs, y = fitted, color = method), size = 0.3, alpha = 0.1, show.legend = FALSE) +
  scale_y_log10() +
  ggtitle("Distribution of observed zeros as fitted by the models (log scale)") +
  xlab("observations #") +
  scale_color_viridis_d() + theme_bw()

# distribution of non zero
data.frame(
  fitted   = c(as.vector(fitted(PLN)[!true_zeros]), as.vector(fitted(ZI)[!true_zeros]),
               as.vector(fitted(ZI_time)[!true_zeros]), as.vector(fitted(ZI_site)[!true_zeros])),
  obs      = rep(1:n_nzeros, 4),
  method   = factor(rep(c("PLN", "ZIPLN (single param)", "ZI (1 par per species)", "ZI (1 par group of site)"), each = n_nzeros))
) %>%
  ggplot() + aes(x = fitted, fill = method, group = method) +
  ggtitle("Distribution of observed non zeros as fitted by the models") +
  geom_rug(show.legend = FALSE)  +
  geom_histogram() + scale_x_log10() + scale_fill_viridis_d() + theme_bw()

# distribution of non zero
data.frame(
  fitted   = c(as.vector(fitted(PLN)[!true_zeros]), as.vector(fitted(ZI)[!true_zeros]),
               as.vector(fitted(ZI_time)[!true_zeros]), as.vector(fitted(ZI_site)[!true_zeros])),
  observed = rep(c(microcosm$Abundance[!true_zeros]), 4),
  method   = factor(rep(c("PLN", "ZIPLN (single param)", "ZI (1 par per species)", "ZI (1 par group of site)"), each = n_nzeros))
) %>%
  ggplot(aes(x = observed, y = fitted)) +
  geom_point(size = .5, alpha =.25 ) +
  facet_wrap( ~ method) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + annotation_logticks()

heatmap(coef(ZI_site_PLN_site, "count"))
heatmap(coef(ZI_site_PLN_site, "zero"))
heatmap(coef(ZI_site, "zero"))

prcomp(ZI_site$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()
prcomp(ZI_site$latent_pos) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()
prcomp(ZI_site_PLN_site$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()
prcomp(ZI_site_PLN_site$latent_pos) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()

prcomp(ZI_PLN_time$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()
prcomp(ZI_time_site$latent) %>% factoextra::fviz_pca_ind(axes = c(1,2), col.ind = microcosm$site) + scale_color_viridis_d()


site_time_palette <- c(
  "_**<span style = 'color:#6baed6;'>Buccal</span>**_"  = "transparent",
  "B -1W" = "#bdd7e7", 
  "B 1M" = "#6baed6", 
  "B 3M" = "#2171b5", 
  "B 7M" = "#053061",
  "_**<span style = 'color:#fb6a4a;'>Nasal</span>**_"  = "transparent",
  "N -1W" = "#fcae91", 
  "N 1M" = "#fb6a4a", 
  "N 3M" = "#cb181d",
  "N 7M" = "#67001F",
  "_**<span style = 'color:#bf812d;'>Vaginal</span>**_"  = "transparent",
  "V -1W" = "#dfc27d", 
  "V 1M" = "#bf812d", 
  "V 3M" = "#8c510a",
  "V 7M" = "#543005",
  "_**<span style = 'color:#74c476;'>Milk</span>**_"  = "transparent",
  "L -1W" = "#bae4b3", 
  "L 1M" = "#74c476", 
  "L 3M" = "#238b45", 
  "L 7M" = "#00441B"
)
site_time_palette_simple <- site_time_palette[-c(1,6,11,16)]
index_palette <- c("minus" = "#67a9cf", 
                   "plus" = "#ef8a62")
site_time_design <- strip_nested(
  background_x = elem_list_rect(fill = c(
    ## Time
    rep("gray", 4), 
    ## Site within time
    site_time_palette_simple[c(1, 5, 9, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16)]))
)