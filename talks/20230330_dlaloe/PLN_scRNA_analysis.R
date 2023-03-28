## ----glimpse Abundance--------------------------------------------------------
data(oaks)
oaks$Abundance %>% as_tibble() %>% 
  dplyr::select(1:10) %>% 
  rmarkdown::paged_table()


## ----glance Abundances, fig.height=6------------------------------------------
log(1 + oaks$Abundance) %>% 
  corrplot::corrplot(is.corr = FALSE,
    addgrid.col = NA,  tl.cex = .5,  cl.pos = "n")


## ----simple PLN offsets, cache = TRUE, results = FALSE------------------------
M01_oaks <- PLN(Abundance ~ 1 + offset(log(Offset)) , oaks)


## ----PLN covariate oaks, cache = TRUE, results = FALSE------------------------
M11_oaks <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), oaks)


## ----PLN regroup oaks modalities, cache = TRUE, results = FALSE---------------
M21_oaks <- PLN(Abundance ~  0 + tree + orientation + offset(log(Offset)), oaks)


## ----PLN covariate oaks results-----------------------------------------------
rbind(M01 = M01_oaks$criteria,
      M11 = M11_oaks$criteria, M21 = M21_oaks$criteria) %>% 
  knitr::kable(format = "html")


## ----oaks matrix plot, fig.width=14, fig.height=2, echo = FALSE---------------
coef(M11_oaks) %>% t() %>% corrplot(method = "color", is.corr = FALSE, tl.cex = 1, cl.pos = "n")


## ----PLNLDA oaks, results = FALSE---------------------------------------------
myLDA_tree <- 
  PLNLDA(Abundance ~ 1 + offset(log(Offset)), grouping = oaks$tree, data = oaks)


## ----plot oaks1, echo = FALSE, fig.height = 6---------------------------------
plot(myLDA_tree, map = "individual")


## ----plot oaks2, echo = FALSE, fig.height = 6---------------------------------
plot(myLDA_tree, map = "variable")


## ----PLNPCA offset, cache = TRUE, results = FALSE-----------------------------
PCA_offset <- 
  PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = oaks, ranks = 1:30)
PCA_offset_BIC <- getBestModel(PCA_offset, "BIC")


## ----PCA offset vizu tree, fig.width=6, fig.height=6, fig.align="center", echo=FALSE----
factoextra::fviz_pca_biplot(
  PCA_offset_BIC, select.var = list(contrib = 10), addEllipses = TRUE, habillage = oaks$tree,
  title = "Biplot (10 most contributing species)"
  ) + labs(col = "tree status") + scale_color_viridis_d()


## ----PCA covariate tree, cache = TRUE, results = FALSE, warning=FALSE, message=FALSE----
PCA_tree <- 
  PLNPCA(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks, ranks = 1:30)


## ----PCA covariate tree plot, echo = FALSE, fig.align="center", fig.width=6, fig.height=6----
PCA_tree %>% getBestModel("BIC") %>% 
factoextra::fviz_pca_biplot(
  select.var = list(contrib = 10), col.ind = oaks$distTOground,
  title = "Biplot (10 most contributing species)"
  ) + labs(col = "distance (cm)") + scale_color_viridis_c()


## ----PLNmixture oaks, cache = TRUE--------------------------------------------
PLN_mixtures <-
   PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = oaks, clusters = 1:3, control_main = list(iterates=0))
myPLN_mix <- getModel(PLN_mixtures, 3)


## ----PLN clustering 1 fake, eval = FALSE--------------------------------------
## myPLN_mix$plot_clustering_pca()


## ----PLN clustering 1, echo = FALSE-------------------------------------------
myPLN_mix$plot_clustering_pca(main = 'clustering memberships in individual factor map')


## ----PLN clustering 2---------------------------------------------------------
myPLN_mix$plot_clustering_data()


## ----PLNnetwork, cache = TRUE, results='hide', message = FALSE----------------
networks <- PLNnetwork(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks)


## ----PLNnetwork criteria, echo = FALSE, fig.align='center', fig.height=6------
plot(networks)


## ----PLNnetwork network, echo = FALSE, fig.align='center', fig.height=6-------
type <- rep(c("bacterial", "fungal", "Alphiltoides"), c(66,114-67,1))
net <- fancy_network(getBestModel(networks, "EBIC"), type)
print(net)

