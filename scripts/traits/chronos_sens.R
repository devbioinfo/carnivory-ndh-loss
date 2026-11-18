suppressMessages({library(ape); library(phytools)})
tr <- read.tree(commandArgs(TRUE)[1]); tr$node.label <- NULL
tr <- root(tr, "Amborella_trichopoda", resolve.root = TRUE); tr$edge.length[tr$edge.length <= 0] <- 1e-6
cal <- read.csv("calib_ML_GBM.csv")
nd <- mapply(function(x, y) getMRCA(tr, c(x, y)), cal$tipA, cal$tipB)
calib <- data.frame(node = nd, age.min = cal$age_mean, age.max = cal$age_mean, soft.bounds = FALSE)
calib <- calib[!duplicated(calib$node), ]
main <- read.tree("f311_GBM_t3_chronogram.nwk")
ages <- function(t) { H <- max(node.depth.edgelength(t)); a <- H - node.depth.edgelength(t); a[(Ntip(t)+1):(Ntip(t)+t$Nnode)] }
key <- function(t) sapply((Ntip(t)+1):(Ntip(t)+t$Nnode), function(n) paste(sort(extract.clade(t, n)$tip.label), collapse="|"))
km <- key(main); am <- ages(main); names(am) <- km
free <- setdiff(km, key(tr)[calib$node - Ntip(tr)])
cat(sprintf("internal nodes %d; fixed %d; smoothed %d; root %.2f Ma\n", main$Nnode, nrow(calib), length(free), max(am)))
out <- list()
for (cfg in list(c("correlated",0.1), c("correlated",10), c("relaxed",1), c("discrete",1))) {
  ch <- chronos(tr, lambda = as.numeric(cfg[2]), model = cfg[1], calibration = calib, quiet = TRUE,
                control = chronos.control(iter.max = 1e4, eval.max = 1e4, nb.rate.cat = 10))
  class(ch) <- "phylo"; a <- ages(ch); names(a) <- key(ch)
  d <- a[free] - am[free]; rel <- abs(d) / am[free]
  cat(sprintf("%s lambda=%s: PHIIC=%.1f  median |dAge| %.2f Myr (%.1f%%), 95th pct %.2f Myr, max %.2f Myr; r=%.4f\n",
      cfg[1], cfg[2], attr(ch, "PHIIC")$PHIIC, median(abs(d)), 100*median(rel), quantile(abs(d), .95), max(abs(d)), cor(a[free], am[free])))
}
