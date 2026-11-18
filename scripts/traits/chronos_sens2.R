suppressMessages({library(ape); library(phytools)})
tr <- read.tree(commandArgs(TRUE)[1]); tr$node.label <- NULL
tr <- root(tr, "Amborella_trichopoda", resolve.root = TRUE); tr$edge.length[tr$edge.length <= 0] <- 1e-6
cal <- read.csv("calib_ML_GBM.csv")
nd <- mapply(function(x, y) getMRCA(tr, c(x, y)), cal$tipA, cal$tipB)
calib <- data.frame(node = nd, age.min = cal$age_mean, age.max = cal$age_mean, soft.bounds = FALSE)
calib <- calib[!duplicated(calib$node), ]
tv <- read.csv("phaseC_traits_cov311.csv", stringsAsFactors = FALSE)
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
tv$ndh <- ifelse(grepl("\\?", tv$ndh_calls), NA, ifelse(nbad >= 3, "lost", "intact"))
pag <- function(ch, g) {
  keep <- tv$tip[!is.na(tv$ndh) & (tv$nutrition == "autotrophic" | tv$nutrition %in% g)]
  t2 <- keep.tip(ch, keep); d <- tv[match(t2$tip.label, tv$tip), ]
  x <- setNames(d$ndh, d$tip); y <- setNames(ifelse(d$nutrition == "autotrophic", "autotrophic", "heterotrophic"), d$tip)
  fx <- fitPagel(t2, x, y, dep.var = "x"); Q <- fx$dependent.Q
  c(LR = 2 * (fx$dependent.logL - fx$independent.logL), P = fx$P,
    ratio = Q["intact|heterotrophic", "lost|heterotrophic"] / Q["intact|autotrophic", "lost|autotrophic"])
}
trees <- list(main_correlated_l1 = read.tree("f311_GBM_t3_chronogram.nwk"))
for (cfg in list(c("relaxed", 1), c("discrete", 1))) {
  ch <- chronos(tr, lambda = as.numeric(cfg[2]), model = cfg[1], calibration = calib, quiet = TRUE,
                control = chronos.control(iter.max = 1e4, eval.max = 1e4, nb.rate.cat = 10))
  cat(cfg[1], "PHIIC", attr(ch, "PHIIC")$PHIIC, "conv msg:", attr(ch, "message"), "\n")
  class(ch) <- "phylo"; ch$edge.length[ch$edge.length <= 0] <- 1e-4
  write.tree(ch, paste0("sens_chronos_", cfg[1], ".nwk")); trees[[paste0(cfg[1], "_l1")]] <- ch
}
for (n in names(trees)) {
  a <- pag(trees[[n]], "carnivorous"); b <- pag(trees[[n]], c("hemiparasitic", "holoparasitic"))
  cat(sprintf("%s  carnivory LR %.1f P %.2e ratio %.0f | parasitism LR %.1f P %.2e ratio %.0f\n", n, a[1], a[2], a[3], b[1], b[2], b[3]))
}
