# Correlated evolution of ndh loss and nutritional strategy (Pagel 1994 test, phytools::fitPagel).
# usage: Rscript pagel.R <chronogram.nwk> <phaseC_traits.csv> <out_prefix> [threshold=3]
# For each heterotrophic strategy, the tree is pruned to that strategy + autotrophs (other heterotrophs removed;
# Byblis removed because its ndh state is unresolved), and an independent model is compared with a model in which
# ndh transition rates depend on the strategy (dep.var = "x"), and with a fully dependent model.
suppressMessages({library(ape); library(phytools)})
a <- commandArgs(TRUE); thr <- if (length(a) >= 4) as.integer(a[4]) else 3
ch <- read.tree(a[1]); tv <- read.csv(a[2], stringsAsFactors = FALSE)
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
tv$ndh <- ifelse(grepl("\\?", tv$ndh_calls), NA, ifelse(nbad >= thr, "lost", "intact"))
groups <- list(carnivory = "carnivorous", parasitism = c("hemiparasitic", "holoparasitic"),
               mycoheterotrophy = c("full_mycoheterotroph", "partial_mycoheterotroph"),
               any_heterotrophy = c("carnivorous", "hemiparasitic", "holoparasitic", "full_mycoheterotroph", "partial_mycoheterotroph"))
res <- list()
for (g in names(groups)) {
  keep <- tv$tip[!is.na(tv$ndh) & (tv$nutrition == "autotrophic" | tv$nutrition %in% groups[[g]])]
  t2 <- keep.tip(ch, keep); d <- tv[match(t2$tip.label, tv$tip), ]
  x <- setNames(d$ndh, d$tip); y <- setNames(ifelse(d$nutrition == "autotrophic", "autotrophic", "heterotrophic"), d$tip)
  fx <- fitPagel(t2, x, y, dep.var = "x")
  fxy <- fitPagel(t2, x, y, dep.var = "xy")
  Q <- fx$dependent.Q
  loss_auto <- Q["intact|autotrophic", "lost|autotrophic"]; loss_het <- Q["intact|heterotrophic", "lost|heterotrophic"]
  res[[g]] <- data.frame(contrast = paste(g, "vs autotrophy"), n_tips = Ntip(t2),
    n_lost_auto = sum(x == "lost" & y == "autotrophic"), n_lost_het = sum(x == "lost" & y == "heterotrophic"),
    n_het = sum(y == "heterotrophic"),
    logL_indep = fx$independent.logL, logL_dep_ndh = fx$dependent.logL, LR_dep_ndh = 2 * (fx$dependent.logL - fx$independent.logL),
    p_dep_ndh = fx$P, AIC_indep = fx$independent.AIC, AIC_dep_ndh = fx$dependent.AIC, AIC_dep_both = fxy$dependent.AIC,
    p_dep_both = fxy$P,
    ndh_loss_rate_autotrophic = loss_auto, ndh_loss_rate_heterotrophic = loss_het, ratio = loss_het / loss_auto)
  print(res[[g]])
}
R <- do.call(rbind, res); write.csv(R, paste0(a[3], "_pagel.csv"), row.names = FALSE)
cat("done\n")
