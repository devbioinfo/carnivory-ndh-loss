suppressMessages({library(ape); library(phytools)})
ch <- read.tree("f311_GBM_t3_chronogram.nwk"); tv <- read.csv("phaseC_traits_cov311.csv", stringsAsFactors = FALSE)
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
tv$ndh <- ifelse(grepl("\\?", tv$ndh_calls), NA, ifelse(nbad >= 3, "lost", "intact"))
keep <- tv$tip[!is.na(tv$ndh) & tv$nutrition %in% c("autotrophic", "carnivorous")]
t2 <- keep.tip(ch, keep); d <- tv[match(t2$tip.label, tv$tip), ]
x <- setNames(d$ndh, d$tip); y <- setNames(ifelse(d$nutrition == "carnivorous", "carnivorous", "autotrophic"), d$tip)
for (dv in c("x", "xy")) { f <- fitPagel(t2, x, y, dep.var = dv)
  Q <- f$dependent.Q; write.csv(Q, paste0("pagelQ_carn_dep", dv, ".csv"))
  write.csv(f$independent.Q, paste0("pagelQ_carn_indep.csv"))
  cat(dv, f$independent.logL, f$dependent.logL, f$P, f$independent.AIC, f$dependent.AIC, "\n") }
