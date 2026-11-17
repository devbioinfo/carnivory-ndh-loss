# Is the carnivory-ndh association explained by habitat (wet / aquatic) or CAM photosynthesis?
# usage: Rscript confounder.R <chronogram.nwk> <phaseC_traits_cov.csv> <out_prefix>
# 1. phylogenetic logistic regression (phylolm::phyloglm, Ives & Garland 2010): ndh lost ~ strategies + habitat + CAM
# 2. Pagel tests of habitat/CAM within autotrophs (do these factors alone predict ndh loss?)
# 3. carnivory vs autotrophy within habitat strata (wet only; not wet; excluding aquatic and CAM taxa)
suppressMessages({library(ape); library(phytools); library(phylolm)})
a <- commandArgs(TRUE); part <- if (length(a) >= 4) as.integer(a[4]) else 0; ch <- read.tree(a[1]); tv <- read.csv(a[2], stringsAsFactors = FALSE); out <- a[3]
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
tv$lost <- ifelse(grepl("\\?", tv$ndh_calls), NA, as.integer(nbad >= 3))
tv <- tv[!is.na(tv$lost), ]
tv$carn <- as.integer(tv$nutrition == "carnivorous")
tv$para <- as.integer(tv$nutrition %in% c("hemiparasitic", "holoparasitic"))
tv$myco <- as.integer(tv$nutrition %in% c("full_mycoheterotroph", "partial_mycoheterotroph"))
t1 <- keep.tip(ch, tv$tip); rownames(tv) <- tv$tip; d <- tv[t1$tip.label, ]
H <- max(nodeHeights(t1)); t1$edge.length <- t1$edge.length / H   # scale tree height to 1 (phyloglm alpha bounds)

if (part == 0) {
fit <- function(f) phyloglm(f, data = d, phy = t1, method = "logistic_MPLE", btol = 30, boot = 0)
mods <- list(full = lost ~ carn + para + myco + CAM + aquatic + wet,
             no_carnivory = lost ~ para + myco + CAM + aquatic + wet,
             habitat_only = lost ~ CAM + aquatic + wet,
             strategies_only = lost ~ carn + para + myco)
F <- lapply(mods, fit)
cs <- summary(F$full)$coefficients
coef_tab <- data.frame(term = rownames(cs), estimate = cs[, 1], SE = cs[, 2], z = cs[, 3], p = cs[, 4], odds_ratio = exp(cs[, 1]))
write.csv(coef_tab, paste0(out, "_phyloglm_coefficients.csv"), row.names = FALSE)
aic <- data.frame(model = names(F), AIC = sapply(F, function(m) m$aic), logLik = sapply(F, function(m) m$logLik), alpha = sapply(F, function(m) m$alpha))
aic$dAIC <- aic$AIC - min(aic$AIC)
write.csv(aic, paste0(out, "_phyloglm_models.csv"), row.names = FALSE)
cat("== phylogenetic logistic regression (full model)\n"); print(coef_tab, digits = 3)
cat("== model comparison\n"); print(aic, digits = 4)
}

pag <- function(keep, yname, ylab, label) {
  t2 <- keep.tip(ch, keep); dd <- tv[t2$tip.label, ]
  x <- setNames(ifelse(dd$lost == 1, "lost", "intact"), dd$tip); y <- setNames(ifelse(dd[[yname]] == 1, ylab, "other"), dd$tip)
  if (length(unique(y)) < 2 || length(unique(x)) < 2) return(NULL)
  f <- fitPagel(t2, x, y, dep.var = "x"); Q <- f$dependent.Q
  data.frame(test = label, n_tips = Ntip(t2), n_factor = sum(y == ylab), n_lost_factor = sum(y == ylab & x == "lost"),
             n_lost_other = sum(y == "other" & x == "lost"), LR = 2 * (f$dependent.logL - f$independent.logL), p = f$P,
             ratio = Q[paste0("intact|", ylab), paste0("lost|", ylab)] / Q["intact|other", "lost|other"])
}
auto <- tv$tip[tv$nutrition == "autotrophic"]
ca <- tv$tip[tv$nutrition %in% c("autotrophic", "carnivorous")]
specs <- list(
  list(auto, "CAM", "CAM", "CAM vs C3, autotrophs only"),
  list(auto, "aquatic", "aquatic", "aquatic vs not, autotrophs only"),
  list(auto, "wet", "wet", "wet habitat vs not, autotrophs only"),
  list(intersect(ca, tv$tip[tv$wet == 1]), "carn", "carnivorous", "carnivory vs autotrophy, wet-habitat taxa only"),
  list(intersect(ca, tv$tip[tv$wet == 0]), "carn", "carnivorous", "carnivory vs autotrophy, non-wet taxa only"),
  list(intersect(ca, tv$tip[tv$aquatic == 0 & tv$CAM == 0]), "carn", "carnivorous", "carnivory vs autotrophy, excluding aquatic and CAM taxa"))
idx <- if (part == 1) 1:3 else if (part == 2) 4:6 else 1:6
P <- do.call(rbind, lapply(specs[idx], function(z) pag(z[[1]], z[[2]], z[[3]], z[[4]])))
write.csv(P, paste0(out, "_stratified_pagel", if (part) paste0("_part", part) else "", ".csv"), row.names = FALSE)
cat("== Pagel tests\n"); print(P, digits = 3)
cat("done\n")
