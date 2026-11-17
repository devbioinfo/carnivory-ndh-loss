suppressMessages({library(ape); library(phytools)})
ch <- read.tree("f311_GBM_t3_chronogram.nwk"); tv <- read.csv("phaseC_traits_cov311.csv", stringsAsFactors = FALSE)
tv <- tv[match(ch$tip.label, tv$tip), ]
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-"))); amb <- grepl("\\?", tv$ndh_calls)
lost <- ifelse(amb, NA, nbad >= 3)
S2 <- c("I", "L"); X2 <- matrix(0, Ntip(ch), 2, dimnames = list(ch$tip.label, S2))
X2[which(lost %in% FALSE), "I"] <- 1; X2[which(lost %in% TRUE), "L"] <- 1; X2[is.na(lost), ] <- 0.5
S4 <- c("I.a", "I.b", "L.a", "L.b"); X4 <- matrix(0, Ntip(ch), 4, dimnames = list(ch$tip.label, S4))
X4[, 1:2] <- X2[, "I"]; X4[, 3:4] <- X2[, "L"]
m <- function(v, S) matrix(v, length(S), length(S), byrow = TRUE, dimnames = list(S, S))
mods <- list(
 "Mk irreversible"            = list(X2, m(c(0,1, 0,0), S2), c(1,0)),
 "Mk reversible (ARD)"        = list(X2, m(c(0,1, 2,0), S2), c(1,0)),
 # hidden rate categories a (slow) / b (fast); category switching rate r (shared)
 "HRM irreversible, 2 rates"  = list(X4, m(c(0,3,1,0,  3,0,0,2,  0,0,0,3,  0,0,3,0), S4), c(.5,.5,0,0)),
 "HRM reversible, 2 rates"    = list(X4, m(c(0,5,1,0,  5,0,0,2,  3,0,0,5,  0,4,5,0), S4), c(.5,.5,0,0)))
res <- do.call(rbind, lapply(names(mods), function(nm) {
  z <- mods[[nm]]; set.seed(1)
  best <- NULL
  for (i in 1:5) { f <- tryCatch(fitMk(ch, z[[1]], model = z[[2]], pi = z[[3]], opt.method = "nlminb", rand_start = i > 1), error = function(e) NULL)
                   if (!is.null(f) && (is.null(best) || logLik(f) > logLik(best))) best <- f }
  r <- best$rates; cat(nm, ":", signif(r, 3), "\n")
  data.frame(model = nm, k = length(r), logLik = as.numeric(logLik(best)), AIC = AIC(best),
             gain_rate = if (grepl("reversible \\(ARD\\)", nm)) r[2] else if (nm == "HRM reversible, 2 rates") paste(signif(r[3:4], 3), collapse = "/") else 0,
             loss_rates = paste(signif(r[1:min(2, length(r))][if (grepl("^Mk", nm)) 1 else 1:2], 3), collapse = "/"))
}))
res$dAIC <- res$AIC - min(res$AIC); res$w <- exp(-0.5 * res$dAIC) / sum(exp(-0.5 * res$dAIC))
print(res, row.names = FALSE); write.csv(res, "hrm_model_comparison.csv", row.names = FALSE)
