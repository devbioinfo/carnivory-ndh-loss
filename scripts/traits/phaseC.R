# Phase C: ndh loss across nutritional strategies on a dated tree.
# usage: Rscript phaseC.R <full_tree.nwk> <calib.csv> <phaseC_traits.csv> <out_prefix> [nsim=200]
#  1. full 295-taxon chronogram: penalised likelihood (ape::chronos, correlated rates), with every node of the
#     MCMCTree dating subset fixed at its posterior mean age (secondary calibrations)
#  2. ndh state (intact = all 11 ndh genes functional; lost = any pseudogenised/absent; Byblis ambiguous):
#     irreversible-loss Mk model vs reversible (AIC), stochastic maps under the irreversible model, root intact
#  3. nutritional strategy: star model (autotrophy <-> each heterotrophic strategy), root autotrophic, stochastic maps
#  4. paired maps -> ndh loss events, the regime in force at each loss, exposure time per regime while ndh intact,
#     loss rate per Myr per regime; timing of ndh loss vs origin of carnivory per carnivorous lineage
suppressMessages({library(ape); library(phytools)})
a <- commandArgs(TRUE)
nsim <- if (length(a) >= 5) as.integer(a[5]) else 200
thr <- if (length(a) >= 6) as.integer(a[6]) else 3   # ndh "lost" = at least thr of 11 genes pseudogenised/absent
set.seed(20261116)
out <- a[4]

## 1. chronogram ------------------------------------------------------------
tr <- read.tree(a[1]); tr$node.label <- NULL
tr <- root(tr, "Amborella_trichopoda", resolve.root = TRUE)
tr$edge.length[tr$edge.length <= 0] <- 1e-6
cal <- read.csv(a[2])
nd <- mapply(function(x, y) getMRCA(tr, c(x, y)), cal$tipA, cal$tipB)
calib <- data.frame(node = nd, age.min = cal$age_mean, age.max = cal$age_mean, soft.bounds = FALSE)
calib <- calib[!duplicated(calib$node), ]
cf <- paste0(out, "_chronogram.nwk")
if (file.exists(cf)) { ch <- read.tree(cf) } else {
ch <- chronos(tr, lambda = 1, model = "correlated", calibration = calib, quiet = TRUE,
              control = chronos.control(iter.max = 1e4, eval.max = 1e4))
class(ch) <- "phylo"
ch$edge.length[ch$edge.length <= 0] <- 1e-4
write.tree(ch, cf) }
H <- max(nodeHeights(ch))
cat(sprintf("chronogram: %d tips, root age %.2f Ma, %d secondary calibrations\n", Ntip(ch), H, nrow(calib)))

## 2. ndh state ----------------------------------------------------------------
tv <- read.csv(a[3], stringsAsFactors = FALSE)
tv <- tv[match(ch$tip.label, tv$tip), ]
stopifnot(!any(is.na(tv$tip)))
X <- matrix(0, Ntip(ch), 2, dimnames = list(ch$tip.label, c("intact", "lost")))
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
amb <- grepl("\\?", tv$ndh_calls)
X[!amb & nbad < thr, "intact"] <- 1
X[!amb & nbad >= thr, "lost"] <- 1
X[amb, ] <- 0.5
cat(sprintf("ndh threshold %d: %d intact, %d lost, %d ambiguous\n", thr, sum(X[, 1] == 1), sum(X[, 2] == 1), sum(amb)))
irr <- matrix(c(0, 1, 0, 0), 2, 2, byrow = TRUE, dimnames = list(c("intact", "lost"), c("intact", "lost")))
f_irr <- fitMk(ch, X, model = irr, pi = c(1, 0))
f_ard <- fitMk(ch, X, model = "ARD", pi = c(1, 0))
mc <- data.frame(model = c("irreversible loss", "reversible (ARD)"),
                 logLik = c(logLik(f_irr), logLik(f_ard)), k = c(1, 2),
                 AIC = c(AIC(f_irr), AIC(f_ard)),
                 loss_rate = c(f_irr$rates[1], f_ard$rates[which(f_ard$index.matrix[1, 2] == 1:2)]))
write.csv(mc, paste0(out, "_ndh_model_comparison.csv"), row.names = FALSE)
print(mc)
smN <- make.simmap(ch, X, Q = as.Qmatrix(f_irr), pi = c(1, 0), nsim = nsim, message = FALSE)

## 3. nutritional strategy ----------------------------------------------------
st <- sort(unique(tv$nutrition)); st <- c("autotrophic", setdiff(st, "autotrophic"))
Y <- to.matrix(setNames(tv$nutrition, tv$tip), st)
star <- matrix(0, length(st), length(st), dimnames = list(st, st)); k <- 0
for (i in 2:length(st)) { k <- k + 1; star[1, i] <- k; k <- k + 1; star[i, 1] <- k }
f_nut <- fitMk(ch, Y, model = star, pi = c(1, rep(0, length(st) - 1)))
smT <- make.simmap(ch, Y, Q = as.Qmatrix(f_nut), pi = c(1, rep(0, length(st) - 1)), nsim = nsim, message = FALSE)

## 4. paired maps: events, exposure, rates, timing -------------------------------
nh <- nodeHeights(ch)
seg <- function(m) {  # list per edge of data.frame(state, start, end) in time from root
  lapply(seq_len(nrow(ch$edge)), function(e) {
    v <- m$maps[[e]]; s <- nh[e, 1] + c(0, cumsum(v)[-length(v)])
    data.frame(state = names(v), start = s, end = s + v, stringsAsFactors = FALSE) })
}
state_at <- function(sg, t) { i <- which(sg$start <= t + 1e-9 & sg$end >= t - 1e-9)[1]; sg$state[i] }
rates <- list(); events <- list()
for (s in seq_len(nsim)) {
  sn <- seg(smN[[s]]); stt <- seg(smT[[s]])
  expo <- setNames(numeric(length(st)), st); loss <- expo
  for (e in seq_along(sn)) {
    A <- sn[[e]]; B <- stt[[e]]
    for (i in seq_len(nrow(A))) {         # exposure: time intact, split by regime
      if (A$state[i] != "intact") next
      for (j in seq_len(nrow(B))) {
        ov <- min(A$end[i], B$end[j]) - max(A$start[i], B$start[j])
        if (ov > 0) expo[B$state[j]] <- expo[B$state[j]] + ov
      }
    }
    if (nrow(A) > 1) for (i in 2:nrow(A)) if (A$state[i - 1] == "intact" && A$state[i] == "lost") {
      t <- A$start[i]; r <- state_at(B, t); loss[r] <- loss[r] + 1
      events[[length(events) + 1]] <- data.frame(sim = s, edge = e, age_Ma = H - t, regime = r,
        tips_below = { d <- ch$edge[e, 2]; tl <- if (d <= Ntip(ch)) ch$tip.label[d] else extract.clade(ch, d)$tip.label; paste0(length(tl), ":", paste(head(sort(tl), 3), collapse = ";")) })
    }
  }
  rates[[s]] <- data.frame(sim = s, regime = st, losses = loss, exposure_Myr = expo, rate_per_Myr = ifelse(expo >= 1, loss / expo, NA))
}
R <- do.call(rbind, rates); E <- do.call(rbind, events)
write.csv(R, paste0(out, "_rates_per_map.csv"), row.names = FALSE)
write.csv(E, paste0(out, "_loss_events.csv"), row.names = FALSE)
q <- function(x) sprintf("%.4f (%.4f-%.4f)", median(x, na.rm = TRUE), quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
S <- do.call(rbind, lapply(st, function(r) { x <- R[R$regime == r, ]
  data.frame(regime = r, losses_median = median(x$losses), exposure_Myr_median = round(median(x$exposure_Myr), 1),
             rate_per_Myr = q(x$rate_per_Myr),
             ratio_vs_autotrophic = q(x$rate_per_Myr / R$rate_per_Myr[R$regime == "autotrophic"])) }))
write.csv(S, paste0(out, "_rates_summary.csv"), row.names = FALSE)
print(S)

## timing of ndh loss vs origin of carnivory, per carnivorous lineage (MRCA of the listed genera)
lin <- list(Lentibulariaceae = c("Pinguicula", "Utricularia", "Genlisea"), Byblis = "Byblis",
            Caryophyllales_carnivores = c("Drosera", "Nepenthes", "Drosophyllum", "Triphyophyllum", "Dionaea", "Aldrovanda"),
            Sarraceniaceae = c("Sarracenia", "Heliamphora", "Darlingtonia"), Roridula = "Roridula", Cephalotus = "Cephalotus",
            Brocchinia = c("Brocchinia_reducta", "Brocchinia_hechtioides"), Triantha = "Triantha")
pick <- function(p) ch$tip.label[Reduce(`|`, lapply(p, function(g) startsWith(ch$tip.label, g)))]
TM <- list()
for (L in names(lin)) {
  tips <- pick(lin[[L]]); if (!length(tips)) next
  node <- if (length(tips) > 1) getMRCA(ch, tips) else which(ch$tip.label == tips)
  desc_edges <- which(ch$edge[, 2] %in% c(node, getDescendants(ch, node)))
  path <- c(); n <- node
  while (n != Ntip(ch) + 1) { e <- which(ch$edge[, 2] == n); path <- c(e, path); n <- ch$edge[e, 1] }
  for (s in seq_len(nsim)) {
    stt <- seg(smT[[s]]); sn <- seg(smN[[s]])
    carn_t <- NA   # earliest entry into carnivory along the path to and within the lineage
    for (e in c(path, desc_edges)) { B <- stt[[e]]; i <- which(B$state == "carnivorous")
      if (length(i)) carn_t <- min(carn_t, B$start[i[1]], na.rm = TRUE) }
    loss_t <- NA   # earliest ndh loss on the stem path or within the lineage
    for (e in c(path, desc_edges)) { A <- sn[[e]]; i <- which(A$state == "lost")
      if (length(i) && (i[1] > 1 || (e %in% desc_edges))) if (i[1] > 1) loss_t <- min(loss_t, A$start[i[1]], na.rm = TRUE) }
    TM[[length(TM) + 1]] <- data.frame(lineage = L, sim = s, carnivory_origin_Ma = H - carn_t, first_ndh_loss_Ma = H - loss_t)
  }
}
T <- do.call(rbind, TM); write.csv(T, paste0(out, "_timing_per_map.csv"), row.names = FALSE)
TS <- do.call(rbind, lapply(split(T, T$lineage), function(x) data.frame(lineage = x$lineage[1],
  carnivory_origin_Ma = q(x$carnivory_origin_Ma), first_ndh_loss_Ma = q(x$first_ndh_loss_Ma),
  p_loss_after_carnivory = round(mean(x$first_ndh_loss_Ma < x$carnivory_origin_Ma, na.rm = TRUE), 3),
  p_no_loss = round(mean(is.na(x$first_ndh_loss_Ma)), 3))))
write.csv(TS, paste0(out, "_timing_summary.csv"), row.names = FALSE)
print(TS)
cat("done\n")
