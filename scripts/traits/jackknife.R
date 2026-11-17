# Leave-one-lineage-out Pagel test (carnivory vs autotrophy), to check that no single carnivorous lineage
# drives the association (Maddison & FitzJohn 2015 pseudoreplication concern).
suppressMessages({library(ape); library(phytools)})
a <- commandArgs(TRUE); ch <- read.tree(a[1]); tv <- read.csv(a[2], stringsAsFactors = FALSE); thr <- 3
nbad <- sapply(strsplit(tv$ndh_calls, ""), function(z) sum(z %in% c("P", "-")))
tv$ndh <- ifelse(grepl("\\?", tv$ndh_calls), NA, ifelse(nbad >= thr, "lost", "intact"))
lin <- list(none = character(0), Lentibulariaceae = c("Pinguicula", "Utricularia", "Genlisea"),
            Utricularia = "Utricularia", Droseraceae = c("Drosera", "Dionaea", "Aldrovanda"),
            Drosophyllum = "Drosophyllum", Nepenthes = "Nepenthes", Sarraceniaceae = c("Sarracenia", "Heliamphora", "Darlingtonia"),
            Triantha = "Triantha", Brocchinia = "Brocchinia")
sel <- if (length(a) >= 4) strsplit(a[4], ",")[[1]] else names(lin)
out <- list()
for (L in sel) {
  drop <- tv$tip[Reduce(`|`, lapply(c(lin[[L]], "ZZZ"), function(g) startsWith(tv$tip, g)))]
  keep <- tv$tip[!is.na(tv$ndh) & tv$nutrition %in% c("autotrophic", "carnivorous") & !(tv$tip %in% drop)]
  t2 <- keep.tip(ch, keep); d <- tv[match(t2$tip.label, tv$tip), ]
  x <- setNames(d$ndh, d$tip); y <- setNames(ifelse(d$nutrition == "carnivorous", "carnivorous", "autotrophic"), d$tip)
  f <- fitPagel(t2, x, y, dep.var = "x"); Q <- f$dependent.Q
  out[[L]] <- data.frame(dropped = L, n_tips = Ntip(t2), n_carn = sum(y == "carnivorous"), n_carn_lost = sum(y == "carnivorous" & x == "lost"),
    LR = 2 * (f$dependent.logL - f$independent.logL), p = f$P,
    ratio = Q["intact|carnivorous", "lost|carnivorous"] / Q["intact|autotrophic", "lost|autotrophic"])
  print(out[[L]]); flush.console()
}
write.csv(do.call(rbind, out), a[3], row.names = FALSE)
