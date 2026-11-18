# Reference-uncertainty propagation for the ES&T Water manuscript, Section 3.2.
#
# Reuses the grab-and-read pipeline of scripts/generate_three_way_comparison.R:
# per-barcode linear calibration on the instrument's own channels, Colilert reference,
# reproducing n = 209 exactly. Reference intervals are the tabulated IDEXX
# Quanti-Tray/2000 95% confidence limits (quantitray_error_intervals.csv).
#
# Reports, at a given operating point:
#   point-label balanced accuracy;
#   expected balanced accuracy, weighting each sample by P(true >= limit) implied
#     by the reference's own interval (log10-normal, sd from the tabulated 95% limits);
#   the reference-tolerant ceiling, crediting a prediction whenever the reference
#     interval straddles the limit and so cannot resolve the category.
#
# Comparison is made at the DECISION CUT and labeled against the ACTION LIMIT:
# 82 CFU reproduces the canonical operating point (sens 0.941 / spec 0.859 / BA 0.900).
#
# Run: Rscript scripts/reference_uncertainty_propagation.R

suppressMessages({library(tidyverse)})
PAPER <- "data"; LIMIT <- 126
bucket <- readRDS(file.path(PAPER,"paper_clean_bucket.rds"))
stream <- readRDS(file.path(PAPER,"paper_clean_stream.rds"))
errors <- read_csv(file.path(PAPER,"quantitray_error_intervals.csv"), show_col_types=FALSE) %>%
  group_by(mpn) %>% summarize(lo=min(lower_confidence_95), hi=max(upper_confidence_95), .groups="drop") %>% arrange(mpn)
fit_lume <- function(df,f) df %>% group_by(barcode) %>% nest() %>%
  mutate(pd=map(data,~mutate(.x, pred=10^predict(lm(f,data=.x), newdata=.x)))) %>% select(barcode,pd) %>% unnest(pd) %>% ungroup()
ab <- function(df) df %>%
  left_join(errors %>% select(mpn_low=mpn, lo), by=join_by(closest(min_e_colilert>=mpn_low))) %>%
  left_join(errors %>% select(mpn_high=mpn, hi), by=join_by(closest(max_e_colilert<=mpn_high)))
d <- bind_rows(ab(fit_lume(bucket, log10(colilert)~mon2_val*floor_temp*tof_mean)),
               ab(fit_lume(stream, log10(colilert)~mon2_val*floor_temp+tof_mean))) %>%
     filter(colilert>0, pred>0, !is.na(lo), !is.na(hi)) %>%
     mutate(obs = colilert>=LIMIT,
            sd_log = (log10(pmax(hi,1e-9))-log10(pmax(lo,1e-9)))/(2*1.959964),
            p_exc  = 1-pnorm((log10(LIMIT)-log10(colilert))/sd_log),
            unresolved = lo<LIMIT & hi>=LIMIT)

cat("finding the published operating point (target sens 0.94 / spec 0.86):\n")
for (cut in seq(76,92,2)) { p <- d$pred>=cut
  cat(sprintf("  cut %2d: sens %.3f spec %.3f BA %.3f\n", cut, mean(p[d$obs]), mean(!p[!d$obs]), (mean(p[d$obs])+mean(!p[!d$obs]))/2)) }

report <- function(cut){
  p <- d$pred>=cut
  ba <- (mean(p[d$obs])+mean(!p[!d$obs]))/2
  se <- sum(d$p_exc*p)/sum(d$p_exc); sp <- sum((1-d$p_exc)*!p)/sum(1-d$p_exc)
  tol <- (p==d$obs)|d$unresolved
  cat(sprintf("\n--- operating point %d CFU (n=%d) ---\n", cut, nrow(d)))
  cat(sprintf("point-label    sens %.3f spec %.3f BA %.4f\n", mean(p[d$obs]), mean(!p[!d$obs]), ba))
  cat(sprintf("expected       sens %.3f spec %.3f BA %.4f   shift %+.4f\n", se, sp, (se+sp)/2, (se+sp)/2-ba))
  cat(sprintf("ref-tolerant   sens %.3f spec %.3f BA %.4f   ceiling\n",
      mean(tol[d$obs]), mean(tol[!d$obs]), (mean(tol[d$obs])+mean(tol[!d$obs]))/2))
}
report(82); report(126)
cat(sprintf("\nunresolved at 126 (tabulated 95%% CI straddles the limit): %.1f%% (%d/%d)\n",
    100*mean(d$unresolved), sum(d$unresolved), nrow(d)))
cat(sprintf("prediction inside the reference 95%% interval: %.1f%%\n", 100*mean(d$pred>=d$lo & d$pred<=d$hi)))
