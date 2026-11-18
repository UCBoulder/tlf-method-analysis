# Reviewer 1 comment 5 (ew-2026-014215): the three grab-and-read IA values in Table 2 come from
# different sample sets (190 / 184 / 134). This scores all three comparisons on ONE set: the
# chamber (bucket) unit readings where Colilert, MF and both per-unit calibrations exist, with every
# value in every comparison at or above the 10 CFU/100 mL LOQ.
# Calibrations are those behind Table 2: per unit, log10(reference) ~ mon2_val*floor_temp*tof_mean
# on the paper's chamber record, fitted once to Colilert and once to MF (non-detects at 0.1).
# Corrected 2026-11-18: an earlier version took the MF calibration from bucket_combo_extended.csv,
# which is not the calibration behind the published MF panel.
suppressMessages(library(tidyverse))
PAPER <- Sys.getenv("PAPER_DATA", "data")
fit_lume<-function(df,f,y) df %>% group_by(barcode) %>% nest() %>%
  mutate(pd=map(data,~mutate(.x,!!y:=10^predict(lm(f,data=.x),newdata=.x)))) %>% select(barcode,pd) %>% unnest(pd) %>% ungroup()
ia<-function(o,p){lo<-log10(o);lp<-log10(p);m<-mean(lo)
  c(n=length(o),IA=1-sum((lo-lp)^2)/sum((abs(lp-m)+abs(lo-m))^2),R2=cor(lo,lp)^2,RMSE=sqrt(mean((lo-lp)^2)),sdref=sd(lo))}
b<-readRDS(file.path(PAPER,"paper_clean_bucket.rds"))%>%mutate(barcode=as.character(barcode))
pc<-fit_lume(b,log10(colilert)~mon2_val*floor_temp*tof_mean,"pred_col")
pm<-fit_lume(b%>%filter(!is.na(mf))%>%mutate(mf2=pmax(mf,0.1)),log10(mf2)~mon2_val*floor_temp*tof_mean,"pred_mf")%>%
  select(round_minute,barcode,pred_mf)
j<-pc%>%inner_join(pm,by=c("round_minute","barcode"))
L<-10; cs<-j%>%filter(colilert>=L,mf>=L,pred_col>=L,pred_mf>=L)
cat(sprintf("COMMON SET: %d unit readings, %d samples, %d units\n",nrow(cs),n_distinct(cs$round_minute),n_distinct(cs$barcode)))
r<-rbind(`Colilert vs candidate (Colilert-calibrated)`=ia(cs$colilert,cs$pred_col),
         `MF vs candidate (MF-calibrated)`=ia(cs$mf,cs$pred_mf),
         `Colilert vs MF (unit readings)`=ia(cs$colilert,cs$mf))
d<-cs%>%distinct(round_minute,.keep_all=TRUE)
r<-rbind(r,`Colilert vs MF (distinct samples)`=ia(d$colilert,d$mf))
print(round(r,3))
