# FINAL four-panel regeneration with recovered Col1 data.
# Col1 MF vs Colilert (n=153) = Dec 2025 dedicated study (recovered Google Sheet
#   'v1.2_results') + Jan 2026 bench bucket (paper_clean_bucket.rds).
# Col2/3/4 from the paper's RDS records (reproduce the published panels).
suppressMessages({library(tidyverse); library(patchwork); library(ggthemes); library(yardstick)})
PAPER<-"data"; VALID<-"data"; SP<-"data"
OUT<-file.path("output","three_way_comparison.png")

bucket<-readRDS(file.path(PAPER,"paper_clean_bucket.rds"))
stream<-readRDS(file.path(PAPER,"paper_clean_stream.rds"))
errors<-read_csv(file.path(PAPER,"quantitray_error_intervals.csv"),show_col_types=FALSE) %>%
  group_by(mpn) %>% summarize(lower_confidence_95=min(lower_confidence_95),upper_confidence_95=max(upper_confidence_95)) %>% arrange(mpn)
fit_lume<-function(df,f) df %>% group_by(barcode) %>% nest() %>%
  mutate(pd=map(data,~mutate(.x,pred=10^predict(lm(f,data=.x),newdata=.x)))) %>% select(barcode,pd) %>% unnest(pd) %>% ungroup()
be<-fit_lume(bucket,log10(colilert)~mon2_val*floor_temp*tof_mean)
se<-fit_lume(stream,log10(colilert)~mon2_val*floor_temp+tof_mean)
ab<-function(df) df %>% left_join(errors%>%select(mpn_low=mpn,lower_confidence_95),by=join_by(closest(min_e_colilert>=mpn_low))) %>%
  left_join(errors%>%select(mpn_high=mpn,upper_confidence_95),by=join_by(closest(max_e_colilert<=mpn_high)))
be<-ab(be); se<-ab(se)

# ---- Col1: recovered Dec study + Jan bench ----
# v12_results.csv (December 2025 MF/Colilert study) is not distributed; without it Col1 holds the
# January chamber pairs only. Its paired values are in data/validation_panels.json (col1).
dec<-if(file.exists(file.path(SP,"v12_results.csv"))) {
  read_csv(file.path(SP,"v12_results.csv"),show_col_types=FALSE) %>%
  filter(method%in%c("colilert","mf"),!is.na(ecoli)) %>%
  group_by(datetime,method) %>% summarise(m=mean(ecoli),lo=min(ecoli),hi=max(ecoli),.groups="drop") %>%
  pivot_wider(names_from=method,values_from=c(m,lo,hi)) %>%
  filter(!is.na(m_colilert),!is.na(m_mf)) %>%
  transmute(colilert=m_colilert,mf=m_mf,min_e_colilert=lo_colilert,max_e_colilert=hi_colilert,min_e_mf=lo_mf,max_e_mf=hi_mf)
} else tibble(colilert=numeric(),mf=numeric(),min_e_colilert=numeric(),max_e_colilert=numeric(),min_e_mf=numeric(),max_e_mf=numeric())
jan<-bucket %>% distinct(round_minute,.keep_all=TRUE) %>% filter(!is.na(colilert),!is.na(mf)) %>%
  transmute(colilert,mf,min_e_colilert,max_e_colilert,min_e_mf,max_e_mf)
mf_vs_col<-bind_rows(dec,jan) %>% filter(colilert>0,mf>0) %>%
  transmute(x=colilert,y=mf,x_lo=min_e_colilert,x_hi=max_e_colilert,y_lo=min_e_mf,y_hi=max_e_mf,
            is_in=abs(log10(mf/colilert))<=log10(3))

lume_vs_col<-bind_rows(
  be%>%transmute(x=colilert,y=pred,x_lo=min_e_colilert,x_hi=max_e_colilert,y_lo=NA_real_,y_hi=NA_real_,is_in=pred>=lower_confidence_95&pred<=upper_confidence_95),
  se%>%transmute(x=colilert,y=pred,x_lo=min_e_colilert,x_hi=max_e_colilert,y_lo=NA_real_,y_hi=NA_real_,is_in=pred>=lower_confidence_95&pred<=upper_confidence_95))
lume_vs_mf_colcal<-be%>%filter(!is.na(mf))%>%transmute(x=mf,y=pred,x_lo=min_e_mf,x_hi=max_e_mf,y_lo=NA_real_,y_hi=NA_real_,is_in=pred>=min_e_mf&pred<=max_e_mf)
# Col4: per-unit MF calibration on the paper's own record (chamber 173 + field 33 = 206), the same
# model forms as the Colilert calibration, non-detects at 0.1. This is the published MF panel
# (Table 2: 206 / 0.906, 184 / 0.881). Corrected 2026-11-18; an earlier version fitted
# bucket_combo_extended.csv instead, which is not the calibration behind the published panel.
bmf<-bucket%>%mutate(barcode=as.character(barcode))%>%filter(!is.na(mf))%>%mutate(mf2=pmax(mf,0.1))
smf<-stream%>%mutate(barcode=as.character(barcode))%>%filter(!is.na(mf))%>%mutate(mf2=pmax(mf,0.1))
mfc<-bind_rows(fit_lume(bmf,log10(mf2)~mon2_val*floor_temp*tof_mean),fit_lume(smf,log10(mf2)~mon2_val*floor_temp+tof_mean))
lume_vs_mf_mfcal<-mfc%>%transmute(x=mf,y=pred,x_lo=min_e_mf,x_hi=max_e_mf,y_lo=NA_real_,y_hi=NA_real_,is_in=pred>=coalesce(min_e_mf,mf)&pred<=coalesce(max_e_mf,mf))

bin<-function(v)factor(cut(v,breaks=c(-Inf,10,100,Inf),labels=c("<10","10-100",">100")),levels=c("<10","10-100",">100"),ordered=TRUE)
st<-function(nm,d){d<-d%>%filter(x>0,y>0,!is.na(x),!is.na(y));f<-lm(log10(d$y)~log10(d$x));dl<-log10(d$y)-log10(d$x)
  dd<-d%>%mutate(o=bin(x),p=bin(y))
  cat(sprintf("%-24s n=%3d R2=%.3f bias=%+.2f LoA[%.2f,%.2f] Acc=%.2f BalAcc=%.2f K=%.2f slope=%.3f\n",
    nm,nrow(d),summary(f)$r.squared,mean(dl),mean(dl)-1.96*sd(dl),mean(dl)+1.96*sd(dl),
    accuracy(dd,o,p)$.estimate,bal_accuracy(dd,o,p)$.estimate,kap(dd,o,p)$.estimate,coef(f)[2]))}
cat("=== FINAL 4-PANEL STATS ===\n")
st("Col1 MF vs Colilert",mf_vs_col); st("Col2 Lume vs Colilert",lume_vs_col)
st("Col3 Lume vs MF(Ccal)",lume_vs_mf_colcal); st("Col4 Lume vs MF(MFcal)",lume_vs_mf_mfcal)

reg<-function(p){d<-p$data%>%filter(x>0,y>0);f<-lm(log10(d$y)~log10(d$x))
  lab<-sprintf("R² = %.3f, %d%% within range, n = %d",summary(f)$r.squared,round(mean(d$is_in,na.rm=TRUE)*100),nrow(d))
  rng<-range(c(d$x,d$y,d$x_lo,d$x_hi,d$y_lo,d$y_hi),na.rm=TRUE);rng[1]<-max(rng[1],0.1)
  ggplot(d,aes(x,y))+{if(any(!is.na(d$x_lo)))geom_errorbarh(aes(xmin=x_lo,xmax=x_hi),height=0,color="gray70",alpha=.5)}+
    {if(any(!is.na(d$y_lo)))geom_errorbar(aes(ymin=y_lo,ymax=y_hi),width=0,color="gray70",alpha=.5)}+
    geom_abline(slope=1,intercept=0,color="red",linetype="dashed")+geom_smooth(method="lm",formula=y~x,se=FALSE,color="#1b6e5c",linewidth=.9)+
    geom_point(aes(color=is_in),alpha=.85,size=1.6)+scale_color_manual(values=c("TRUE"="#1b9e77","FALSE"="#e08a3c"),labels=c("TRUE"="Range overlap","FALSE"="Non-overlapping"),name=NULL,na.translate=FALSE)+
    scale_x_log10(limits=rng,labels=scales::comma)+scale_y_log10(limits=rng,labels=scales::comma)+coord_fixed()+
    labs(title=p$title,subtitle=lab,x=p$xlab,y=p$ylab)+theme_few(base_size=10)+
    theme(plot.subtitle=element_text(size=8),legend.position=c(.98,.02),legend.justification=c(1,0),legend.background=element_rect(fill=alpha("white",.7),color=NA),legend.key.size=unit(.3,"cm"),legend.text=element_text(size=7))}
ba<-function(p){d<-p$data%>%filter(x>0,y>0)%>%mutate(ml=(log10(x)+log10(y))/2,dl=log10(y)-log10(x));b<-mean(d$dl);s<-sd(d$dl);lo<-b-1.96*s;hi<-b+1.96*s
  ggplot(d,aes(10^ml,dl))+annotate("rect",xmin=-Inf,xmax=Inf,ymin=lo,ymax=hi,fill="#1b9e77",alpha=.15)+
    geom_hline(yintercept=0,color="red",linetype="dashed")+geom_hline(yintercept=b,color="#1b6e5c",linewidth=.6)+geom_hline(yintercept=c(lo,hi),color="#1b6e5c",linetype="dotted",linewidth=.6)+
    geom_point(color="#1b6e5c",alpha=.75,size=1.4)+scale_x_log10(labels=scales::comma)+scale_y_continuous(limits=c(-1.75,1.75))+
    labs(subtitle=sprintf("Bias = %+.2f log₁₀  |  LoA [%.2f, %.2f]",b,lo,hi),x="Mean (CFU/100 mL)",y="Difference (log₁₀)")+theme_few(base_size=10)+theme(plot.subtitle=element_text(size=8))}
cm<-function(p){d<-p$data%>%filter(!is.na(x),!is.na(y))%>%mutate(obs_bin=bin(x),pred_bin=bin(y));cmx<-d%>%count(obs_bin,pred_bin,.drop=FALSE)
  lab<-sprintf("Acc: %.2f  |  Bal Acc: %.2f  |  Kappa: %.2f",accuracy(d,obs_bin,pred_bin)$.estimate,bal_accuracy(d,obs_bin,pred_bin)$.estimate,kap(d,obs_bin,pred_bin)$.estimate)
  ggplot(cmx,aes(obs_bin,pred_bin,fill=n))+geom_tile(color="white",linewidth=1)+geom_text(aes(label=n),color="white",fontface="bold",size=4.2)+
    scale_fill_viridis_c(begin=.15,end=.85,option="magma",guide="none")+scale_y_discrete(limits=rev)+coord_fixed()+
    labs(subtitle=lab,x=p$xlab,y=p$ylab)+theme_few(base_size=10)+theme(plot.subtitle=element_text(size=8),panel.grid=element_blank())}
pan<-list(list(data=mf_vs_col,xlab="Colilert (MPN/100 mL)",ylab="MF (CFU/100 mL)",title="MF vs. Colilert"),
  list(data=lume_vs_col,xlab="Colilert (MPN/100 mL)",ylab="Lume (CFU/100 mL)",title="Lume vs. Colilert"),
  list(data=lume_vs_mf_colcal,xlab="MF (CFU/100 mL)",ylab="Lume (CFU/100 mL)",title="Lume vs. MF (Colilert-trained)"),
  list(data=lume_vs_mf_mfcal,xlab="MF (CFU/100 mL)",ylab="Lume (CFU/100 mL)",title="Lume vs. MF (MF-trained)"))
fig<-wrap_plots(lapply(pan,reg),nrow=1)/wrap_plots(lapply(pan,ba),nrow=1)/wrap_plots(lapply(pan,cm),nrow=1)
ggsave(OUT,plot=fig,width=16,height=12,dpi=150)
cat("saved",OUT,"\n")
