### Load required packages ###

library(dplyr)
library(ggplot2)
library(irr)
library(vcd)
library(epibasix)
library(Epi)

### Creating a dataset with raters of interest ###

irr_lisate<-df_lisate%>%
  select(rater_tif, rater_cdic)

irr_recomb<-df_recomb%>%
  select(rater_tif, rater_cdic)

irr_iics<-df_iics%>%
  select(rater_tif, rater_cdic)

### Estimating Kappas ###

tab_lisate<-table(irr_lisate)

tab_lisate

tab_recomb<-table(irr_recomb)

tab_iics<-table(irr_iics)

print(tab_iics)
print(tab_lisate)
print(tab_recomb)

### Compute the kappa ###

kappa_lisate<-Kappa(tab_lisate)

kappa_recomb<-Kappa(tab_recomb)

kappa_iics<- Kappa(tab_iics)

print(kappa_iics)
print(kappa_recomb)
print(kappa_lisate)

### Check confidence intervals ###

confint(kappa_lisate)

confint(kappa_recomb)

confint(kappa_iics)

### Agreement polots ###

agreementplot(tab_iics)

### Determining prevalence in the cohort ###

## Prevalence with RDT TIF + ELISA CDIC + IHA when available
### This is the prevalence used in the end for prevalence calculation#


df_master%>%
  group_by(def_elisa_cdic_rdt_iha)%>%
  summarise(preval_total = n())

(26/(269+26))*100

### Prevalence of 8.813% and 5 indeterminate participants ##

metadata<-merge(metadata, df_master, by= "sample")

metadata%>%
  group_by(dept)%>%
  summarise(n())

### Performance and agreement between different tests ###

### The filtered dataframe excluding IND for RDT+ELISA-CDIC- IHA is below ###
### And identified as _3 ###


df_master_filt_3<-df_master%>%
  filter(def_elisa_cdic_rdt_iha != "IND")

metadata_filt_3<-metadata%>%
  filter(def_elisa_cdic_rdt_iha != "IND")

### Performance evaluation using the subsample of 52 participants ###


## summaries ##

df_master%>%
  group_by(iha_lcsp)%>%
  summarise(n())

df_iha_filt<-df_iha%>%
  filter(iha_lcsp != "NO_SAMP")

df_master%>%
  group_by(iha_lcsp)%>%
  summarise(n())

### summaries result by test ###

df_iha_filt%>%
  group_by(def_elisa_cdic_rdt_iha)%>%
  summarise(n())


### Agreement single RDT vs RDT+ELISA_cdic+IHA

irr_1rdt_3<-df_iha_filt%>%
  select(rdt_tif, def_elisa_cdic_rdt_iha)

tab_1rdt_3<-table(irr_1rdt_3)

tab_1rdt_3

kappa_1rdt_3<-Kappa(tab_1rdt_3)

kappa_1rdt_3

confint(kappa_1rdt_3)

### Agreement single ELISA_tif vs RDT+ELISA_cdic+IHA

irr_ELISA_tif_3<-df_iha_filt%>%
  select(elisa_tif, def_elisa_cdic_rdt_iha)

tab_ELISA_tif_3<-table(irr_ELISA_tif_3)

tab_ELISA_tif_3

kappa_ELISA_tif_3<-Kappa(tab_ELISA_tif_3)

kappa_ELISA_tif_3

confint(kappa_ELISA_tif_3)

## Agreement single ELISA_cdic vs RDT+ELISA_cdic_IHA ##


irr_ELISACDIC_3<-df_iha_filt%>%
  select(elisa_cdic, def_elisa_cdic_rdt_iha)

tab_ELISACDIC_3<-table(irr_ELISACDIC_3)

tab_ELISACDIC_3

kappa_ELISACDIC_3<-Kappa(tab_ELISACDIC_3)

kappa_ELISACDIC_3

confint(kappa_ELISACDIC_3)

### Agreement elisa_rdt_tif vs RDT+ELISA_cdic+IHA

irr_algoTIF<-df_iha_filt%>%
  select(def_elisa_rdt_tif_iha, def_elisa_cdic_rdt_iha)

tab_algoTIF<-table(irr_algoTIF)

tab_algoTIF

kappa_algoTIF<-Kappa(tab_algoTIF)

kappa_algoTIF

confint(kappa_algoTIF)


## sensitivity and specificity RDT vs RDT_ELISA_cdic_TIF ##

tab_1rdt_3

sensspec_1rdt_3<-as.table(rbind(c(24,2), c(2,24)))

rdt_perf_summ_3<-sensSpec(sensspec_1rdt_3, alpha = 0.05, CL = TRUE)

summary(rdt_perf_summ_3)


## sensitivity and specificity ELISA_tif vs RDT_ELISA_cdic_TIF

tab_ELISA_tif_3

sensspec_ELISA_tif_3<-as.table(rbind(c(23,1),c(3,25)))

ELISA_tif_perf_summ_3<-sensSpec(sensspec_ELISA_tif_3, alpha = 0.05, CL=TRUE)

summary(ELISA_tif_perf_summ_3)

## sensitivity and specificity ELISA_cdic vs ELISA_ELISA_cdic_TIF

tab_ELISACDIC_3

sensspec_ELISACDIC_3<-as.table(rbind(c(26,1),c(0,25)))

ELISACDIC_perf_summ_3<-sensSpec(sensspec_ELISACDIC_3, alpha = 0.05, CL=TRUE)

summary(ELISACDIC_perf_summ_3)

## sensitivity and specificity def_elisa_rdt_TIF (algoTIF) vs RDT_ELISA_cdic_TIF ##

tab_algoTIF

sensspec_algoTIF<-as.table(rbind(c(23,1), c(3,25)))

algoTIF_perf<-sensSpec(sensspec_algoTIF, alpha = 0.05, CL = TRUE)

summary(algoTIF_perf)

## Separating based on type of ELISA ##

df_iha_lys<-df_iha_filt%>%
  filter(elisa_type=="CONVENCIONAL")

df_iha_recomb<-df_iha_filt%>%
  filter(elisa_type=="RECOMBINANTE")

df_iha_iics<-df_iha_filt%>%
  filter(elisa_type=="IICS")

## performance comparing lysate tests ##

## RDT ##

irr_1rdt_lys<-df_iha_lys%>%
  select(rdt_tif, def_elisa_cdic_rdt_iha)

tab_1rdt_lys<-table(irr_1rdt_lys)

tab_1rdt_lys

sensspec_1rdt_lys<-as.table(rbind(c(3,0), c(0,6)))

rdt_perf_summ_lys<-sensSpec(sensspec_1rdt_lys, alpha = 0.05, CL = TRUE)

summary(rdt_perf_summ_lys)

#ELISA TIF ##

irr_ELISA_tif_lys<-df_iha_lys%>%
  select(elisa_tif, def_elisa_cdic_rdt_iha)

tab_ELISA_tif_lys<-table(irr_ELISA_tif_lys)

tab_ELISA_tif_lys

sensspec_ELISA_tif_lys<-as.table(rbind(c(3,0), c(0,6)))

ELISA_tif_perf_summ_lys<-sensSpec(sensspec_ELISA_tif_lys, alpha = 0.05, CL = TRUE)

summary(ELISA_tif_perf_summ_lys)

## ELISA CDIC ##

irr_ELISACDIC_lys<-df_iha_lys%>%
  select(elisa_cdic, def_elisa_cdic_rdt_iha)

tab_ELISACDIC_lys<-table(irr_ELISACDIC_lys)

tab_ELISACDIC_lys

sensspec_ELISACDIC_lys<-as.table(rbind(c(3,0), c(0,6)))

ELISACDIC_perf_summ_lys<-sensSpec(sensspec_ELISACDIC_lys, alpha = 0.05, CL = TRUE)

summary(ELISACDIC_perf_summ_lys)

## performance comparing recombinant tests ##

## RDT ##

irr_1rdt_recomb<-df_iha_recomb%>%
  select(rdt_tif, def_elisa_cdic_rdt_iha)

tab_1rdt_recomb<-table(irr_1rdt_recomb)

tab_1rdt_recomb

sensspec_1rdt_recomb<-as.table(rbind(c(7,1), c(0,7)))

rdt_perf_summ_recomb<-sensSpec(sensspec_1rdt_recomb, alpha = 0.05, CL = TRUE)

summary(rdt_perf_summ_recomb)

## ELISA TIF ##

irr_ELISA_tif_recomb<-df_iha_recomb%>%
  select(elisa_tif, def_elisa_cdic_rdt_iha)

tab_ELISA_tif_recomb<-table(irr_ELISA_tif_recomb)

tab_ELISA_tif_recomb

sensspec_ELISA_tif_recomb<-as.table(rbind(c(7,1), c(0,7)))

ELISA_tif_perf_summ_recomb<-sensSpec(sensspec_ELISA_tif_recomb, alpha = 0.05, CL = TRUE)

summary(ELISA_tif_perf_summ_recomb)

## ELISA CDIC ##

irr_ELISACDIC_recomb<-df_iha_recomb%>%
  select(elisa_cdic, def_elisa_cdic_rdt_iha)

tab_ELISACDIC_recomb<-table(irr_ELISACDIC_recomb)

tab_ELISACDIC_recomb

sensspec_ELISACDIC_recomb<-as.table(rbind(c(7,0), c(0,8)))

ELISACDIC_perf_summ_recomb<-sensSpec(sensspec_ELISACDIC_recomb, alpha = 0.05, CL = TRUE)

summary(ELISACDIC_perf_summ_recomb)

## performance comparing recombinant iics tests ##

## RDT ##

irr_1rdt_iics<-df_iha_iics%>%
  select(rdt_tif, def_elisa_cdic_rdt_iha)

tab_1rdt_iics<-table(irr_1rdt_iics)

tab_1rdt_iics

sensspec_1rdt_iics<-as.table(rbind(c(14,1), c(2,11)))

rdt_perf_summ_iics<-sensSpec(sensspec_1rdt_iics, alpha = 0.05, CL = TRUE)

summary(rdt_perf_summ_iics)

## ELISA TIF ##

irr_ELISA_tif_iics<-df_iha_iics%>%
  select(elisa_tif, def_elisa_cdic_rdt_iha)

tab_ELISA_tif_iics<-table(irr_ELISA_tif_iics)

tab_ELISA_tif_iics

sensspec_ELISA_tif_iics<-as.table(rbind(c(13,0), c(3,12)))

ELISA_tif_perf_summ_iics<-sensSpec(sensspec_ELISA_tif_iics, alpha = 0.05, CL = TRUE)

summary(ELISA_tif_perf_summ_iics)

## ELISA CDIC ##

irr_ELISACDIC_iics<-df_iha_iics%>%
  select(elisa_cdic, def_elisa_cdic_rdt_iha)

tab_ELISACDIC_iics<-table(irr_ELISACDIC_iics)

tab_ELISACDIC_iics

sensspec_ELISACDIC_iics<-as.table(rbind(c(16,1), c(0,11)))

ELISACDIC_perf_summ_iics<-sensSpec(sensspec_ELISACDIC_iics, alpha = 0.05, CL = TRUE)

summary(ELISACDIC_perf_summ_iics)

### Mcnemar tests ###

## McNemar RDT_TIF vs RDT+ELISA_cdic_IHA

tab_1rdt_3

mcnemar.test(tab_1rdt_3, correct = FALSE)

## McNemar RDT-TIF vs ELISA_CDIC ##

irr_rdt_elisa_cdic<-df_iha_filt%>%
  select(rdt_tif, elisa_cdic)

tab_rdt_cdic<-table(irr_rdt_elisa_cdic)

mcnemar.test(tab_rdt_cdic, correct = FALSE)

# McNemar ELISA_TIF vs RDT+ELISA_cdic_IHA

tab_ELISA_tif_3

mcnemar.test(tab_ELISA_tif_3, correct = FALSE)

## McNemar ELISA_TIF vs ELISA_CDIC ##

irr_elisa_tif_cdic<-df_iha_filt%>%
  select(elisa_tif, elisa_cdic)

tab_elisa_tif_cdic<-table(irr_elisa_tif_cdic)

mcnemar.test(tab_elisa_tif_cdic, correct = FALSE)

## McNemar RDT vs ELISA TIF ##

irr_rdt_tif_elisa_tif<-df_iha_filt%>%
  select(rdt_tif, elisa_tif)

tab_rdt_tif_elisa_tif<-table(irr_rdt_tif_elisa_tif)

mcnemar.test(tab_rdt_tif_elisa_tif, correct = FALSE)

## McNemar ELISa_CDIC vs RDT+ELISA_cdic_IHA ###

tab_ELISACDIC_3

mcnemar.test(tab_ELISACDIC_3, correct = FALSE)

##McNemar algoTIF vs RDT+ELISA_cdic+IFA##

tab_algoTIF

mcnemar.test(tab_algoTIF, correct = FALSE)

### Plots ###

## Plots using ggplot ##

library(dplyr)
library(ggplot2)
library(irr)
library(vcd)
library(epibasix)
library(exact2x2)
library(ggthemes)
library(ggalluvial)

## Using RDT+ELISA_cdic+IHA as GS from the 52 subsample ###

## dataset is summary_df_3 ###

## sensitivity ##

ggplot(data = summary_df_3, aes(x= test, y= sens, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_sens_dow, ymax= ci_sens_up), width=0.1)+
  ylim(10,100)+
  ylab("Sensitivity")+
  xlab("Test")+
  theme_bw()

## specificity ##

ggplot(data = summary_df_3, aes(x= test, y= spec, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_spec_dow, ymax= ci_spec_up), width=0.1)+
  ylim(10,100)+
  ylab("Specificity")+
  xlab("Test")+
  theme_bw()


## Agreement ##

ggplot(data = summary_df_3, aes(x= test, y= kappa, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_kappa_dow, ymax= ci_kappa_up), width=0.1)+
  ylim(0.1,1)+
  ylab("Kappa coefficient")+
  xlab("Test")+
  theme_bw()

###Plotting sensitivity and specificity of RDT_TIF and algorithm TIF ###

### Sensitivity ###

ggplot(data = summary_algorithms, aes(x= test, y= sens, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_sens_dow, ymax= ci_sens_up), width=0.1)+
  ylim(50,100)+
  ylab("Sensitivity")+
  xlab("Test")+
  theme_bw()

ggsave("sens_algo.png", units = "in", width = 5, height = 4, dpi = 300)

### Specificity ###

ggplot(data = summary_algorithms, aes(x= test, y= spec, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_spec_dow, ymax= ci_spec_up), width=0.1)+
  ylim(50,100)+
  ylab("Specificity")+
  xlab("Test")+
  theme_bw()

ggsave("spec_algo.png", units = "in", width = 5, height = 4, dpi = 300)

### Plotting sensitivity and specificity of individual tests ###


## sensitivity ##

ggplot(data = ind_elisa, aes(x= test, y= sens, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_sens_dow, ymax= ci_sens_up), width=0.1)+
  ylim(10,100)+
  ylab("Sensitivity")+
  xlab("Test")+
  theme_bw()

ggsave("sens_ELISAs.png", units = "in", width = 5, height = 4, dpi = 300)

## specificity ##

ggplot(data = ind_elisa, aes(x= test, y= spec, color= test))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymin = ci_spec_dow, ymax= ci_spec_up), width=0.1)+
  ylim(10,100)+
  ylab("Specificity")+
  xlab("Test")+
  theme_bw()

ggsave("spec_ELISAs.png", units = "in", width = 5, height = 4, dpi = 300)

### Alluvial plot ###

ggplot(data = alluvial_iha,
       aes(axis1 = rdt_tif, axis2 = elisa_tif,
           axis3 = elisa_cdic, axis4 = iha_lcsp,
           axis5 = inf, y = total)) +
  scale_x_discrete(limits = c("RDT CSTIF", "ELISA CSTIF",
                              "ELISA CDIC", "IHA LCSP", 
                              "Diagnosis"), expand = c(.1, .05)) +
  geom_flow()+
  scale_y_continuous(breaks = seq(0,55, by = 5))+
  xlab("") +
  ylab("Samples")+
  geom_alluvium(aes(fill = inf), color = "black") +
  geom_stratum(alpha = 0.2, color = "black", size = 1.1) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_manual(values = c("#619CFF", "#F8766D"))+
  theme_minimal()

ggsave("alluvial_julio.png", units = "in", width = 8.5, height = 5, dpi = 300)


