##### Fig 3. Refined genomic classification of uveal melanomas
### Load packages & files
library(tidyverse)
library(ggsurvfit)
library(survival)

df_st1 <- read_tsv("data/SuppTable1.tsv")

df_prog_final <- df_st1 %>% 
  select(sampleID,
         is_BAP1,chr3CN,is_BAP1_LOH,is_BAP1_aberration,
         chr8qCN,is_8q_gain,is_1q_gain,
         is_metastasis,relapse_free_days,TCGA) %>% 
  mutate(mmmmTCGA = ifelse(is_BAP1_aberration == 1,
                           ifelse(is_1q_gain == 1 & chr8qCN >= 4, "D","C"),
                           ifelse(is_1q_gain == 1 & chr8qCN >= 4, "C",
                                  ifelse(is_1q_gain == 1 | chr8qCN > 2, "B", "A")))) %>%
  mutate(mmmmmmTCGA = ifelse(is_BAP1_aberration + is_1q_gain + is_8q_gain == 3, "D",
                             ifelse(is_BAP1_aberration + is_1q_gain + is_8q_gain == 2, "C",
                                    ifelse(is_BAP1_aberration + is_1q_gain + is_8q_gain == 1, "B", "A"))))

  


### Fig3. Genomic classifications
## TCGA
fig3c <- crossing(tibble(TCGA = c("A","B","C","D")), tibble(is_metastasis = c(0,1))) %>% 
  left_join(df_prog_final %>% count(TCGA,is_metastasis)) %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>% 
  mutate(is_metastasis = factor(is_metastasis,levels=c(1,0))) %>% 
  group_by(TCGA) %>% mutate(n_group = sum(n)) %>% mutate(prop_group = n/n_group) %>% ungroup() %>% 
  mutate(prop_group = ifelse(is_metastasis == 0, "", str_c(round(prop_group*100),"%"))) %>% 
  ggplot(aes(TCGA,n,fill=is_metastasis))+
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values = c("darkred","grey"))+
  xlab("TCGA subtype")+ylab("Proportion of patients with metastasis")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size=22))+
  ggtitle("Proportion of metastatic samples")+
  geom_text(aes(x=TCGA,y=1,label = prop_group),vjust=-0.5,size=5)

print(fig3c)

fig3e <- survfit2(Surv(relapse_free_days, is_metastasis) ~ TCGA, data = df_prog_final %>% filter(relapse_free_days>0)) %>% 
  ggsurvfit() + 
  labs(x = "Days", y="Relapse-free survival probability")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.position = "top",
        legend.text = element_text(size = 16))

print(fig3e)

logrank_tcga <- survdiff(Surv(relapse_free_days, is_metastasis) ~ TCGA, data = df_prog_final %>% filter(relapse_free_days>0))
logrank_tcga$pvalue # 0.1754545



## mmmmTCGA
# crossing(tibble(mmmmTCGA = c("A","B","C","D")), tibble(is_metastasis = c(0,1))) %>% 
#   left_join(df_prog_final %>% count(mmmmTCGA,is_metastasis)) %>% 
#   mutate(n = ifelse(is.na(n), 0, n)) %>% 
#   mutate(is_metastasis = factor(is_metastasis,levels=c(1,0))) %>% 
#   group_by(mmmmTCGA) %>% mutate(n_group = sum(n)) %>% mutate(prop_group = n/n_group) %>% ungroup() %>% 
#   mutate(prop_group = ifelse(is_metastasis == 0, "", str_c(round(prop_group*100),"%"))) %>% 
#   ggplot(aes(mmmmTCGA,n,fill=is_metastasis))+
#   geom_bar(stat="identity", position=position_fill())+
#   scale_fill_manual(values = c("darkred","grey"))+
#   xlab("mmmmTCGA subtype")+ylab("Proportion of patients with metastasis")+
#   theme_classic()+
#   theme(axis.title = element_text(size=20),
#         axis.text  = element_text(size=20),
#         legend.position = "top",
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 16),
#         plot.title = element_text(size=22))+
#   ggtitle("Proportion of metastatic samples")+
#   geom_text(aes(x=mmmmTCGA,y=1,label = prop_group),vjust=-0.5,size=5)

fig3f <- survfit2(Surv(relapse_free_days, is_metastasis) ~ mmmmTCGA, data = df_prog_final %>% filter(relapse_free_days>0)) %>% 
  ggsurvfit() + 
  labs(x = "Days", y="Relapse-free survival probability")+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.position = "top",
        legend.text = element_text(size = 16))+
  ylim(0,1)

print(fig3g)

logrank_mmmmtcga <- survdiff(Surv(relapse_free_days, is_metastasis) ~ mmmmTCGA, data = df_prog_final %>% filter(relapse_free_days>0))
logrank_mmmmtcga$pvalue # 1.164318e-08



## mmmmmmTCGA
# crossing(tibble(mmmmmmTCGA = c("A","B","C","D")), tibble(is_metastasis = c(0,1))) %>% 
#   left_join(df_prog_final %>% count(mmmmmmTCGA,is_metastasis)) %>% 
#   mutate(n = ifelse(is.na(n), 0, n)) %>% 
#   mutate(is_metastasis = factor(is_metastasis,levels=c(1,0))) %>% 
#   group_by(mmmmmmTCGA) %>% mutate(n_group = sum(n)) %>% mutate(prop_group = n/n_group) %>% ungroup() %>% 
#   mutate(prop_group = ifelse(is_metastasis == 0, "", str_c(round(prop_group*100),"%"))) %>% 
#   ggplot(aes(mmmmmmTCGA,n,fill=is_metastasis))+
#   geom_bar(stat="identity", position=position_fill())+
#   scale_fill_manual(values = c("darkred","grey"))+
#   xlab("mmmmmmTCGA subtype")+ylab("Proportion of patients with metastasis")+
#   theme_classic()+
#   theme(axis.title = element_text(size=20),
#         axis.text  = element_text(size=20),
#         legend.position = "top",
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 16),
#         plot.title = element_text(size=22))+
#   ggtitle("Proportion of metastatic samples")+
#   geom_text(aes(x=mmmmmmTCGA,y=1,label = prop_group),vjust=-0.5,size=5)
# 
# survfit2(Surv(relapse_free_days, is_metastasis) ~ mmmmmmTCGA, data = df_prog_final %>% filter(relapse_free_days>0)) %>% 
#   ggsurvfit() + 
#   labs(x = "Days", y="Relapse-free survival probability")+
#   theme(axis.title = element_text(size=20),
#         axis.text = element_text(size=20),
#         legend.position = "top",
#         legend.text = element_text(size = 16))+
#   ylim(0,1)
# 
# logrank_mmmmmmTCGA <- survdiff(Surv(relapse_free_days, is_metastasis) ~ mmmmmmTCGA, data = df_prog_final %>% filter(relapse_free_days>0))
# logrank_mmmmmmTCGA$pvalue #  1.481453e-09
