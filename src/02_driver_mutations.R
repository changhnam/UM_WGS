##### Fig 2. Mutation-related genomic abererations in uveal melanomas
### Load packages & files
library(tidyverse)

df_st1 <- read_tsv("data/SuppTable1.tsv")

### Fig 2a. Frequency of alterations based on metastatic risk in our cohort
col_names <- df_st1 %>% select(contains("is_")) %>% select(-metastasis_organ) %>% colnames()

df_freq <- tibble()
for(cn in col_names){
  if(cn != "is_LOY"){
    df_freq <- bind_rows(df_freq,
                         tibble(alteration = cn,
                                n_ref = nrow(df_st1[df_st1[cn] == 0,]),
                                n_alt = nrow(df_st1[df_st1[cn] == 1,]),
                                n_ref_met = nrow(df_st1[df_st1[cn] == 0 & df_st1["is_metastasis"] == 1,]),
                                n_alt_met = nrow(df_st1[df_st1[cn] == 1 & df_st1["is_metastasis"] == 1,]),
                                pval = df_st1[c(cn,"is_metastasis")] %>% table() %>% chisq.test() %>% .$p.value))

  } else {
    df_freq <- bind_rows(df_freq,
                         tibble(alteration = cn,
                                n_ref = nrow(df_st1[df_st1[cn] == 0 & df_st1["sex"] == "M",]),
                                n_alt = nrow(df_st1[df_st1[cn] == 1 & df_st1["sex"] == "M",]),
                                n_ref_met = nrow(df_st1[df_st1[cn] == 0 & df_st1["sex"] == "M" & df_st1["is_metastasis"] == 1,]),
                                n_alt_met = nrow(df_st1[df_st1[cn] == 1 & df_st1["sex"] == "M" & df_st1["is_metastasis"] == 1,]),
                                pval = df_st1 %>% filter(sex == "M") %>% .[c(cn,"is_metastasis")] %>% table() %>% chisq.test() %>% .$p.value))
  }
}

df_freq <- df_freq %>%
  mutate(n_tot = n_ref + n_alt) %>%
  mutate(freq_alt = n_alt/n_tot*100) %>%
  mutate(freq_alt_in_met = n_alt_met/(n_ref_met+n_alt_met)*100) %>%
  mutate(freq_alt_in_nomet = (n_alt-n_alt_met)/(n_ref-n_ref_met+n_alt-n_alt_met)*100)

alt_order <- c("is_GNAQ","is_GNA11","is_SF3B1","is_EIF1AX","is_BAP1",
               "is_8q_gain","is_6p_gain","is_1q_gain",
               "is_6q_loss","is_1p_loss","is_3q_loss","is_3p_loss",
               "is_LOY","is_WGD","is_monosomy3","is_BAP1_LOH","is_BAP1_aberration")

fig2a1 <- df_freq %>% 
  filter(!alteration %in% c("is_metastasis","is_dead","is_kataegis","is_cgr")) %>% 
  mutate(alteration = factor(alteration, levels=rev(alt_order))) %>% 
  ggplot(aes(alteration,freq_alt_in_met))+
  geom_bar(fill = "#c03b3f", stat="identity")+
  coord_flip()+
  geom_text(aes(label = str_c(round(freq_alt_in_met,1)), hjust=+1.1), size=12)+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=24))+
  scale_x_discrete(position = "top")+
  scale_y_reverse()+
  ylim(100,0)

fig2a2 <- df_freq %>% 
  filter(!alteration %in% c("is_metastasis","is_dead","is_kataegis","is_cgr")) %>% 
  mutate(alteration = factor(alteration, levels=rev(alt_order))) %>% 
  ggplot(aes(alteration,freq_alt_in_nomet))+
  geom_bar(fill = "#c1e1bf", stat="identity")+
  coord_flip()+
  geom_text(aes(label = str_c(round(freq_alt_in_nomet,1)), hjust=-0.1), size=12)+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=24))+
  ylim(0,100)

print(fig2a1)
print(fig2a2)




### Fig 2c. Frequency of alterations between our cohort and ICGC cohort
fig2c <- df_freq %>% 
  filter(alteration %in% alt_order) %>% 
  mutate(alteration = factor(alteration, levels=rev(alt_order))) %>% 
  ggplot(aes(alteration,freq_alt))+
  geom_bar(fill = "#7bcfd7", stat="identity")+
  coord_flip()+
  geom_text(aes(label = str_c(round(freq_alt,1),"%"), hjust=-0.2), size=8)+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=16))+
  ylim(0,100)

print(fig2c)
