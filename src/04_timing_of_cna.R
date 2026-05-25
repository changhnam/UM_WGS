##### Fig 4. Timing of prognostic copy gains in uveal melanomas
### Load packages & files
library(tidyverse)
library(ComplexHeatmap)

df_st1 <- read_tsv("data/SuppTable1.tsv")
df_timing <- read_tsv("data/timing_summary.tsv")

df_timing <- df_timing %>% 
  left_join(df_timing %>% filter(arm_name == "mrca") %>% rename(t_mrca = amp_time) %>% select(sampleID,t_mrca),by=c("sampleID")) %>% 
  left_join(df_timing %>% filter(arm_name == "diagnosis") %>% rename(age = amp_time) %>% select(sampleID,age),by=c("sampleID")) %>% 
  mutate(diff_time_mrca = t_mrca - amp_time, diff_time_age = age - amp_time,
         rela_time_mrca = amp_time/t_mrca, rela_time_age = amp_time/age)

# Fig 4c. Landscape of timing of CNAs
sample_order_age <- df_timing %>% filter(arm_name == "diagnosis") %>% arrange(amp_time) %>% pull(sampleID)
copy_gain_palette <- c("#6dbf51","#f5bf38","#be272d","#ff6600","#408dcb","#333334")
names(copy_gain_palette) <- c("6p","8q","1q","3p","mrca","diagnosis")

fig4c <- df_timing %>%
  select(sampleID,arm_name,amp_time) %>%
  filter(arm_name %in% c("6p","8q","1q","mrca","diagnosis")) %>% 
  mutate(arm_name = factor(arm_name, levels=c("6p","8q","1q","mrca","diagnosis"))) %>% 
  mutate(sampleID = factor(sampleID, levels = sample_order_age)) %>%
  ggplot(aes(sampleID,amp_time))+geom_point(aes(color=arm_name),size=3)+
  coord_flip()+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  xlab("")+ylab("Estimated time (years)")+
  scale_color_manual(values = copy_gain_palette)

print(fig4c)

df_hm_age <- df_st1 %>%
  select(sampleID,is_metastasis,is_BAP1,is_8q_gain,is_6p_gain,is_1q_gain,chr3CN,is_monosomy3,is_BAP1_LOH,is_BAP1_aberration) %>% 
  select(sampleID, is_metastasis, is_BAP1_aberration) %>% 
  filter(sampleID %in% sample_order_age) %>% 
  mutate(sampleID = factor(sampleID, levels = sample_order_age)) %>% 
  arrange(sampleID)

mat_age <- df_hm_age %>% as.matrix() %>% t()
colnames(mat_age) <- mat_age[1,]
mat_age <- mat_age[-1,]

colors = structure(c("red3", "white"), names = c("1", "0"))
fig4c_bottom <- Heatmap(mat_age, col=colors, rect_gp = gpar(col = "white", lwd = 2), border=TRUE, row_names_side = "left")

print(fig4c_bottom)





# Fig 4d. Timing of initial CNAs
df_timing_earliest <- df_timing %>% 
  group_by(sampleID) %>% mutate(rank = rank(amp_time)) %>% ungroup() %>% 
  filter(rank == 1 & arm_name %in% c("6p","8q","1q"))

fig4d <- df_timing_earliest %>% 
  mutate(is_metastasis = factor(is_metastasis)) %>% 
  ggplot(aes(amp_time,age))+
  geom_point(aes(fill = arm_name, color = is_metastasis), shape = 21, size=4, stroke = 1.2)+
  xlim(0,50)+ylim(20,85)+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  scale_color_manual(values = c("white","black"))+
  scale_fill_manual(values = copy_gain_palette)+
  scale_shape_manual(values = c(19,17))+
  geom_vline(xintercept = 20)

print(fig4d)





# Fig 4e. Initial CNAs by age groups
fig4e <- df_timing_earliest %>% 
  mutate(amp_time_group = ifelse(amp_time < 20, "young", "old")) %>% 
  count(amp_time_group,arm_name) %>% 
  mutate(arm_name = factor(arm_name, levels = c("1q","8q","6p"))) %>% 
  mutate(amp_time_group = factor(amp_time_group, levels = c("young","old"))) %>% 
  ggplot(aes(amp_time_group,n,fill=arm_name))+
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  scale_fill_manual(values = copy_gain_palette)+
  ylab("Proportion of initial CNA")
  
print(fig4e)





# Fig 4f. Total CNAs by age groups
fig4f <- df_timing_earliest %>% 
  mutate(amp_time_group = ifelse(amp_time < 20, "young", "old")) %>% 
  select(sampleID,amp_time_group) %>% 
  left_join(df_timing %>% select(sampleID,arm_name) %>% filter(arm_name %in% c("6p","8q","1q")), by=c("sampleID")) %>% 
  count(amp_time_group, arm_name) %>% 
  mutate(arm_name = factor(arm_name, levels = c("1q","8q","6p"))) %>% 
  mutate(amp_time_group = factor(amp_time_group, levels = c("young","old"))) %>% 
  ggplot(aes(amp_time_group,n, fill=arm_name))+
  geom_bar(stat="identity", position="fill")+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  scale_fill_manual(values = copy_gain_palette)+
  ylab("Proportion of total CNA")

print(fig4f)





# Fig 4g. Latency from initial CNAs to age at Dx
fig4g <- df_timing_earliest %>% 
  mutate(is_metastasis = factor(is_metastasis)) %>% 
  mutate(arm_name = factor(arm_name, levels = c("8q","6p","1q"))) %>% 
  ggplot(aes(x = arm_name, y = diff_time_age))+
  geom_point(aes(fill = arm_name, color = is_metastasis, shape = is_metastasis), shape = 21, size = 4, stroke = 1.2)+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  scale_color_manual(values = c("white","black"))+
  scale_fill_manual(values = copy_gain_palette)
  
print(fig4g)





# Fig 4h. Impact of initial CNAs on BAP1 aberration and metastasis
sample_order_earliest <- df_timing_earliest %>% 
  arrange(factor(arm_name, levels = c("6p","1q","8q")),
          amp_time) %>% 
  pull(sampleID)

fig4h <- df_timing %>% 
  filter(arm_name %in% c("6p","1q","8q","mrca","diagnosis")) %>% 
  filter(sampleID %in% sample_order_earliest) %>% 
  mutate(sampleID = factor(sampleID,levels=sample_order_earliest)) %>% 
  arrange(sampleID) %>% 
  ggplot(aes(sampleID,amp_time))+
  geom_point(aes(color=arm_name), size=5)+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text.x = element_text(size=16,angle=90),
        axis.text.y = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  xlab("")+ylab("Estimated time (years)")+
  scale_color_manual(values = copy_gain_palette)

print(fig4h)

df_hm_earliest <- df_st1 %>%
  select(sampleID,is_metastasis,is_BAP1,is_8q_gain,is_6p_gain,is_1q_gain,chr3CN,is_monosomy3,is_BAP1_LOH) %>% 
  mutate(is_BAP1_aberration = ifelse(is_BAP1 == 1 | is_BAP1_LOH == 1, 1, 0)) %>%
  select(sampleID, is_BAP1_aberration, is_metastasis) %>% 
  filter(sampleID %in% sample_order_earliest) %>% 
  mutate(sampleID = factor(sampleID, levels = sample_order_earliest)) %>% 
  arrange(sampleID)

mat_earliest <- df_hm_earliest %>% as.matrix() %>% t()
colnames(mat_earliest) <- mat_earliest[1,]
mat_earliest <- mat_earliest[-1,]

colors = structure(c("red3", "white"), names = c("1", "0"))
fig4h_bottom <- Heatmap(mat_earliest, col=colors, rect_gp = gpar(col = "white", lwd = 2), border=TRUE, row_names_side = "left")

print(fig4h_bottom)





### Fig 4i. Initial CNA ~ BAP1 aberration
fig4i <- df_timing_earliest %>%
  left_join(df_st1 %>% select(sampleID,is_BAP1_aberration),by=c("sampleID")) %>% 
  mutate(is_BAP1_aberration = factor(is_BAP1_aberration)) %>% 
  count(arm_name, is_BAP1_aberration) %>% 
  mutate(arm_name = factor(arm_name, levels=c("6p","1q","8q"))) %>% 
  ggplot(aes(arm_name, n, fill = is_BAP1_aberration))+
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size=22))+
  scale_fill_manual(values = c("grey","firebrick3"))+
  xlab("earliest_6p_gain")+
  ylab("Proportion of patients with BAP1 pathology")

print(fig4i)





### Fig 4j. Initial CNA ~ metastasis
fig4j <- df_timing_earliest %>% 
  mutate(arm_name_new = ifelse(arm_name == "6p", "6p", "non-6p")) %>% 
  mutate(is_metastasis = factor(is_metastasis)) %>% 
  count(arm_name, is_metastasis) %>% 
  mutate(arm_name = factor(arm_name, levels=c("6p","1q","8q"))) %>% 
  ggplot(aes(arm_name, n, fill = is_metastasis))+
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size=22))+
  scale_fill_manual(values = c("grey","darkred"))+
  xlab("earliest_6p_gain")+
  ylab("Proportion of patients with metastasis")

print(fig4j)





### Fig 4k. 6p gain (initial vs. non-initial) ~ BAP1 aberration
fig4k <- df_timing %>% filter(arm_name == "6p") %>% 
  select(sampleID) %>% 
  left_join(df_timing_earliest %>% select(sampleID, arm_name),by=c("sampleID")) %>% 
  left_join(df_st1 %>% select(sampleID,is_BAP1_aberration), by=c("sampleID")) %>% 
  mutate(arm_name_new = ifelse(arm_name == "6p", "initial_6p", "non_initial_6p")) %>% 
  count(arm_name_new, is_BAP1_aberration) %>% 
  mutate(is_BAP1_aberration = factor(is_BAP1_aberration)) %>% 
  ggplot(aes(arm_name_new, n, fill = is_BAP1_aberration))+
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size=22))+
  scale_fill_manual(values = c("grey","firebrick3"))+
  xlab("timing_of_6p_gain")+
  ylab("Proportion of patients with BAP1 pathology")

print(fig4k)





### Fig 4l. 6p gain (initial vs. non-initial) ~ metastasis
fig4l <- df_timing %>% filter(arm_name == "6p") %>% 
  select(sampleID, is_metastasis) %>% 
  left_join(df_timing_earliest %>% select(sampleID, arm_name),by=c("sampleID")) %>% 
  mutate(arm_name_new = ifelse(arm_name == "6p", "initial_6p", "non_initial_6p")) %>% 
  count(arm_name_new, is_metastasis) %>% 
  mutate(is_metastasis = factor(is_metastasis)) %>% 
  ggplot(aes(arm_name_new, n, fill = is_metastasis))+
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=20),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size=22))+
  scale_fill_manual(values = c("grey","darkred"))+
  xlab("timing_of_6p_gain")+
  ylab("Proportion of patients with metastasis")

print(fig4l)
