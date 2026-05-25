##### Fig1. Mutational landscape of uveal melanomas
### Load packages & files
library(tidyverse)
library(ComplexHeatmap)

df_st1 <- read_tsv("data/SuppTable1.tsv")

### Fig 1a. SNV burden & signature
sample_order_manual <- c("SEV-UM-014","SEV-UM-021","SEV-UM-036","SEV-UM-008","SEV-UM-001","SEV-UM-030",
                         "SEV-UM-005","SEV-UM-022",
                         "SEV-UM-040","SEV-UM-031","SEV-UM-020",
                         "SEV-UM-025","SEV-UM-006","SEV-UM-003","SEV-UM-009",
                         "SEV-UM-033","SEV-UM-015","SEV-UM-019","SEV-UM-024","SEV-UM-017","SEV-UM-011",
                         "SEV-UM-012",
                         "SEV-UM-037","SEV-UM-039","SEV-UM-027",
                         "SEV-UM-023","SEV-UM-032","SEV-UM-010","SEV-UM-038","SEV-UM-029",
                         "SEV-UM-028","SEV-UM-004","SEV-UM-013","SEV-UM-007",
                         "SEV-UM-035","SEV-UM-018","SEV-UM-002","SEV-UM-034","SEV-UM-016","SEV-UM-026")

df_sbs <- df_st1 %>% select(sampleID,n_snv,contains("prop_SBS")) %>% 
  mutate(`prop_SBS2/13` = prop_SBS2 + prop_SBS13) %>% select(-c(prop_SBS2, prop_SBS13)) %>% 
  pivot_longer(cols=contains("prop"), names_to = "Signature", values_to = "prop") %>% 
  mutate(Signature = str_remove(Signature,"prop_"), sig_count = n_snv*prop/100) %>% 
  mutate(Signature = factor(Signature, levels = c("SBS2/13","SBS18","SBS1","SBS5/40"))) %>% 
  arrange(Signature)

fig1a <- df_sbs %>% 
  mutate(sampleID = factor(sampleID, levels = sample_order_manual)) %>% 
  ggplot(aes(sampleID,sig_count,fill=Signature))+geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.title = element_text(size=16), 
        axis.text.x = element_text(size=16,angle=45,hjust=1),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))+
  xlab("")+ylab("Number of SNVs")+ylim(0,4000)+
  scale_fill_manual(values = c("#cc3333","#1e81b0","#36705c","#e6ab02"))

print(fig1a)





### Fig 1c. Driver & Clinical outcomes
# /home/users/changhyunnam/Projects/09_Uveal_melanoma/03_VCF/14_CNV/summary/hm_cnv_tcga_prognosis.pdf

# sample_order_hm <- c("SEV-UM-008","SEV-UM-001","SEV-UM-036","SEV-UM-021","SEV-UM-030","SEV-UM-033","SEV-UM-006","SEV-UM-040",
#                      "SEV-UM-019","SEV-UM-024","SEV-UM-025","SEV-UM-031","SEV-UM-014","SEV-UM-005","SEV-UM-022","SEV-UM-003",
#                      "SEV-UM-029","SEV-UM-010","SEV-UM-023","SEV-UM-028","SEV-UM-032","SEV-UM-037","SEV-UM-039","SEV-UM-027",
#                      "SEV-UM-038","SEV-UM-011","SEV-UM-015","SEV-UM-012","SEV-UM-017","SEV-UM-009","SEV-UM-004","SEV-UM-013",
#                      "SEV-UM-007","SEV-UM-016","SEV-UM-035","SEV-UM-034","SEV-UM-018","SEV-UM-026","SEV-UM-002","SEV-UM-020")

# sample_order_manual_new <- c("SEV-UM-014","SEV-UM-021","SEV-UM-036","SEV-UM-008","SEV-UM-001","SEV-UM-030",
#                              "SEV-UM-005","SEV-UM-022",
#                              "SEV-UM-040","SEV-UM-031","SEV-UM-020",
#                              "SEV-UM-025","SEV-UM-006","SEV-UM-003","SEV-UM-009",
#                              "SEV-UM-033","SEV-UM-015","SEV-UM-019","SEV-UM-024","SEV-UM-017","SEV-UM-011",
#                              "SEV-UM-012",
#                              "SEV-UM-037","SEV-UM-039","SEV-UM-027",
#                              "SEV-UM-023","SEV-UM-032","SEV-UM-010","SEV-UM-038","SEV-UM-029",
#                              "SEV-UM-028","SEV-UM-004","SEV-UM-013","SEV-UM-007",
#                              "SEV-UM-035","SEV-UM-018","SEV-UM-002","SEV-UM-034","SEV-UM-016","SEV-UM-026")

sample_order_hm <- sample_order_manual

# (a) driver
driver_list = c("GNAQ","GNA11","SF3B1","EIF1AX","BAP1")

df_driver <- df_st1 %>% 
  select(sampleID,contains(str_c("vartype_",driver_list))) %>% 
  rename_with(~ str_remove(.,"vartype_"), contains(c("vartype_"))) %>% 
  mutate(across(contains(driver_list), ~ replace_na(.,"None"))) %>% 
  mutate(across(contains(driver_list), ~ str_remove(.,";.+"))) %>% 
  mutate(sampleID = factor(sampleID, levels = sample_order_hm)) %>% arrange(sampleID)

mat_driver <- df_driver %>% as.matrix() %>% t()
colnames(mat_driver)=mat_driver[1,]
mat_driver <- mat_driver[-1,]

colors = structure(c("#0071bc","#c1272d","#fbb03b","#458b74","#7a378b","gray"), names = c("missense", "nonsense","frameshift","structural","splicing","None"))
hm_driver <- Heatmap(mat_driver,col=colors,rect_gp = gpar(col = "white", lwd = 2),border=TRUE, row_names_side="left")

# (b) CNV
df_cnv_table <- df_st1 %>% 
  select(sampleID,contains(c("gain","loss","LOY","WGD"))) %>% 
  rename_with(~str_remove(.,"is_"), contains(c("gain","loss","LOY","WGD"))) %>% 
  mutate(across(contains(c("gain","WGD")), ~ str_replace(.,"1","Gain"))) %>%
  mutate(across(contains(c("gain","WGD")), ~ str_replace(.,"0",""))) %>%
  mutate(across(contains(c("loss","LOY")), ~ str_replace(.,"1","Loss"))) %>%
  mutate(across(contains(c("loss","LOY")), ~ str_replace(.,"0","")))

mat_cnv <- df_cnv_table %>%
  pivot_longer(cols = c(`8q_gain`:WGD), names_to = "CN_event", values_to = "occur") %>%
  pivot_wider(names_from = sampleID, values_from = occur, values_fill = "") %>%
  select(CN_event,all_of(sample_order_hm)) %>%
  as.matrix()

rownames(mat_cnv) <- mat_cnv[,1]
mat_cnv <- mat_cnv[,-1]

colors = structure(c("#961E23","#162242","#CCCCCC"),
                   names = c("Gain","Loss","None"))

hm_cnv <- Heatmap(mat_cnv,col=colors,rect_gp = gpar(col = "white", lwd = 2),
                  cluster_rows = FALSE, show_column_dend = FALSE, row_dend_reorder = FALSE, border=TRUE,
                  row_names_side = "left") # show_heatmap_legend = FALSE,

# (c) chr3CN
mat_chr3CN <- df_st1 %>% select(sampleID,chr3CN) %>% 
  pivot_wider(names_from = sampleID, values_from = chr3CN) %>% 
  select(all_of(sample_order_hm)) %>% as.matrix()

rownames(mat_chr3CN) <- "chr3CN"

colors = c(colors, structure(c("#6492C1","#9AC4E9","white","#F2EE71","#FFCCCC"),
                             names = c("1","2, partial LOH", "2", "2-3, LOH", "4")))

hm_chr3CN <- Heatmap(mat_chr3CN,col=colors,rect_gp = gpar(col = "white", lwd = 2),
                     cluster_rows = FALSE, show_column_dend = FALSE, row_dend_reorder = FALSE, border=TRUE,
                     row_names_side = "left") # show_heatmap_legend = FALSE,

# (d) BAP1 LOH
mat_BAP1_LOH <- df_st1 %>% select(sampleID,is_BAP1_LOH) %>% 
  pivot_wider(names_from = sampleID, values_from = is_BAP1_LOH) %>% 
  select(all_of(sample_order_hm)) %>% as.matrix()

rownames(mat_BAP1_LOH) <- "BAP1_LOH"

colors = structure(c("white","red3"), names = c("0", "1"))

hm_BAP1_LOH <- Heatmap(mat_BAP1_LOH,col=colors,rect_gp = gpar(col = "white", lwd = 2),
                       cluster_rows = FALSE, show_column_dend = FALSE, row_dend_reorder = FALSE, border=TRUE,
                       row_names_side = "left") # show_heatmap_legend = FALSE,

# (e) TCGA
mat_TCGA <- df_st1 %>% select(sampleID,TCGA) %>%
  pivot_wider(names_from = sampleID, values_from = TCGA) %>%
  select(all_of(sample_order_hm)) %>% as.matrix()

rownames(mat_TCGA) <- "TCGA"

colors = structure(c("#C7DF95","#FFCC33","#FF6600","#993333"), names = c("A", "B", "C", "D"))

hm_TCGA <- Heatmap(mat_TCGA, col=colors,rect_gp = gpar(col = "white", lwd = 2),
                   cluster_rows = FALSE, show_column_dend = FALSE, row_dend_reorder = FALSE,border=TRUE,
                   show_heatmap_legend = TRUE,row_names_side = "left",
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", mat_TCGA[i, j]), x, y, gp = gpar(fontsize = 12))
                   })

# (f) prognosis
mat_prognosis <- df_st1 %>% select(sampleID,is_metastasis,is_dead) %>%
  arrange(factor(sampleID,levels=sample_order_hm)) %>%
  as.matrix() %>% t()

colnames(mat_prognosis) <- mat_prognosis[1,]
mat_prognosis <- mat_prognosis[-1,]

colors = structure(c("white","red3"), names = c("0", "1"))
hm_prognosis <- Heatmap(mat_prognosis, col=colors, rect_gp = gpar(col = "white", lwd = 2),
                        cluster_rows = FALSE, show_column_dend = FALSE, row_dend_reorder = FALSE,border=TRUE,
                        show_heatmap_legend = TRUE,row_names_side = "left")


fig1c <- hm_driver %v% hm_cnv %v% hm_chr3CN %v% hm_BAP1_LOH %v% hm_TCGA %v% hm_prognosis
print(fig1c)





### Fig. 1d. CNV co-occurence & exclusiveness

common_cnv_list <- df_st1 %>% select(contains(c("gain","loss","LOY"))) %>% colnames()

df_cnv_cor_common <- tibble()

for(cnv1 in common_cnv_list){
  cnv1_idx = which(common_cnv_list == cnv1)
  if(cnv1_idx == length(common_cnv_list)){next}

  for(cnv2 in common_cnv_list[(cnv1_idx+1):length(common_cnv_list)]){
    df_cnv12 <- df_st1[c("sampleID","sex",cnv1,cnv2)]
    df_cnv12 <- if(cnv1 == "is_LOY" | cnv2 == "LOY") {df_cnv12 %>% filter(sex == "M")} else {df_cnv12}
    df_cnv12["cnv1"] <- df_st1[cnv1]
    df_cnv12["cnv2"] <- df_st1[cnv2]
    df_cnv12 <- df_cnv12 %>% select(sampleID,cnv1,cnv2)
    chisq_cnv12 <- df_cnv12 %>% select(cnv1,cnv2) %>% table() %>% chisq.test(.)

    df_cnv_cor_common <- bind_rows(df_cnv_cor_common,
                                   tibble(cnv1 = cnv1, cnv2 = cnv2,
                                          case_00 = df_cnv12 %>% filter(cnv1 == 0, cnv2 == 0) %>% nrow(),
                                          case_01 = df_cnv12 %>% filter(cnv1 == 0, cnv2 == 1) %>% nrow(),
                                          case_10 = df_cnv12 %>% filter(cnv1 == 1, cnv2 == 0) %>% nrow(),
                                          case_11 = df_cnv12 %>% filter(cnv1 == 1, cnv2 == 1) %>% nrow(),
                                          phi = (case_11*case_00-case_01*case_10)/((case_00+case_01)*(case_10+case_11)*(case_00+case_10)*(case_01+case_11))**(1/2),
                                          pval = chisq_cnv12$p.value))
  }
}

df_cnv_cor_common <- df_cnv_cor_common %>% mutate(qval = p.adjust(pval,method="fdr")) %>% arrange(qval)

df_cnv_cor_common_v2 <- bind_rows(df_cnv_cor_common, df_cnv_cor_common %>% rename(cnv1 = cnv2, cnv2 = cnv1))

fig1d <- df_cnv_cor_common_v2 %>%
  mutate(cnv1 = factor(cnv1, levels=common_cnv_list),
         cnv2 = factor(cnv2, levels=common_cnv_list)) %>%
  bind_rows(., tibble(cnv1 = "dummy_neg", cnv2 = "dummy_neg", phi = -1, pval = 1)) %>%
  bind_rows(., tibble(cnv1 = "dummy_pos", cnv2 = "dummy_pos", phi = 1, pval = 1)) %>%
  mutate(signif = ifelse(qval < 0.05, "*", "")) %>%
  ggplot(aes(cnv1,cnv2))+
  geom_tile(aes(fill=phi),color="black")+
  geom_text(aes(label=signif),size=12,vjust=0.75)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text.x  = element_text(size=12,angle=90),
        axis.text.y  = element_text(size=12),
        axis.line = element_blank())

print(fig1d)
