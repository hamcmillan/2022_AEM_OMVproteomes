library(tidyverse)
library(dplyr)
library(readxl)
library(missForest)
library(EnhancedVolcano)


# Read abundance data, get only log-transformed columns (will normalize differently to compare),
# filter for only those with significant interaction effects
abundance <- read_tsv(file="processed_data/total_processed_abundance_2wayANOVA.txt") %>% 
  filter(significant_interaction == TRUE) %>% 
  select(1:62, "significant_interaction")

#how many missing values per row of each condition
pf_complete_na <- rowSums(is.na(dplyr::select(abundance, gi_number, description, pf_complete_1, pf_complete_2, pf_complete_3)))
pf_minimal_na <- rowSums(is.na(dplyr::select(abundance, gi_number, description, pf_minimal_1, pf_minimal_2, pf_minimal_3)))
pst_complete_na <- rowSums(is.na(dplyr::select(abundance, gi_number, description, pst_complete_1, pst_complete_2, pst_complete_3)))
pst_minimal_na <- rowSums(is.na(dplyr::select(abundance, gi_number, description, pst_minimal_1, pst_minimal_2, pst_minimal_3)))

#append these lists to main table
#if 2 or more missing values, add to new data table "differential"
abundance <- mutate(abundance, pf_comp_na = pf_complete_na, pf_min_na = pf_minimal_na,
                    pst_comp_na = pst_complete_na, pst_min_na = pst_minimal_na)
differential <- data.frame(filter(abundance, pf_comp_na >= 2 | pf_min_na >= 2 | pst_comp_na >= 2 | pst_min_na >= 2))

#save/write differentially expressed to new directory
write_tsv(differential, "processed_data/total_differential_abundance_OnlySigInteraction.tsv")

# ##  normalize abundance data
# # test for normality
# for(i in 15:23){
#   print(shapiro.test(pull(abundance, i)))
# }
# 
# # log2 transform, add transformed to df, "test" contains results of transformed shapiro
# newname <- character()
# test <- numeric()
# for(i in 15:23){
#   newname <- c(newname, paste("log_", colnames(abundance)[[i]], sep = "_"))
#   abundance <- mutate(abundance,
#                       "{newname[[i-14]]}" := log2(abundance[[i]]))
#   test <- c(test, print(shapiro.test(pull(abundance, newname[[i-14]]))))
# }
# 
# # pst_abundance_go <- mutate(pst_abundance_go, log_pst_complete_1 = log2(pst_complete_1))
# # shapiro.test(pull(pst_abundance_go, "log_pst_complete_1"))
# 
# # pst_abundance_go <- mutate(pst_abundance_go, avg_pst_complete_1 = pst_complete_1/summary(pst_abundance_go$pst_complete_1)[["Mean"]])
# 
# # plot overlayed normalized distributions 
# log_for_plot <- pst_abundance_go %>% 
#   dplyr::select(31:39) %>% 
#   pivot_longer(cols = 1:9,names_to = "log_sample", values_to = "values")
# 
# log_for_plot %>% 
#   ggplot(aes(x=values, color=log_sample)) +
#   geom_density()

# mean-normalize Pf complete to Pst complete and Pf minimal to Pst minimal, leave out the pooled samples
abundance <- mutate(abundance, avg_log_pst_complete = 
                      abundance %>% 
                      dplyr::select(49:51) %>% 
                      rowMeans())

abundance <- mutate(abundance, avg_log_pst_minimal = 
                      abundance %>% 
                      dplyr::select(52:54) %>% 
                      rowMeans())

abundance <- mutate(abundance, mc_log_pf_complete_1 = abundance$log__pf_complete_1 / summary(abundance$avg_log_pst_complete)[["Mean"]])
abundance <- mutate(abundance, mc_log_pf_complete_2 = abundance$log__pf_complete_2 / summary(abundance$avg_log_pst_complete)[["Mean"]])
abundance <- mutate(abundance, mc_log_pf_complete_3 = abundance$log__pf_complete_3 / summary(abundance$avg_log_pst_complete)[["Mean"]])

abundance <- mutate(abundance, mc_log_pf_minimal_1 = abundance$log__pf_minimal_1 / summary(abundance$avg_log_pst_minimal)[["Mean"]])
abundance <- mutate(abundance, mc_log_pf_minimal_2 = abundance$log__pf_minimal_2 / summary(abundance$avg_log_pst_minimal)[["Mean"]])
abundance <- mutate(abundance, mc_log_pf_minimal_3 = abundance$log__pf_minimal_3 / summary(abundance$avg_log_pst_minimal)[["Mean"]])

abundance <- mutate(abundance, mc_log_pst_complete_1 = abundance$log__pst_complete_1 / summary(abundance$avg_log_pst_complete)[["Mean"]])
abundance <- mutate(abundance, mc_log_pst_complete_2 = abundance$log__pst_complete_2 / summary(abundance$avg_log_pst_complete)[["Mean"]])
abundance <- mutate(abundance, mc_log_pst_complete_3 = abundance$log__pst_complete_3 / summary(abundance$avg_log_pst_complete)[["Mean"]])

abundance <- mutate(abundance, mc_log_pst_minimal_1 = abundance$log__pst_minimal_1 / summary(abundance$avg_log_pst_minimal)[["Mean"]])
abundance <- mutate(abundance, mc_log_pst_minimal_2 = abundance$log__pst_minimal_2 / summary(abundance$avg_log_pst_minimal)[["Mean"]])
abundance <- mutate(abundance, mc_log_pst_minimal_3 = abundance$log__pst_minimal_3 / summary(abundance$avg_log_pst_minimal)[["Mean"]])

mc_log_for_plot <- abundance %>% 
  dplyr::select(65:76) %>% 
  pivot_longer(cols = 1:12, names_to = "mc_log_sample", values_to = "values")

mc_log_for_plot %>% 
  ggplot(aes(x=values, color=mc_log_sample))+
  geom_density()

#impute missing values
missing <- abundance %>% 
  dplyr::select(65:76) %>% 
  pivot_longer(cols = 1:12, names_to = "mc_log_sample", values_to = "values") %>% 
  mutate_all(~replace_na(.,0))

missing %>% 
  ggplot(aes(x=values, color=mc_log_sample))+
  geom_density()

hm_missing <- abundance %>% 
  dplyr::select(65:76) %>% 
  mutate_all(~replace_na(.,0))
color <- ifelse(hm_missing == 0, "white", "black")
heatmap(as.matrix(hm_missing), col = color, labRow = abundance$description)

# hm_missing <- pst_abundance_go %>% 
#   select(41:46) %>%
#   mutate_all(~replace(.,!is.na(.),1)) %>% 
#   mutate_all(~replace_na(.,0))
# heatmap(as.matrix(hm_missing))  


blanks_abund <- abundance %>% 
  dplyr::select(65:76) %>% 
  as.matrix()

imp_abund <- missForest(blanks_abund)

imp_for_plot <- imp_abund$ximp %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "imp_mc_log_sample", values_to = "values")

imp_for_plot %>% 
  ggplot(aes(x=values, color=imp_mc_log_sample))+
  geom_density()

# add imputed data to dataframe, test for differential expression
imp_for_tests <- as.data.frame(imp_abund$ximp)

abundance <- mutate(abundance, imp_mc_log_pf_complete_1 = imp_for_tests$mc_log_pf_complete_1)
abundance <- mutate(abundance, imp_mc_log_pf_complete_2 = imp_for_tests$mc_log_pf_complete_2)
abundance <- mutate(abundance, imp_mc_log_pf_complete_3 = imp_for_tests$mc_log_pf_complete_3)

abundance <- mutate(abundance, imp_mc_log_pf_minimal_1 = imp_for_tests$mc_log_pf_minimal_1)
abundance <- mutate(abundance, imp_mc_log_pf_minimal_2 = imp_for_tests$mc_log_pf_minimal_2)
abundance <- mutate(abundance, imp_mc_log_pf_minimal_3 = imp_for_tests$mc_log_pf_minimal_3)

abundance <- mutate(abundance, imp_mc_log_pst_complete_1 = imp_for_tests$mc_log_pst_complete_1)
abundance <- mutate(abundance, imp_mc_log_pst_complete_2 = imp_for_tests$mc_log_pst_complete_2)
abundance <- mutate(abundance, imp_mc_log_pst_complete_3 = imp_for_tests$mc_log_pst_complete_3)

abundance <- mutate(abundance, imp_mc_log_pst_minimal_1 = imp_for_tests$mc_log_pst_minimal_1)
abundance <- mutate(abundance, imp_mc_log_pst_minimal_2 = imp_for_tests$mc_log_pst_minimal_2)
abundance <- mutate(abundance, imp_mc_log_pst_minimal_3 = imp_for_tests$mc_log_pst_minimal_3)

# test for differential expression
imp_for_t <- as.matrix(imp_for_tests)

complete_p <- double()
minimal_p <- double()
for(i in 1:nrow(imp_for_t)){
  pf_complete = c(imp_for_t[i, 1:3])
  pf_minimal = c(imp_for_t[i, 4:6])
  pst_complete = c(imp_for_t[i, 7:9])
  pst_minimal = c(imp_for_t[i, 10:12])
  complete_p <- c(complete_p, print(t.test(pf_complete, pst_complete)$p.value))
  minimal_p <- c(minimal_p, print(t.test(pf_minimal, pst_minimal)$p.value))
}

# add p-values to dataframe, calculate and add adjusted p-values
# add logical column, TRUE = adj significant, FALSE = adj not significant
abundance <- mutate(abundance, complete_p = complete_p, minimal_p = minimal_p)
abundance <- abundance %>% 
  mutate(adj_complete_p = p.adjust(.$complete_p, method = "BH")) %>% 
  mutate(adj_minimal_p = p.adjust(.$minimal_p, method = "BH"))

abundance <- abundance %>% 
  mutate(significant_comp = if_else(adj_complete_p <= 0.05, "TRUE", "FALSE")) %>%
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(adj_minimal_p <= 0.05, "TRUE", "FALSE")) %>%
  mutate(significant_min = as.logical(significant_min))

# volcano plot
# calculate fold change in each media and plot 
avg_imp_mc_log_pf_complete <- double()
avg_imp_mc_log_pf_minimal <- double()
avg_imp_mc_log_pst_complete <- double()
avg_imp_mc_log_pst_minimal <- double()
for(i in 1:nrow(abundance)){
  pf_complete = c(abundance[i, 77:79])
  pf_minimal = c(abundance[i, 80:82])
  pst_complete = c(abundance[i, 83:85])
  pst_minimal = c(abundance[i, 86:88])
  avg_imp_mc_log_pf_complete <- c(avg_imp_mc_log_pf_complete, print(mean(as.double(pf_complete))))
  avg_imp_mc_log_pf_minimal <- c(avg_imp_mc_log_pf_minimal, print(mean(as.double(pf_minimal))))
  avg_imp_mc_log_pst_complete <- c(avg_imp_mc_log_pst_complete, print(mean(as.double(pst_complete))))
  avg_imp_mc_log_pst_minimal <- c(avg_imp_mc_log_pst_minimal, print(mean(as.double(pst_minimal))))
}

abundance <- mutate(abundance, avg_imp_mc_log_pf_complete = avg_imp_mc_log_pf_complete)
abundance <- mutate(abundance, avg_imp_mc_log_pf_minimal = avg_imp_mc_log_pf_minimal)
abundance <- mutate(abundance, avg_imp_mc_log_pst_complete = avg_imp_mc_log_pst_complete)
abundance <- mutate(abundance, avg_imp_mc_log_pst_minimal = avg_imp_mc_log_pst_minimal)

vp <- abundance %>% 
  dplyr::select("accession", contains("imp_mc_log"), "adj_complete_p", "adj_minimal_p", 
                "significant_comp", "significant_min") %>% 
  mutate(comp_fc = avg_imp_mc_log_pf_complete/avg_imp_mc_log_pst_complete) %>% 
  mutate(min_fc = avg_imp_mc_log_pf_minimal/avg_imp_mc_log_pst_minimal) %>% 
  mutate(log_adj_comp_p = -log10(adj_complete_p)) %>% 
  mutate(log_adj_min_p = -log10(adj_minimal_p)) %>% 
  mutate(log2_comp_fc = log2(comp_fc)) %>% 
  mutate(log2_min_fc = log2(min_fc))

abundance <- mutate(abundance, comp_fc = vp$comp_fc)
abundance <- mutate(abundance, min_fc = vp$min_fc)
abundance <- mutate(abundance, nlog10_adj_comp_p = vp$log_adj_comp_p)
abundance <- mutate(abundance, nlog10_adj_min_p = vp$log_adj_min_p)
abundance <- mutate(abundance, log2_comp_fc = vp$log2_comp_fc)
abundance <- mutate(abundance, log2_min_fc = vp$log2_min_fc)

# vp %>% 
#   ggplot(aes(x = log2_fc, y = log_p_adj, color = significant_comp_vs_min))+
#   geom_point()

# comp volcano with t test values
EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_comp_fc", 
                y = "adj_complete_p",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 3),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
                )

# min volcano with t test values
EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_min_fc", 
                y = "adj_minimal_p",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 3),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
                )

# write new file with complete abundance, go terms, and significance testing
write_tsv(abundance, "processed_data/total_processed_abundance_t_test.txt")

# get list of significant difference by t test 
t_aov_sigdif_comp <- abundance %>% 
  dplyr::select("gi_number", "accession", "description", contains("imp_mc_log"), "adj_complete_p",
                "adj_minimal_p", "significant_comp", "significant_min", contains("log2")) %>% 
  filter(significant_comp == TRUE)

t_aov_sigenrich_pf_comp <- t_aov_sigdif_comp %>% 
  filter(avg_imp_mc_log_pf_complete > avg_imp_mc_log_pst_complete)
t_aov_sigenrich_pst_comp <- t_aov_sigdif_comp %>% 
  filter(avg_imp_mc_log_pst_complete > avg_imp_mc_log_pf_complete)
t_aov_enrich_and_fc_pf_comp <- t_aov_sigenrich_pf_comp %>% 
  filter(log2_comp_fc >= 0.5)
t_aov_enrich_and_fc_pst_comp <- t_aov_sigenrich_pst_comp %>% 
  filter(abs(log2_comp_fc) >= 0.5)


#save/write differentially expressed to new directory
write_tsv(t_aov_sigenrich_pf_comp, "processed_data/pf_sigenrich_comp_t_and_2wayANOVA.tsv")
write_tsv(t_aov_sigenrich_pst_comp, "processed_data/pst_sigenrich_comp_t_and_2wayANOVA.tsv")
write_tsv(t_aov_enrich_and_fc_pf_comp, "processed_data/pf_enrich_and_fc_comp_t_and_2wayANOVA.tsv")
write_tsv(t_aov_enrich_and_fc_pst_comp, "processed_data/pst_enrich_and_fc_comp_t_and_2wayANOVA.tsv")


################ compare only those with increased packaging in minimal in Pst #################
# here we ask of those that by 2-way anova have a significant interaction effect
# and have increased packaging in Pst minimal compared to Pst complete,
# which of these have the same baseline abundance in Pf complete and Pst complete
# (ie. which may be induced in minimal for virulence functions)
# and which of these have different baselines in Pf complete and Pst complete
# (ie. are perhaps induced in minimal for nutrient acquisition/starvation response)

# Read abundance data from the 2-way anova test,
# filter for only those with significant interaction effects
# and those that have pst_fc >1 (display increased packaging in minimal media)
# and are significantly different by 2-way anova plus tukey in Pst minimal vs complete
# make also a list with a fold change higher than log2 0.5
# keep only accession identifier to merge, pst_fc from anova-normalized set
total_abundance_go <- read_tsv(file="processed_data/total_processed_abundance_2wayANOVA.txt")

total_sig_int_and_high_min <- total_abundance_go %>%  
  filter(significant_interaction == TRUE & 
           adj_tukey_p_pst_comp_v_min <= 0.05 &
           pst_fc > 1)

total_sig_int_and_high_min_with_fc <- total_sig_int_and_high_min %>% 
  filter(log2_pst_fc >= 0.5)

total_merge <- total_sig_int_and_high_min %>% 
  select("accession", "pf_fc", "pst_fc", "log2_pf_fc", "log2_pst_fc", "adj_tukey_p_pst_comp_v_min", 
         "group_pst_min", "group_pst_comp", "group_pf_comp", "group_pf_min")

total_fc_merge <- total_sig_int_and_high_min_with_fc %>% 
  select("accession", "pf_fc", "pst_fc", "log2_pf_fc", "log2_pst_fc", "adj_tukey_p_pst_comp_v_min", 
         "group_pst_min", "group_pst_comp", "group_pf_comp", "group_pf_min")


# combine the sets for the merge with the abundance data
# we want to know, of the set that has a significant interaction effect media x species
# and is significantly higher in Pst min vs Pst comp,
# which of these have the same baseline in Pf comp and Pst comp

higher_min <- inner_join(abundance, total_merge, by = "accession")
fc_higher_min <- inner_join(abundance, total_fc_merge, by = "accession")

# comp volcano with t test values filter for higher in Pst min
EnhancedVolcano(higher_min, 
                lab = higher_min$accession, 
                x = "log2_comp_fc", 
                y = "adj_complete_p",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 3),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)

# comp volcano with t test values filter for fc higher in Pst min
#svglite::svglite("fc_enrich_Pst_min_Pfcomp_vs_Pstcomp.svg")
EnhancedVolcano(fc_higher_min, 
                lab = fc_higher_min$accession, 
                x = "log2_comp_fc", 
                y = "adj_complete_p",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 3),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)
#dev.off()

# write new file with the subsets of abundance data
write_tsv(higher_min, "processed_data/enrich_pst_min_processed_abundance_t_test_for_baseline.txt")
write_tsv(fc_higher_min, "processed_data/fc_enrich_pst_min_processed_abundance_t_test_for_baseline.txt")

# get lists from higher_min subset:
# 1) sig interaction effect, sig higher in Pst min, no difference between Pst comp and Pf comp (same baseline, maybe involved in virulence)
# 2) sig interaction effect, sig higher in Pst min, higher in Pf comp than Pst comp (higher baseline in Pf, maybe contribute to Pf robustness)
# 3) sig interaction effect, sig higher in Pst min, lower in Pf comp than Pst comp (lower baseline in Pf, maybe involved in pf-specific response?)
same_baseline <- higher_min %>% 
  filter(significant_comp == FALSE)
higher_baseline <- higher_min %>% 
  filter(significant_comp == TRUE & comp_fc > 1)
lower_baseline <- higher_min %>% 
  filter(significant_comp == TRUE & comp_fc < 1)

fc_higher_baseline <- higher_baseline %>% 
  filter(significant_comp == TRUE & log2_comp_fc >= 0.5)
fc_lower_baseline <- lower_baseline %>% 
  filter(significant_comp == TRUE & abs(log2_comp_fc) >= 0.5)

#save/write to new directory
write_tsv(same_baseline, "processed_data/enrich_Pst_same_Pf_baseline.tsv")
write_tsv(higher_baseline, "processed_data/enrich_Pst_higher_Pf_baseline.tsv")
write_tsv(lower_baseline, "processed_data/enrich_Pst_lower_Pf_baseline.tsv")
write_tsv(fc_higher_baseline, "processed_data/enrich_Pst_fc_higher_Pf_baseline.tsv")
write_tsv(fc_lower_baseline, "processed_data/enrich_Pst_fc_lower_Pf_baseline.tsv")

# repeat get lists from fc_higher_min subset:
# 1) sig interaction effect, sig higher in Pst min, no difference between Pst comp and Pf comp (same baseline, maybe involved in virulence)
# 2) sig interaction effect, sig higher in Pst min, higher in Pf comp than Pst comp (higher baseline in Pf, maybe contribute to Pf robustness)
# 3) sig interaction effect, sig higher in Pst min, lower in Pf comp than Pst comp (lower baseline in Pf, maybe involved in pf-specific response?)
fcpst_same_baseline <- fc_higher_min %>% 
  filter(significant_comp == FALSE)
fcpst_higher_baseline <- fc_higher_min %>% 
  filter(significant_comp == TRUE & comp_fc > 1)
fcpst_lower_baseline <- fc_higher_min %>% 
  filter(significant_comp == TRUE & comp_fc < 1)

fcpst_fc_higher_baseline <- fcpst_higher_baseline %>% 
  filter(significant_comp == TRUE & log2_comp_fc >= 0.5)
fcpst_fc_lower_baseline <- fcpst_lower_baseline %>% 
  filter(significant_comp == TRUE & abs(log2_comp_fc) >= 0.5)

#save/write to new directory
write_tsv(fcpst_same_baseline, "processed_data/fc_enrich_Pst_same_Pf_baseline.tsv")
write_tsv(fcpst_higher_baseline, "processed_data/fc_enrich_Pst_higher_Pf_baseline.tsv")
write_tsv(fcpst_lower_baseline, "processed_data/fc_enrich_Pst_lower_Pf_baseline.tsv")
write_tsv(fcpst_fc_higher_baseline, "processed_data/fc_enrich_Pst_fc_higher_Pf_baseline.tsv")
write_tsv(fcpst_fc_lower_baseline, "processed_data/fc_enrich_Pst_fc_lower_Pf_baseline.tsv")


# # lists of interesting categories within the significantly enriched in minimal set
# sig_pathogenesis <- pst_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("pathogenesis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_photosynthesis <- pst_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("photosynthesis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_iron <- pst_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("iron", .$description)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_symbiosis <- pst_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("symbiosis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)

# heatmap to explore data
hm_total <- abundance %>% 
  dplyr::select(40:45, 49:54)

rownames(hm_total) <- abundance$accession

drop_na(hm_total)

# heatmap(as.matrix(hm_total), labRow = pst_abundance_go$description)

############### make pretty heatmap ##################
library(pheatmap)
pheatmap(as.matrix(hm_total), na_col = "white")

hm_total_scale <- scale(as.matrix(hm_total))
pheatmap(na.omit(hm_total_scale), na_col = "white", cutree_rows = 10, fontsize_row = 0.5, show_rownames = FALSE)

# use silent = TRUE to suppress the plot
hm_clustering <- pheatmap(na.omit(hm_total_scale), cutree_rows = 10, silent = TRUE, fontsize_row = 0.5)

# get the dendrogram
hm_clustering$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)

# get list of clusters
hm_clustering_list <- cbind(hm_total_scale, cluster = cutree(hm_clustering$tree_row, k = 10))

# heatmap with replicates grouped (using imputed values)
hm_grouped <- select(abundance, "avg_imp_mc_log_pf_complete", "avg_imp_mc_log_pf_minimal",
                     "avg_imp_mc_log_pst_complete", "avg_imp_mc_log_pst_minimal")
rownames(hm_grouped) <- abundance$accession
hm_grouped <- scale(as.matrix(hm_grouped))
pheatmap(hm_grouped, show_rownames = FALSE, cluster_cols = TRUE, cutree_rows = 10, cellwidth = 75)

hm_grouped_clust <- pheatmap(hm_grouped, show_rownames = FALSE, cluster_cols = TRUE, cutree_rows = 10, cellwidth = 75)
hm_grouped_clust_list <- cbind(hm_grouped, cluster = cutree(hm_grouped_clust$tree_row, k = 10))


# PCA to explore data
# library(limma)

total_pca <- abundance %>%
  dplyr::select(40:60) %>% # log transformed data
  drop_na() %>%
  prcomp(center = TRUE, scale. = TRUE) #scale. = TRUE


color <- case_when(grepl("pf_complete", rownames(total_pca$rotation)) ~ "pf_complete", 
                   grepl("pf_minimal", rownames(total_pca$rotation)) ~ "pf_minimal",
                   grepl("pf_pool", rownames(total_pca$rotation)) ~ "pf_pool",
                   grepl("pst_complete", rownames(total_pca$rotation)) ~ "pst_complete", 
                   grepl("pst_minimal", rownames(total_pca$rotation)) ~ "pst_minimal",
                   grepl("pst_pool", rownames(total_pca$rotation)) ~ "pst_pool",
                   grepl("total_pool", rownames(total_pca$rotation)) ~ "total_pool")


x_lab = paste("PC1", " (", (summary(total_pca)$importance[2,1])*100, "%)", sep = "")
y_lab = paste("PC2", " (", (summary(total_pca)$importance[2,2])*100, "%)", sep = "")
ggplot(as.data.frame(total_pca$rotation), aes(x=PC1*100, y=PC2*100, color=color)) +
  geom_point(shape=19, size=2) +
  #coord_fixed() +
  labs(title=NULL,
       x=x_lab,
       y=y_lab) +
  theme_classic()


############### export graphs ##################
# setEPS()
# postscript("filename", width = 4, height = 2)
# ggplot(or whatever you use to make the plot)
# dev.off()

############### KEGG pathway analysis #####################
# Resources for contingency tables and fisher test:
# https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
# https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
# library(BiocManager)
# library(KEGGREST)
# library(KEGGprofile)
# library(biomaRt)


# using only identified proteins as the background set (lines 296-475)
total_k_terms <- read_excel("raw_data/kegg_Pst_total.xlsx") %>% 
  dplyr::select("gi_number", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

kegg <- pst_abundance_go %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
         contains("p_value"), contains("significant"), "fc_min_vs_comp", "log2_fc_min_vs_comp") %>% 
  inner_join(total_k_terms, by = "gi_number")

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make contingency_tables

diff_abund <- sum(grepl(TRUE, kegg$significant_comp_vs_min))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(k_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, k_term_count[i,2])
  number_not_ko <- c(number_not_ko, sum(grepl(TRUE, !is.na(kegg$ko))) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_comp_vs_min) & grepl(k_term_count[i,1], kegg$ko)))
  contingency_table[[i]] <- matrix(c(
    as.numeric(number_diff_and_ko[i]),
    diff_abund - as.numeric(number_diff_and_ko[i]),
    as.numeric(number_in_ko[i]) - as.numeric(number_diff_and_ko[i]),
    as.numeric(number_not_ko[i] - (diff_abund - as.numeric(number_diff_and_ko[i])))),
    ncol = 2)
}

# run Fisher's exact test

fisher_p <- double()
for(i in 1:length(contingency_table)){
  fisher_p <- c(fisher_p, fisher.test(contingency_table[[i]])$p.value)
}

k_term_fisher <- k_term_count %>% 
  filter(ko != is.na(ko)) %>% 
  mutate(fisher_p_value = fisher_p) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))


# use p. putida as comparison 
# https://academic.oup.com/mbe/article/28/1/483/982857#77836793
# https://www.mdpi.com/2073-4425/11/2/139
# https://link.springer.com/article/10.1007/s10658-017-1401-8

# convert GI to NCBI ID using:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=GIs separated by commas&rettype=acc

#just kidding...can upload FASTA file to blast KOALA https://www.kegg.jp/blastkoala/ and get list of k terms


k_term_nona <- total_k_terms %>% 
  mutate(mapped_to_k = grepl("", total_k_terms$ko)) %>% 
  filter(.$mapped_to_k == TRUE) %>% 
  dplyr::select(ko)

kegg_db <- read_tsv("raw_data/kegg_hierarchy.txt")

k_term_for_db <- character()
for(i in 1:nrow(kegg_db)){
  k_term_for_db <- c(k_term_for_db, gsub(pattern = "(D      )(.*)(  .*)", replacement = "\\2", x = kegg_db[[i,4]]))
}

kegg_db <- mutate(kegg_db, ko = k_term_for_db)

kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  full_join(kegg, kegg_db, by = "ko")


## evaluate highest kegg level A
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

diff_abund <- sum(grepl(TRUE, kegg$significant_comp_vs_min))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(k_path_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, k_path_count[i,2])
  number_not_path <- c(number_not_path, sum(grepl(TRUE, !is.na(kegg_with_path$A))) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_comp_vs_min) & grepl(k_path_count[i,1], kegg_with_path$ko)))
  contingency_table_path[[i]] <- matrix(c(
    as.numeric(number_diff_and_path[i]),
    diff_abund - as.numeric(number_diff_and_path[i]),
    as.numeric(number_in_path[i]) - as.numeric(number_diff_and_path[i]),
    as.numeric(number_not_path[i] - (diff_abund - as.numeric(number_diff_and_path[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_path <- double()
for(i in 1:length(contingency_table_path)){
  fisher_p_path <- c(fisher_p_path, fisher.test(contingency_table_path[[i]])$p.value)
}

k_path_fisher <- k_path_count %>% 
  filter(A != is.na(A)) %>% 
  mutate(fisher_p_value = fisher_p_path) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

## evaluate kegg level B
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

diff_abund <- sum(grepl(TRUE, kegg$significant_comp_vs_min))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(k_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, k_pathB_count[i,2])
  number_not_pathB <- c(number_not_pathB, sum(grepl(TRUE, !is.na(kegg_with_path$B))) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_comp_vs_min) & grepl(k_pathB_count[i,1], kegg_with_path$ko)))
  contingency_table_pathB[[i]] <- matrix(c(
    as.numeric(number_diff_and_pathB[i]),
    diff_abund - as.numeric(number_diff_and_pathB[i]),
    as.numeric(number_in_pathB[i]) - as.numeric(number_diff_and_pathB[i]),
    as.numeric(number_not_pathB[i] - (diff_abund - as.numeric(number_diff_and_pathB[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_pathB <- double()
for(i in 1:length(contingency_table_pathB)){
  fisher_p_pathB <- c(fisher_p_pathB, fisher.test(contingency_table_pathB[[i]])$p.value)
}

k_pathB_fisher <- k_pathB_count %>% 
  filter(B != is.na(B)) %>% 
  mutate(fisher_p_value = fisher_p_pathB) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

## evaluate kegg level C
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

diff_abund <- sum(grepl(TRUE, kegg$significant_comp_vs_min))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(k_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, k_pathC_count[i,2])
  number_not_pathC <- c(number_not_pathC, sum(grepl(TRUE, !is.na(kegg_with_path$C))) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_comp_vs_min) & grepl(k_pathC_count[i,1], kegg_with_path$ko)))
  contingency_table_pathC[[i]] <- matrix(c(
    as.numeric(number_diff_and_pathC[i]),
    diff_abund - as.numeric(number_diff_and_pathC[i]),
    as.numeric(number_in_pathC[i]) - as.numeric(number_diff_and_pathC[i]),
    as.numeric(number_not_pathC[i] - (diff_abund - as.numeric(number_diff_and_pathC[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_pathC <- double()
for(i in 1:length(contingency_table_pathC)){
  fisher_p_pathC <- c(fisher_p_pathC, fisher.test(contingency_table_pathC[[i]])$p.value)
}

k_pathC_fisher <- k_pathC_count %>% 
  filter(C != is.na(C)) %>% 
  mutate(fisher_p_value = fisher_p_pathC) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))




######## KEGG using Pst database as background (lines 480 - 698) ################
######## **Running this will overwrite variables from previous chunk
abundance <- read_tsv("processed_data/total_processed_abundance_t_test.txt")

total_k_terms <- read_excel("raw_data/kegg_combined_total.xlsx") %>% 
  dplyr::select("gi_number", "accession", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

# there's some homology between Pst and Pf that is accounted for in the DB but not in the "combined total"
# innerjoin to eliminate homologous entries that will artifactually inflate k terms
db_total_k_terms <- read_tsv("raw_data/kegg_combined_DB_hits.txt") %>% 
  dplyr::select("accession", "ko")

total_k_terms <- inner_join(total_k_terms, db_total_k_terms, by = "accession") %>% 
  select("gi_number", "ko.x") %>% 
  rename(ko = ko.x)

# do the analysis with those that have significant interaction term, and are significantly different in complete
# media between the species
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_complete_p", "significant_comp", "comp_fc", "log2_comp_fc") %>% 
  inner_join(total_k_terms, by = "gi_number")

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# db_total_k_terms <- read_tsv("raw_data/kegg_combined_DB_hits.txt") %>% 
#   dplyr::select("accession", "ko") 

db_k_term_count <- db_total_k_terms %>% 
  group_by(ko) %>% 
  summarize(N=n())

combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

# make contingency_tables

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(grepl(TRUE, !is.na(kegg$ko))) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_comp) & grepl(combo_term_count[i,1], kegg$ko)))
  contingency_table[[i]] <- matrix(c(
    as.numeric(number_diff_and_ko[i]),
    diff_abund - as.numeric(number_diff_and_ko[i]),
    as.numeric(number_in_ko[i]) - as.numeric(number_diff_and_ko[i]),
    as.numeric(number_not_ko[i] - (diff_abund - as.numeric(number_diff_and_ko[i])))),
    ncol = 2)
}

# run Fisher's exact test

fisher_p <- double()
for(i in 1:length(contingency_table)){
  fisher_p <- c(fisher_p, fisher.test(contingency_table[[i]])$p.value)
}

k_term_fisher <- combo_term_count %>% 
  filter(ko != is.na(ko)) %>% 
  mutate(fisher_p_value = fisher_p) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

k_term_graph <- k_term_fisher %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10)

k_term_graph %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), x = N.x))+
  geom_col(fill = "gray26")+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  geom_text(aes(label = signif(adj_fisher_p, 3)), hjust = -0.1, size = 3)+
  #xlim(0,300)+
  theme_classic()

k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  #xlim(0,300)+
  theme_classic()


# use p. putida as comparison 
# https://academic.oup.com/mbe/article/28/1/483/982857#77836793
# https://www.mdpi.com/2073-4425/11/2/139
# https://link.springer.com/article/10.1007/s10658-017-1401-8

# convert GI to NCBI ID using:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=GIs separated by commas&rettype=acc

#just kidding...can upload FASTA file to blast KOALA https://www.kegg.jp/blastkoala/ and get list of k terms


k_term_nona <- total_k_terms %>% 
  mutate(mapped_to_k = grepl("", total_k_terms$ko)) %>% 
  filter(.$mapped_to_k == TRUE) %>% 
  dplyr::select(ko)

db_k_term_nona <- db_total_k_terms %>% 
  mutate(mapped_to_k = grepl("", db_total_k_terms$ko)) %>% 
  filter(.$mapped_to_k == TRUE) %>% 
  dplyr::select(ko)


kegg_db <- read_tsv("raw_data/kegg_hierarchy.txt")

k_term_for_db <- character()
for(i in 1:nrow(kegg_db)){
  k_term_for_db <- c(k_term_for_db, gsub(pattern = "(D      )(.*)(  .*)", replacement = "\\2", x = kegg_db[[i,4]]))
}

kegg_db <- mutate(kegg_db, ko = k_term_for_db)

kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  full_join(kegg, kegg_db, by = "ko")

db_kegg_with_path <- kegg_db %>% 
  filter(ko %in% db_k_term_nona$ko) %>% 
  full_join(db_k_term_nona, kegg_db, by = "ko")


## evaluate highest kegg level A
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

db_k_path_count <- db_kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(grepl(TRUE, !is.na(db_kegg_with_path$A))) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_comp) & grepl(combo_path_count[i,1], kegg_with_path$A)))
  contingency_table_path[[i]] <- matrix(c(
    as.numeric(number_diff_and_path[i]),
    diff_abund - as.numeric(number_diff_and_path[i]),
    as.numeric(number_in_path[i]) - as.numeric(number_diff_and_path[i]),
    as.numeric(number_not_path[i] - (diff_abund - as.numeric(number_diff_and_path[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_path <- double()
for(i in 1:length(contingency_table_path)){
  fisher_p_path <- c(fisher_p_path, fisher.test(contingency_table_path[[i]])$p.value)
}

k_path_fisher <- combo_path_count %>% 
  filter(A != is.na(A)) %>% 
  mutate(fisher_p_value = fisher_p_path) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

k_path_fisher %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), x = N.x))+
  geom_col(fill = "gray26")+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  geom_text(aes(label = signif(adj_fisher_p, 3)), hjust = -0.1, size = 3)+
  xlim(0,2000)+
  theme_classic()

k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  xlim(0,2000)+
  theme_classic()

## evaluate kegg level B
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

db_k_pathB_count <- db_kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(grepl(TRUE, !is.na(db_kegg_with_path$B))) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_comp) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
  contingency_table_pathB[[i]] <- matrix(c(
    as.numeric(number_diff_and_pathB[i]),
    diff_abund - as.numeric(number_diff_and_pathB[i]),
    as.numeric(number_in_pathB[i]) - as.numeric(number_diff_and_pathB[i]),
    as.numeric(number_not_pathB[i] - (diff_abund - as.numeric(number_diff_and_pathB[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_pathB <- double()
for(i in 1:length(contingency_table_pathB)){
  fisher_p_pathB <- c(fisher_p_pathB, fisher.test(contingency_table_pathB[[i]])$p.value)
}

k_pathB_fisher <- combo_pathB_count %>% 
  filter(B != is.na(B)) %>% 
  mutate(fisher_p_value = fisher_p_pathB) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

k_pathB_graph <- k_pathB_fisher %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15)

k_pathB_graph %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), x = N.x))+
  geom_col(fill = "gray26")+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  geom_text(aes(label = signif(adj_fisher_p, 3)), hjust = -0.1, size = 3)+
  xlim(0,1000)+
  theme_classic()

k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  xlim(0,1000)+
  theme_classic()

## evaluate kegg level C
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

db_k_pathC_count <- db_kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(grepl(TRUE, !is.na(db_kegg_with_path$C))) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_comp) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
  contingency_table_pathC[[i]] <- matrix(c(
    as.numeric(number_diff_and_pathC[i]),
    diff_abund - as.numeric(number_diff_and_pathC[i]),
    as.numeric(number_in_pathC[i]) - as.numeric(number_diff_and_pathC[i]),
    as.numeric(number_not_pathC[i] - (diff_abund - as.numeric(number_diff_and_pathC[i])))),
    ncol = 2)
}

# run Fisher's exact test
fisher_p_pathC <- double()
for(i in 1:length(contingency_table_pathC)){
  fisher_p_pathC <- c(fisher_p_pathC, fisher.test(contingency_table_pathC[[i]])$p.value)
}

k_pathC_fisher <- combo_pathC_count %>% 
  filter(C != is.na(C)) %>% 
  mutate(fisher_p_value = fisher_p_pathC) %>% 
  mutate(adj_fisher_p = p.adjust(fisher_p_value, method = "BH"))

k_pathC_graph <- k_pathC_fisher %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15)

k_pathC_graph %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), x = N.x))+
  geom_col(fill = "gray26")+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  geom_text(aes(label = signif(adj_fisher_p, 3)), hjust = -0.1, size = 3)+
  xlim(0,600)+
  theme_classic()

k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  xlim(0,800)+
  theme_classic()


# if george ever wants to look:
# 
# kegg_db <- read_tsv("raw_data/ko00001.keg")
# kegg_a <- character()
# kegg_b <- character()
# kegg_c <- character()
# k_term_with_path <- tibble(ncol(4), colnames(c("a", "b", "c", "ko")))
# kegg_a_for_df <- character()
# kegg_b_for_df <- character()
# kegg_c_for_df <- character()
# kegg_ko_for_df <- character()
# for(i in 1:nrow(kegg_db)){
#   case_when(grepl("^A") ~ kegg_a <- c(kegg_a, kegg_db[[i]]),
#             grepl("^B") ~ kegg_b <- c(kegg_b, kegg_db[[i]]),
#             grepl("^C") ~ kegg_c <- c(kegg_c, kegg_db[[i]]),
#             grepl("^D") ~ 
#               if(gsub(pattern = "(D      )(.*)(  .*)", replacement = "\\2", x = kegg_db[[i,1]]) %in% k_term_nona$ko)
#                  { kegg_a_for_df <- c(kegg_a_for_df, kegg_a[[nrow(kegg_a)]]);
#                     kegg_b_for_df <- c(kegg_b_for_df, kegg_b[[nrow(kegg_b)]]);
#                     kegg_c_for_df <- c(kegg_c_for_df, kegg_c[[nrow(kegg_c)]]);
#                     kegg_ko_for_df <- c(kegg_ko_for_df, gsub(pattern = "(D      )(.*)(  .*)", replacement = "\\2", x = kegg_db[[i,1]]))}
#   )
# }
# 
# 
# 
# for(i in 1:nrow(kegg_db)){
#   case_when(starts_with("A") ~ kegg_a <- c(kegg_a, kegg_db[[i]]),
#             starts_with("B") ~ kegg_b <- c(kegg_b, kegg_db[[i]]),
#             starts_with("C") ~ kegg_c <- c(kegg_c, kegg_db[[i]]),
#             starts_with("D") ~ 
#               if(gsub(pattern = "(D      )(.*)(  .*)", replacement = "\\2", x = kegg_db[[i,1]]) == k_term_nona[[i,1]]
#                  ~ k_term_with_path$a <- c(k_term_with_path$a, kegg_a[[i]]) &
#                  k_term_with_path$b <- c(k_term_with_path$b, kegg_b[[i]]) &
#                  k_term_with_path$c <- c(k_term_with_path$c, kegg_c[[i]]) &
#                  k_term_with_path$ko <- c(k_term_with_path$ko, k_term_nona[[i.1]])))
# }



################## Virulent Predict ########################
# sig interaction effect, enriched Pst minimal, higher Pf baseline
higher_vir_pred <- read_tsv("raw_data/virulentpred_enrich_Pst_higher_Pf_baseline.txt")

higher_vir_pred <- higher_vir_pred %>% 
  filter(.[,3] == "Virulent") %>% 
  arrange(-.[4])

write_tsv(higher_vir_pred, "processed_data/virulentpred_enrich_Pst_higher_Pf_baseline.tsv")

# sig interaction effect, enriched Pst minimal, same Pf baseline
same_vir_pred <- read_tsv("raw_data/virulentpred_enrich_Pst_same_Pf_baseline.txt")

same_vir_pred <- same_vir_pred %>% 
  filter(.[,3] == "Virulent") %>% 
  arrange(-.[4])

write_tsv(same_vir_pred, "processed_data/virulentpred_enrich_Pst_same_Pf_baseline.tsv")

# sig interaction effect, fc enriched Pst minimal, higher Pf baseline
fc_and_higher_vir_pred <- read_tsv("raw_data/virulentpred_fc_enrich_Pst_higher_Pf_baseline.txt")

fc_and_higher_vir_pred <- fc_and_higher_vir_pred %>% 
  filter(.[,3] == "Virulent") %>% 
  arrange(-.[4])

write_tsv(fc_and_higher_vir_pred, "processed_data/virulentpred_fc_enrich_Pst_higher_Pf_baseline.tsv")

# sig interaction effect, fc enriched Pst minimal, same Pf baseline
fc_and_same_vir_pred <- read_tsv("raw_data/virulentpred_fc_enrich_Pst_same_Pf_baseline.txt")

fc_and_same_vir_pred <- fc_and_same_vir_pred %>% 
  filter(.[,3] == "Virulent") %>% 
  arrange(-.[4])

write_tsv(fc_and_same_vir_pred, "processed_data/virulentpred_fc_enrich_Pst_same_Pf_baseline.tsv")


