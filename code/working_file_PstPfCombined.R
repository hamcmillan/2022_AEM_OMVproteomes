library(tidyverse)
library(dplyr)
library(readxl)
library(missForest)
library(EnhancedVolcano)


#Read abundance data, filter contaminants
total_abundance_go <- read_excel(path="raw_data/ProteinExpressionCombinedPstPf_RawAbundance.xlsx", 
                              sheet=5,
                              col_types=c(GI_Number = "numeric", Accession = "text",
                                          Description = "text", Contaminant = "logical",
                                          Coverage = "numeric", Num_Peptides = "numeric",
                                          Num_PSMs = "numeric", Num_Unique_Peptides = "numeric",
                                          Num_AAs = "numeric", MW_kda = "numeric",
                                          calc_pI = "numeric", Score_MSFragger = "numeric",
                                          Num_Peptides_by_MSFragger = "numeric",
                                          Num_Razor_Peptides = "numeric", 
                                          Pf_Complete_1 = "numeric", Pf_Complete_2 = "numeric",
                                          Pf_Complete_3 = "numeric", Pf_Minimal_1 = "numeric",
                                          Pf_Minimal_2 = "numeric", Pf_Minimal_3 = "numeric",
                                          Pf_Pool_1 = "numeric", Pf_Pool_2 = "numeric", 
                                          Pf_Pool_3 = "numeric", Pst_Complete_1 = "numeric", 
                                          Pst_Complete_2 = "numeric", Pst_Complete_3 = "numeric", 
                                          Pst_Minimal_1 = "numeric", Pst_Minimal_2 = "numeric", 
                                          Pst_Minimal_3 = "numeric", Pst_Pool_1 = "numeric", 
                                          Pst_Pool_2 = "numeric", Pst_Pool_3 = "numeric",
                                          Total_Pool_1 = "numeric", Total_Pool_2 = "numeric",
                                          Total_Pool_3 = "numeric")) %>% 
  rename_all(.funs=tolower) %>% 
  filter(contaminant == FALSE)

#how many missing values per row of each condition
pf_complete_na <- rowSums(is.na(dplyr::select(total_abundance_go, gi_number, description, pf_complete_1, pf_complete_2, pf_complete_3)))
pf_minimal_na <- rowSums(is.na(dplyr::select(total_abundance_go, gi_number, description, pf_minimal_1, pf_minimal_2, pf_minimal_3)))
pst_complete_na <- rowSums(is.na(dplyr::select(total_abundance_go, gi_number, description, pst_complete_1, pst_complete_2, pst_complete_3)))
pst_minimal_na <- rowSums(is.na(dplyr::select(total_abundance_go, gi_number, description, pst_minimal_1, pst_minimal_2, pst_minimal_3)))

#append these lists to main table
#if 2 or more missing values, add to new data table "differential"
total_abundance_go <- mutate(total_abundance_go, pf_comp_na = pf_complete_na, pf_min_na = pf_minimal_na,
                             pst_comp_na = pst_complete_na, pst_min_na = pst_minimal_na)
differential <- data.frame(filter(total_abundance_go, pf_comp_na >= 2 | pf_min_na >= 2 | pst_comp_na >= 2 | pst_min_na >= 2))

#save/write differentially expressed to new directory
write_tsv(differential, "processed_data/total_differential_abundance.tsv")

##  normalize abundance data
# test for normality
for(i in 15:35){
  print(shapiro.test(pull(total_abundance_go, i)))
}

# log2 transform, add transformed to df, "test" contains results of transformed shapiro
newname <- character()
test <- numeric()
for(i in 15:35){
  newname <- c(newname, paste("log_", colnames(total_abundance_go)[[i]], sep = "_"))
  total_abundance_go <- mutate(total_abundance_go,
                            "{newname[[i-14]]}" := log2(total_abundance_go[[i]]))
  test <- c(test, print(shapiro.test(pull(total_abundance_go, newname[[i-14]]))))
}

# pst_abundance_go <- mutate(pst_abundance_go, log_pst_complete_1 = log2(pst_complete_1))
# shapiro.test(pull(pst_abundance_go, "log_pst_complete_1"))

# pst_abundance_go <- mutate(pst_abundance_go, avg_pst_complete_1 = pst_complete_1/summary(pst_abundance_go$pst_complete_1)[["Mean"]])

# plot overlayed normalized distributions 
log_for_plot <- total_abundance_go %>% 
  dplyr::select(40:60) %>% 
  pivot_longer(cols = 1:20,names_to = "log_sample", values_to = "values")

log_for_plot %>% 
  ggplot(aes(x=values, color=log_sample)) +
  geom_density()

# mean-normalize the minimal media samples to their complete media condition, leave out pools
total_abundance_go <- mutate(total_abundance_go, avg_log_pf_complete = 
                            total_abundance_go %>% 
                            dplyr::select(40:42) %>% 
                            rowMeans())

total_abundance_go <- mutate(total_abundance_go, avg_log_pst_complete = 
                               total_abundance_go %>% 
                               dplyr::select(49:51) %>% 
                               rowMeans())

total_abundance_go <- mutate(total_abundance_go, mc_log_pf_complete_1 = total_abundance_go$log__pf_complete_1 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pf_complete_2 = total_abundance_go$log__pf_complete_2 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pf_complete_3 = total_abundance_go$log__pf_complete_3 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])

total_abundance_go <- mutate(total_abundance_go, mc_log_pf_minimal_1 = total_abundance_go$log__pf_minimal_1 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pf_minimal_2 = total_abundance_go$log__pf_minimal_2 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pf_minimal_3 = total_abundance_go$log__pf_minimal_3 / summary(total_abundance_go$avg_log_pf_complete)[["Mean"]])

total_abundance_go <- mutate(total_abundance_go, mc_log_pst_complete_1 = total_abundance_go$log__pst_complete_1 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pst_complete_2 = total_abundance_go$log__pst_complete_2 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pst_complete_3 = total_abundance_go$log__pst_complete_3 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])

total_abundance_go <- mutate(total_abundance_go, mc_log_pst_minimal_1 = total_abundance_go$log__pst_minimal_1 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pst_minimal_2 = total_abundance_go$log__pst_minimal_2 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])
total_abundance_go <- mutate(total_abundance_go, mc_log_pst_minimal_3 = total_abundance_go$log__pst_minimal_3 / summary(total_abundance_go$avg_log_pst_complete)[["Mean"]])

mc_log_for_plot <- total_abundance_go %>% 
  dplyr::select(63:74) %>% 
  pivot_longer(cols = 1:12, names_to = "mc_log_sample", values_to = "values")

mc_log_for_plot %>% 
  ggplot(aes(x=values, color=mc_log_sample))+
  geom_density()

#impute missing values
missing <- total_abundance_go %>% 
  dplyr::select(63:74) %>% 
  pivot_longer(cols = 1:12, names_to = "mc_log_sample", values_to = "values") %>% 
  mutate_all(~replace_na(.,0))

missing %>% 
  ggplot(aes(x=values, color=mc_log_sample))+
  geom_density()

hm_missing <- total_abundance_go %>% 
  dplyr::select(63:74) %>% 
  mutate_all(~replace_na(.,0))
color <- ifelse(hm_missing == 0, "white", "black")
heatmap(as.matrix(hm_missing), col = color, labRow = total_abundance_go$description)

# hm_missing <- pst_abundance_go %>% 
#   select(41:46) %>%
#   mutate_all(~replace(.,!is.na(.),1)) %>% 
#   mutate_all(~replace_na(.,0))
# heatmap(as.matrix(hm_missing))  


blanks_total_abund_go <- total_abundance_go %>% 
  dplyr::select(63:74) %>% 
  as.matrix()

imp_total_abund_go <- missForest(blanks_total_abund_go)

imp_for_plot <- imp_total_abund_go$ximp %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "imp_mc_log_sample", values_to = "values")

imp_for_plot %>% 
  ggplot(aes(x=values, color=imp_mc_log_sample))+
  geom_density()

# add imputed data to dataframe, test for differential expression
imp_for_tests <- as.data.frame(imp_total_abund_go$ximp)

total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_complete_1 = imp_for_tests$mc_log_pf_complete_1)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_complete_2 = imp_for_tests$mc_log_pf_complete_2)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_complete_3 = imp_for_tests$mc_log_pf_complete_3)

total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_minimal_1 = imp_for_tests$mc_log_pf_minimal_1)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_minimal_2 = imp_for_tests$mc_log_pf_minimal_2)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pf_minimal_3 = imp_for_tests$mc_log_pf_minimal_3)

total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_complete_1 = imp_for_tests$mc_log_pst_complete_1)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_complete_2 = imp_for_tests$mc_log_pst_complete_2)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_complete_3 = imp_for_tests$mc_log_pst_complete_3)

total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_minimal_1 = imp_for_tests$mc_log_pst_minimal_1)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_minimal_2 = imp_for_tests$mc_log_pst_minimal_2)
total_abundance_go <- mutate(total_abundance_go, imp_mc_log_pst_minimal_3 = imp_for_tests$mc_log_pst_minimal_3)


# test for differential expression
# helpful ANOVA and Tukey links: 
# https://www.scribbr.com/statistics/two-way-anova/
# https://www.scribbr.com/statistics/anova-in-r/
# https://stats.stackexchange.com/questions/31547/how-to-obtain-the-results-of-a-tukey-hsd-post-hoc-test-in-a-table-showing-groupe
# https://online.stat.psu.edu/stat485/lesson/12/12.7
#imp_for_anova <- as.matrix(imp_for_tests)
imp_for_anova <- as_tibble(imp_for_tests)


species <- c("pf", "pf", "pf", "pf", "pf", "pf", "pst", "pst", "pst", "pst", "pst", "pst")
media <- c("complete", "complete", "complete", "minimal", "minimal", "minimal", "complete", "complete", "complete", "minimal", "minimal", "minimal")
test_table <- tibble()
p_species <- double()
p_media <- double()
p_speciesxmedia <- double()
for(i in 1:nrow(imp_for_anova)){
  test_table <- imp_for_anova[i,] %>% 
    pivot_longer(cols = everything(), names_to = "imp_mc_log_sample", values_to = "values") %>% 
    mutate(species = species, media = media)
  interaction <- summary(aov(values ~ species*media, data = test_table))
  p_species <- c(p_species, interaction[[1]][["Pr(>F)"]][[1]])
  p_media <- c(p_media, interaction[[1]][["Pr(>F)"]][[2]])
  p_speciesxmedia <- c(p_speciesxmedia, interaction[[1]][["Pr(>F)"]][[3]])
}

# add p-values to dataframe, calculate and add adjusted p-values
# add logical column, TRUE = adj significant, FALSE = adj not significant based on interaction term
# only move forward with those that have a significant interaction effect, ie respond differently to media 
# in Pst than they do in Pf
total_abundance_go <- mutate(total_abundance_go, 
                             aov_p_species = p_species,
                             aov_p_media = p_media,
                             aov_p_speciesxmedia = p_speciesxmedia)
total_abundance_go <- total_abundance_go %>% 
  mutate(adj_aov_p_species = p.adjust(.$aov_p_species, method = "BH"),
         adj_aov_p_media = p.adjust(.$aov_p_media, method = "BH"),
         adj_aov_p_speciesxmedia = p.adjust(.$aov_p_speciesxmedia, method = "BH"))

total_abundance_go <- total_abundance_go %>% 
  mutate(significant_interaction = if_else(adj_aov_p_speciesxmedia <= 0.05, "TRUE", "FALSE")) %>%
  mutate(significant_interaction = as.logical(significant_interaction))

# volcano plot
# calculate fold change in each species and plot - there will likely be different enriched
# proteins based on the 2-way ANOVA than there were in the t-test
avg_imp_mc_log_pf_complete <- double()
avg_imp_mc_log_pf_minimal <- double()
avg_imp_mc_log_pst_complete <- double()
avg_imp_mc_log_pst_minimal <- double()
for(i in 1:nrow(total_abundance_go)){
  pf_complete = c(total_abundance_go[i, 75:77])
  pf_minimal = c(total_abundance_go[i, 78:80])
  pst_complete = c(total_abundance_go[i, 81:83])
  pst_minimal = c(total_abundance_go[i, 84:86])
  avg_imp_mc_log_pf_complete <- c(avg_imp_mc_log_pf_complete, print(mean(as.double(pf_complete))))
  avg_imp_mc_log_pf_minimal <- c(avg_imp_mc_log_pf_minimal, print(mean(as.double(pf_minimal))))
  avg_imp_mc_log_pst_complete <- c(avg_imp_mc_log_pst_complete, print(mean(as.double(pst_complete))))
  avg_imp_mc_log_pst_minimal <- c(avg_imp_mc_log_pst_minimal, print(mean(as.double(pst_minimal))))
}

total_abundance_go <- mutate(total_abundance_go, avg_imp_mc_log_pf_complete = avg_imp_mc_log_pf_complete)
total_abundance_go <- mutate(total_abundance_go, avg_imp_mc_log_pf_minimal = avg_imp_mc_log_pf_minimal)
total_abundance_go <- mutate(total_abundance_go, avg_imp_mc_log_pst_complete = avg_imp_mc_log_pst_complete)
total_abundance_go <- mutate(total_abundance_go, avg_imp_mc_log_pst_minimal = avg_imp_mc_log_pst_minimal)

vp <- total_abundance_go %>% 
  dplyr::select("accession", contains("imp_mc_log"), "adj_aov_p_speciesxmedia", "significant_interaction") %>% 
  mutate(pf_fc = avg_imp_mc_log_pf_minimal/avg_imp_mc_log_pf_complete) %>% 
  mutate(pst_fc = avg_imp_mc_log_pst_minimal/avg_imp_mc_log_pst_complete) %>% 
  mutate(log_p_adj = -log10(adj_aov_p_speciesxmedia)) %>% 
  mutate(log2_pf_fc = log2(pf_fc)) %>% 
  mutate(log2_pst_fc = log2(pst_fc))

total_abundance_go <- mutate(total_abundance_go, pf_fc = vp$pf_fc)
total_abundance_go <- mutate(total_abundance_go, pst_fc = vp$pst_fc)
total_abundance_go <- mutate(total_abundance_go, nlog10_p_adj = vp$log_p_adj)
total_abundance_go <- mutate(total_abundance_go, log2_pf_fc = vp$log2_pf_fc)
total_abundance_go <- mutate(total_abundance_go, log2_pst_fc = vp$log2_pst_fc)

# vp %>% 
#   ggplot(aes(x = log2_fc, y = log_p_adj, color = significant_comp_vs_min))+
#   geom_point()

# pf volcano with 2-way anova values
EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_pf_fc", 
                y = "adj_aov_p_speciesxmedia",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 5),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)

# pst volcano with 2-way anova values
EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_pst_fc", 
                y = "adj_aov_p_speciesxmedia",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 5),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)

# run tukey post-hoc test on data to see which groups are different for each protein
library(agricolae)

imp_for_anova <- as_tibble(imp_for_tests)


species <- c("pf", "pf", "pf", "pf", "pf", "pf", "pst", "pst", "pst", "pst", "pst", "pst")
media <- c("complete", "complete", "complete", "minimal", "minimal", "minimal", "complete", "complete", "complete", "minimal", "minimal", "minimal")
test_table <- tibble()
tukey_p_pst_comp_v_min <- double()
connecting_letters <- tibble()
for(i in 1:nrow(imp_for_anova)){
  test_table <- imp_for_anova[i,] %>% 
    pivot_longer(cols = everything(), names_to = "imp_mc_log_sample", values_to = "values") %>% 
    mutate(species = species, media = media)
  test_term <- with(test_table, interaction(species, media))
  test_anova <- aov(values ~ test_term, data = test_table)
  test_tukey <- HSD.test(test_anova, "test_term", group = TRUE)
  test_tukey2 <- TukeyHSD(test_anova)
  tukey_p_pst_comp_v_min <- c(tukey_p_pst_comp_v_min, test_tukey2$test_term[5,4])
  connecting_letters <- rbind(connecting_letters, test_tukey$groups$groups)
}

colnames(connecting_letters)<-c("group_pst_min", "group_pst_comp", "group_pf_comp", "group_pf_min")

# add tukey p values and connecting letters to the full data frame
total_abundance_go <- mutate(total_abundance_go, tukey_p_pst_comp_v_min = tukey_p_pst_comp_v_min)
total_abundance_go <- mutate(total_abundance_go, adj_tukey_p_pst_comp_v_min = 
                               p.adjust(total_abundance_go$tukey_p_pst_comp_v_min, method = "BH"))
total_abundance_go <- mutate(total_abundance_go, 
                             group_pst_min = connecting_letters$group_pst_min, 
                             group_pst_comp = connecting_letters$group_pst_comp, 
                             group_pf_comp = connecting_letters$group_pf_comp, 
                             group_pf_min = connecting_letters$group_pf_min)


# write new file with complete abundance, go terms, and significance testing
write_tsv(total_abundance_go, "processed_data/total_processed_abundance_2wayANOVA.txt")

# get list of significant interaction by 2 way ANOVA
aov_sigdif <- total_abundance_go %>% 
  dplyr::select("gi_number", "accession", "description", contains("imp_mc_log"), "adj_aov_p_speciesxmedia", 
                "significant_interaction", "log2_pf_fc", "log2_pst_fc") %>% 
  filter(significant_interaction == TRUE)

aov_sigenrich_pf_min <- aov_sigdif %>% 
  filter(avg_imp_mc_log_pf_minimal > avg_imp_mc_log_pf_complete)
aov_sigenrich_pf_comp <- aov_sigdif %>% 
  filter(avg_imp_mc_log_pf_complete > avg_imp_mc_log_pf_minimal)
enrich_and_fc_pf_min <- aov_sigenrich_pf_min %>% 
  filter(log2_pf_fc >= 0.5)
enrich_and_fc_pf_comp <- aov_sigenrich_pf_comp %>% 
  filter(abs(log2_pf_fc) >= 0.5)

aov_sigenrich_pst_min <- aov_sigdif %>% 
  filter(avg_imp_mc_log_pst_minimal > avg_imp_mc_log_pst_complete)
aov_sigenrich_pst_comp <- aov_sigdif %>% 
  filter(avg_imp_mc_log_pst_complete > avg_imp_mc_log_pst_minimal)
enrich_and_fc_pst_min <- aov_sigenrich_pst_min %>% 
  filter(log2_pst_fc >= 0.5)
enrich_and_fc_pst_comp <- aov_sigenrich_pst_comp %>% 
  filter(abs(log2_pst_fc) >= 0.5)

#save/write differentially expressed to new directory
write_tsv(aov_sigenrich_pf_comp, "processed_data/pf_sigenrich_comp_2wayANOVA_norm2comp.tsv")
write_tsv(aov_sigenrich_pf_min, "processed_data/pf_sigenrich_min_2wayANOVA_norm2comp.tsv")
write_tsv(enrich_and_fc_pf_comp, "processed_data/pf_enrich_and_fc_comp_2wayANOVA_norm2comp.tsv")
write_tsv(enrich_and_fc_pf_min, "processed_data/pf_enrich_and_fc_min_2wayANOVA_norm2comp.tsv")

write_tsv(aov_sigenrich_pst_comp, "processed_data/pst_sigenrich_comp_2wayANOVA_norm2comp.tsv")
write_tsv(aov_sigenrich_pst_min, "processed_data/pst_sigenrich_min_2wayANOVA_norm2comp.tsv")
write_tsv(enrich_and_fc_pst_comp, "processed_data/pst_enrich_and_fc_comp_2wayANOVA_norm2comp.tsv")
write_tsv(enrich_and_fc_pst_min, "processed_data/pst_enrich_and_fc_min_2wayANOVA_norm2comp.tsv")


# lists of interesting categories within the significantly enriched in minimal set
# sig_pathogenesis <- pf_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("pathogenesis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_photosynthesis <- pf_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("photosynthesis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_iron <- pf_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("iron", .$description)) %>% 
#   arrange(adj_p_value_comp_vs_min)
# 
# sig_symbiosis <- pf_abundance_go %>% 
#   filter(significant_comp_vs_min == TRUE) %>% 
#   filter(grepl("symbiosis", .$biological_process)) %>% 
#   arrange(adj_p_value_comp_vs_min)

# heatmap to explore data
hm_total <- total_abundance_go %>% 
  dplyr::select(40:45, 49:54)

rownames(hm_total) <- total_abundance_go$accession

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
hm_grouped <- select(total_abundance_go, "avg_imp_mc_log_pf_complete", "avg_imp_mc_log_pf_minimal",
                     "avg_imp_mc_log_pst_complete", "avg_imp_mc_log_pst_minimal")
rownames(hm_grouped) <- total_abundance_go$accession
hm_grouped <- scale(as.matrix(hm_grouped))
pheatmap(hm_grouped, show_rownames = FALSE, cluster_cols = TRUE, cutree_rows = 10, cellwidth = 75)

hm_grouped_clust <- pheatmap(hm_grouped, show_rownames = FALSE, cluster_cols = TRUE, cutree_rows = 10, cellwidth = 75)
hm_grouped_clust_list <- cbind(hm_grouped, cluster = cutree(hm_grouped_clust$tree_row, k = 10))


# PCA to explore data
#library(limma)

total_pca <- total_abundance_go %>%
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


############### within sample significance testing #################
# nothing showed up as significant between samples with t-test, so we will
# use a z score instead to ask what is significantly more abundant in each condition,
# follow up with qualitative comparisons between samples and species


# median-normalize for each condition, leave out the pooled samples

pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_complete_1 = pf_abundance_go$log__pf_complete_1 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])
pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_complete_2 = pf_abundance_go$log__pf_complete_2 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])
pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_complete_3 = pf_abundance_go$log__pf_complete_3 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])

pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_minimal_1 = pf_abundance_go$log__pf_minimal_1 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])
pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_minimal_2 = pf_abundance_go$log__pf_minimal_2 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])
pf_abundance_go <- mutate(pf_abundance_go, mednorm_log_pf_minimal_3 = pf_abundance_go$log__pf_minimal_3 / summary(pf_abundance_go$avg_log_pf_complete)[["Median"]])

mednorm_log_for_plot <- pf_abundance_go %>% 
  dplyr::select(61:66) %>% 
  pivot_longer(cols = 1:6, names_to = "mednorm_log_sample", values_to = "values")

mednorm_log_for_plot %>% 
  ggplot(aes(x=values, color=mednorm_log_sample))+
  geom_density()

#impute missing values
missing <- pf_abundance_go %>% 
  dplyr::select(61:66) %>% 
  pivot_longer(cols = 1:6, names_to = "mednorm_log_sample", values_to = "values") %>% 
  mutate_all(~replace_na(.,0))

missing %>% 
  ggplot(aes(x=values, color=mednorm_log_sample))+
  geom_density()

hm_missing <- pf_abundance_go %>% 
  dplyr::select(61:66) %>% 
  mutate_all(~replace_na(.,0))
color <- ifelse(hm_missing == 0, "white", "black")
heatmap(as.matrix(hm_missing), col = color, labRow = pf_abundance_go$description)

# hm_missing <- pst_abundance_go %>% 
#   select(41:46) %>%
#   mutate_all(~replace(.,!is.na(.),1)) %>% 
#   mutate_all(~replace_na(.,0))
# heatmap(as.matrix(hm_missing))  


blanks_pf_abund_go <- pf_abundance_go %>% 
  dplyr::select(61:66) %>% 
  as.matrix()

imp_pf_abund_go <- missForest(blanks_pf_abund_go)

imp_for_plot <- imp_pf_abund_go$ximp %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "imp_mednorm_log_sample", values_to = "values")

imp_for_plot %>% 
  ggplot(aes(x=values, color=imp_mednorm_log_sample))+
  geom_density()

# add imputed data to dataframe, test for differential expression
imp_for_tests <- as.data.frame(imp_pf_abund_go$ximp)

pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_complete_1 = imp_for_tests$mednorm_log_pf_complete_1)
pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_complete_2 = imp_for_tests$mednorm_log_pf_complete_2)
pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_complete_3 = imp_for_tests$mednorm_log_pf_complete_3)

pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_minimal_1 = imp_for_tests$mednorm_log_pf_minimal_1)
pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_minimal_2 = imp_for_tests$mednorm_log_pf_minimal_2)
pf_abundance_go <- mutate(pf_abundance_go, imp_mednorm_log_pf_minimal_3 = imp_for_tests$mednorm_log_pf_minimal_3)

temp <- pf_abundance_go

# z scores for each column
z_scores <- pf_abundance_go %>% 
  select(67:72) %>% # 31-39 log_pf samples, 41-46 mc_log samples, 47-52 imp_mc_log samples, 61-66 mednorm_log samples, 67-72 imp_mednorm_log samples
  as_data_frame()

z_scores <- sapply(z_scores, function(z_scores) ((z_scores-median(z_scores, na.rm = TRUE))/sd(z_scores, na.rm = TRUE)))

z_scores <- rename(as_data_frame(z_scores),
                   z_imp_mednorm_log_pf_complete_1 = imp_mednorm_log_pf_complete_1,
                   z_imp_mednorm_log_pf_complete_2 = imp_mednorm_log_pf_complete_2,
                   z_imp_mednorm_log_pf_complete_3 = imp_mednorm_log_pf_complete_3,
                   z_imp_mednorm_log_pf_minimal_1 = imp_mednorm_log_pf_minimal_1,
                   z_imp_mednorm_log_pf_minimal_2 = imp_mednorm_log_pf_minimal_2,
                   z_imp_mednorm_log_pf_minimal_3 = imp_mednorm_log_pf_minimal_3,
)

pf_abundance_go <- cbind(pf_abundance_go, z_scores)

# p-value for each z-score
# https://replicationindex.com/2019/03/31/one-tail-or-two-tails-that-is-the-question/
# https://www.statology.org/p-value-of-z-score-r/
zp <- z_scores

zp <- sapply(zp, function(zp) ((pnorm(-abs(zp)))*2))

zp <- rename(as_data_frame(zp),
             zp_imp_mednorm_log_pf_complete_1 = z_imp_mednorm_log_pf_complete_1,
             zp_imp_mednorm_log_pf_complete_2 = z_imp_mednorm_log_pf_complete_2,
             zp_imp_mednorm_log_pf_complete_3 = z_imp_mednorm_log_pf_complete_3,
             zp_imp_mednorm_log_pf_minimal_1 = z_imp_mednorm_log_pf_minimal_1,
             zp_imp_mednorm_log_pf_minimal_2 = z_imp_mednorm_log_pf_minimal_2,
             zp_imp_mednorm_log_pf_minimal_3 = z_imp_mednorm_log_pf_minimal_3,
)

pf_abundance_go <- cbind(pf_abundance_go, zp)

# make list of significant p-values in at least 2 of 3 replicates for each condition

#how many significant values per row of each condition
pf_comp_sig <- rowSums(zp[1:3] <= 0.05, na.rm = TRUE)
pf_min_sig <- rowSums(zp[4:6] <= 0.05, na.rm = TRUE)

#append these lists to main table
#if 2 or more significant values, add to new data table "sig_within"
pf_abundance_go <- cbind(pf_abundance_go, pf_comp_sig)
pf_abundance_go <- cbind(pf_abundance_go, pf_min_sig)
sig_within <- data.frame(filter(pf_abundance_go, pf_comp_sig >= 2 | pf_min_sig >= 2))

#save/write differentially expressed to new directory
write_tsv(sig_within, "processed_data/pf_sig_z_abundance.tsv")

# volcano plot
# 31-39 log_pf samples, 41-46 mc_log samples, 47-52 imp_mc_log samples, 61-66 mednorm_log samples, 67-72 imp_mednorm_log samples
avg_imp_mednorm_log_complete <- double()
avg_imp_mednorm_log_minimal <- double()
for(i in 1:nrow(pf_abundance_go)){
  complete = c(pf_abundance_go[i, 67:69])
  minimal = c(pf_abundance_go[i, 70:72])
  avg_imp_mednorm_log_complete <- c(avg_imp_mednorm_log_complete, print(mean(as.double(complete))))
  avg_imp_mednorm_log_minimal <- c(avg_imp_mednorm_log_minimal, print(mean(as.double(minimal))))
}

pf_abundance_go <- mutate(pf_abundance_go, avg_imp_mednorm_log_complete = avg_imp_mednorm_log_complete)
pf_abundance_go <- mutate(pf_abundance_go, avg_imp_mednorm_log_minimal = avg_imp_mednorm_log_minimal)

# compute z scores and p values for averaged data
pf_abundance_go <- pf_abundance_go %>% 
  mutate(z_avg_imp_mednorm_log_complete = (avg_imp_mednorm_log_complete-median(avg_imp_mednorm_log_complete))/sd(avg_imp_mednorm_log_complete)) %>% 
  mutate(zp_avg_imp_mednorm_log_complete = 2* pnorm(-abs(z_avg_imp_mednorm_log_complete))) %>% #mean = mean(z_avg_imp_mednorm_log_complete), sd = sd(z_avg_imp_mednorm_log_complete)
  mutate(z_avg_imp_mednorm_log_minimal = (avg_imp_mednorm_log_minimal-median(avg_imp_mednorm_log_minimal))/sd(avg_imp_mednorm_log_minimal)) %>% 
  mutate(zp_avg_imp_mednorm_log_minimal = 2* pnorm(-abs(z_avg_imp_mednorm_log_minimal))) #mean = mean(z_avg_imp_mednorm_log_minimal), sd = sd(z_avg_imp_mednorm_log_minimal)

vp <- pf_abundance_go %>% 
  mutate(fc_comp = avg_imp_mednorm_log_complete/median(avg_imp_mednorm_log_complete)) %>% 
  mutate(log_p_comp = -log10(zp_avg_imp_mednorm_log_complete)) %>% 
  mutate(log2_fc_comp = log2(fc_comp)) %>% 
  mutate(fc_min = avg_imp_mednorm_log_minimal/median(avg_imp_mednorm_log_minimal)) %>% 
  mutate(log_p_min = -log10(zp_avg_imp_mednorm_log_minimal)) %>% 
  mutate(log2_fc_min = log2(fc_min))

pf_abundance_go <- mutate(pf_abundance_go, fc_comp_vs_median = vp$fc_comp)
pf_abundance_go <- mutate(pf_abundance_go, nlog10_p_avg_mednorm_comp = vp$log_p_comp)
pf_abundance_go <- mutate(pf_abundance_go, log2_fc_comp_vs_median = vp$log2_fc_comp)

pf_abundance_go <- mutate(pf_abundance_go, fc_min_vs_median = vp$fc_min)
pf_abundance_go <- mutate(pf_abundance_go, nlog10_p_avg_mednorm_min = vp$log_p_min)
pf_abundance_go <- mutate(pf_abundance_go, log2_fc_min_vs_median = vp$log2_fc_min)


# vp %>% 
#   ggplot(aes(x = log2_fc, y = log_p_adj, color = significant_comp_vs_min))+
#   geom_point()

# this is less useful since we're not comparing two samples
EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_fc_comp", 
                y = "zp_avg_imp_mednorm_log_complete",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 4),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)

EnhancedVolcano(vp, 
                lab = vp$accession, 
                x = "log2_fc_min", 
                y = "zp_avg_imp_mednorm_log_minimal",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                xlim = c(-1, 1),
                ylim = c(0, 4),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                #labhjust = c(1.0:18),
                #labvjust = c(1.1:1.28)
)

# perhaps more useful, just visualization of ordered list
vp_graph <- vp %>% 
  filter(pf_comp_sig >= 2) %>% 
  arrange(zp_avg_imp_mednorm_log_complete) %>% 
  slice(1:50)
vp_graph %>% 
  ggplot(aes(x=zp_avg_imp_mednorm_log_complete, y=reorder(description, fc_comp)))+
  geom_col(fill = "gray26")+
  labs(x= "P-value (Fold Change vs Median)",
       y= "Protein",
       title = "Pf Complete")+
  geom_text(aes(label = signif(fc_comp, 3)), hjust = -0.1, size = 3)+
  #xlim(0,0.003)+
  theme_classic()

vp_graph <- vp %>% 
  filter(pf_min_sig >= 2) %>% 
  arrange(zp_avg_imp_mednorm_log_minimal) %>% 
  slice(1:50)
vp_graph %>% 
  ggplot(aes(x=zp_avg_imp_mednorm_log_minimal, y=reorder(description, fc_min)))+
  geom_col(fill = "gray26")+
  labs(x= "Log10 P-value (Fold Change vs Median)",
       y= "Protein",
       title = "Pf Minimal")+
  geom_text(aes(label = signif(fc_min, 3)), hjust = -0.1, size = 3)+
  #xlim(0,0.004)+
  theme_classic()


# write new file with complete abundance, go terms, and significance testing
write_tsv(pf_abundance_go, "processed_data/pf_processed_abundance_go.txt")

# get list of significant difference by z score comp and min vs median
z_sigdif <- pf_abundance_go %>% 
  dplyr::select("gi_number", "accession", "description", contains("imp_mednorm_log_pf"), 
                "pf_comp_sig", "pf_min_sig", contains("avg_imp_mednorm"), "fc_comp_vs_median", 
                contains("nlog10_p_avg_mednorm"), "log2_fc_comp_vs_median", "fc_min_vs_median", 
                "log2_fc_min_vs_median") %>% 
  filter(pf_comp_sig >= 2 | pf_min_sig >= 2)

z_sigenrich_min <- z_sigdif %>% 
  filter(avg_imp_mednorm_log_minimal > median(avg_imp_mednorm_log_minimal))
z_sigenrich_comp <- z_sigdif %>% 
  filter(avg_imp_mednorm_log_complete > median(avg_imp_mednorm_log_complete))
z_enrich_and_fc_min <- z_sigenrich_min %>% 
  filter(log2_fc_min_vs_median >= 0.5)
z_enrich_and_fc_comp <- z_sigenrich_comp %>% 
  filter((log2_fc_comp_vs_median) >= 0.5)

# write to files
write_tsv(z_sigenrich_comp, "processed_data/pf_z_sigenrich_comp.tsv")
write_tsv(z_sigenrich_min, "processed_data/pf_z_sigenrich_min.tsv")
write_tsv(z_enrich_and_fc_comp, "processed_data/pf_z_enrich_and_fc_comp.tsv")
write_tsv(z_enrich_and_fc_min, "processed_data/pf_z_enrich_and_fc_min.tsv")


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
total_k_terms <- read_excel("raw_data/kegg_Pf_total.xlsx") %>% 
  dplyr::select("gi_number", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

kegg <- pf_abundance_go %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pf_complete_1", 
                "imp_mednorm_log_pf_complete_2", "imp_mednorm_log_pf_complete_3",
                "imp_mednorm_log_pf_minimal_1", "imp_mednorm_log_pf_minimal_2",
                "imp_mednorm_log_pf_minimal_3", "avg_imp_mednorm_log_complete",
                "avg_imp_mednorm_log_minimal", "zp_avg_imp_mednorm_log_complete",
                "zp_avg_imp_mednorm_log_minimal", "fc_comp_vs_median", "fc_min_vs_median",
                "log2_fc_comp_vs_median", "log2_fc_min_vs_median", "pf_comp_sig", "pf_min_sig") %>% 
  mutate(significant_comp = if_else(pf_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pf_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  mutate(significant_fc_comp = as.logical(if_else(significant_comp == TRUE & fc_comp_vs_median > 1, "TRUE", "FALSE"))) %>% 
  mutate(significant_fc_min = as.logical(if_else(significant_min == TRUE & fc_min_vs_median > 1, "TRUE", "FALSE"))) %>% 
  inner_join(total_k_terms, by = "gi_number")

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make contingency_tables

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(k_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, k_term_count[i,2])
  number_not_ko <- c(number_not_ko, sum(grepl(TRUE, !is.na(kegg$ko))) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_fc_comp) & grepl(k_term_count[i,1], kegg$ko)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(k_path_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, k_path_count[i,2])
  number_not_path <- c(number_not_path, sum(grepl(TRUE, !is.na(kegg_with_path$A))) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(k_path_count[i,1], kegg_with_path$ko)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(k_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, k_pathB_count[i,2])
  number_not_pathB <- c(number_not_pathB, sum(grepl(TRUE, !is.na(kegg_with_path$B))) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(k_pathB_count[i,1], kegg_with_path$ko)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(k_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, k_pathC_count[i,2])
  number_not_pathC <- c(number_not_pathC, sum(grepl(TRUE, !is.na(kegg_with_path$C))) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(k_pathC_count[i,1], kegg_with_path$ko)))
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

total_k_terms <- read_excel("raw_data/kegg_Pf_total.xlsx") %>% 
  dplyr::select("gi_number", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

kegg <- pf_abundance_go %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pf_complete_1", 
                "imp_mednorm_log_pf_complete_2", "imp_mednorm_log_pf_complete_3",
                "imp_mednorm_log_pf_minimal_1", "imp_mednorm_log_pf_minimal_2",
                "imp_mednorm_log_pf_minimal_3", "avg_imp_mednorm_log_complete",
                "avg_imp_mednorm_log_minimal", "zp_avg_imp_mednorm_log_complete",
                "zp_avg_imp_mednorm_log_minimal", "fc_comp_vs_median", "fc_min_vs_median",
                "log2_fc_comp_vs_median", "log2_fc_min_vs_median", "pf_comp_sig", "pf_min_sig") %>% 
  mutate(significant_comp = if_else(pf_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pf_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  mutate(significant_fc_comp = as.logical(if_else(significant_comp == TRUE & fc_comp_vs_median > 1, "TRUE", "FALSE"))) %>% 
  mutate(significant_fc_min = as.logical(if_else(significant_min == TRUE & fc_min_vs_median > 1, "TRUE", "FALSE"))) %>% 
  inner_join(total_k_terms, by = "gi_number")

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

db_total_k_terms <- read_tsv("raw_data/kegg_Pf_db_hits.txt") %>% 
  dplyr::select("accession", "ko") 

db_k_term_count <- db_total_k_terms %>% 
  group_by(ko) %>% 
  summarize(N=n())

combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

# make contingency_tables
######## sig fc comp first ###########
diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(grepl(TRUE, !is.na(kegg$ko))) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_fc_comp) & grepl(combo_term_count[i,1], kegg$ko)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(grepl(TRUE, !is.na(db_kegg_with_path$A))) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(combo_path_count[i,1], kegg_with_path$A)))
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
  xlim(0,1500)+
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
  xlim(0,1500)+
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(grepl(TRUE, !is.na(db_kegg_with_path$B))) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(grepl(TRUE, !is.na(db_kegg_with_path$C))) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_fc_comp) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
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
  xlim(0,400)+
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
  #xlim(0,800)+
  theme_classic()

# make contingency_tables
######## sig fc min now ###########
diff_abund <- sum(grepl(TRUE, kegg$significant_fc_min))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(grepl(TRUE, !is.na(kegg$ko))) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_fc_min) & grepl(combo_term_count[i,1], kegg$ko)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_min))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(grepl(TRUE, !is.na(db_kegg_with_path$A))) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_fc_min) & grepl(combo_path_count[i,1], kegg_with_path$A)))
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
  xlim(0,1500)+
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
  xlim(0,1500)+
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_min))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(grepl(TRUE, !is.na(db_kegg_with_path$B))) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_fc_min) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
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

diff_abund <- sum(grepl(TRUE, kegg$significant_fc_min))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(grepl(TRUE, !is.na(db_kegg_with_path$C))) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_fc_min) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
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
  xlim(0,400)+
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
  #xlim(0,800)+
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

vir_pred <- read_tsv("raw_data/virulentpred_Pst_enriched_minimal.txt")

vir_pred <- vir_pred %>% 
  filter(.[,3] == "Virulent") %>% 
  arrange(-.[4])

vir_pred_graph <- vir_pred %>% 
  slice(1:15)

write_tsv(vir_pred, "processed_data/pst_virpred_enriched_min.tsv")
