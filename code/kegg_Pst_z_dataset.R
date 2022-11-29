library(tidyverse)
library(dplyr)
library(readxl)

######## needed for all comparisons ########
db_total_k_terms <- read_tsv("raw_data/kegg_Pst_db_hits.txt") %>% 
  dplyr::select("accession", "ko")

db_k_term_count <- db_total_k_terms %>% 
  group_by(ko) %>% 
  summarize(N=n())

total_k_terms <- read_excel("raw_data/kegg_Pst_total.xlsx") %>% 
  dplyr::select("gi_number", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

# use KEGG hierarchy to assign pathways to each term
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

db_kegg_with_path <- kegg_db %>% 
  filter(ko %in% db_k_term_nona$ko) %>% 
  full_join(db_k_term_nona, kegg_db, by = "ko")

db_k_path_count <- db_kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

db_k_pathB_count <- db_kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

db_k_pathC_count <- db_kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

######## KEGG - sig by z score, all sig, sig in at least two of three replicates ################
abundance <- read_tsv("processed_data/pst_sig_z_abundance.tsv")

###### significant in complete ######

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>% 
  filter(significant_comp == TRUE)

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
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

kterm_sig_z_Pstcomp <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sig_z_Pstcomp_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
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

pathA_sig_z_Pstcomp <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sig_z_Pstcomp_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
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

pathB_sig_z_Pstcomp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sig_z_Pstcomp_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_z_Pstcomp <- temp %>% 
  add_column(category = NA, .before = "term") %>% 
  replace_na(list(category = "Level B"))


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
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

pathC_sig_z_Pstcomp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>%
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sig_z_Pstcomp_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>%
  slice(1:15)

kegg_z_Pstcomp <- add_row(kegg_z_Pstcomp, temp) %>% 
  replace_na(list(category = "Level C"))

# graph KEGG B and C in one plot
keggBC_sig_z_Pstcomp <- kegg_z_Pstcomp %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(term, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Pathway Enrichment - Level B and C",
       subtitle = "Combined Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Includes Significantly Higher and Lower than the median",
       y= "KEGG Pathway")+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(fct_relevel(category, "Level B", "Level C") ~ ., scales = "free_y", space = "free")+
  #xlim(0,550)+
  theme_classic()+
  theme(plot.title.position = "plot")


###### significant in minimal ######

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>%  
  filter(significant_min == TRUE)
  
  k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_min) & grepl(combo_term_count[i,1], kegg$ko)))
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

kterm_sig_z_Pstmin <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sig_z_Pstmin_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_path_count[i,1], kegg_with_path$A)))
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

pathA_sig_z_Pstmin <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sig_z_Pstmin_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
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

pathB_sig_z_Pstmin <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sig_z_Pstmin_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
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

pathC_sig_z_Pstmin <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sig_z_Pstmin_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Includes Significantly Higher and Lower than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")



######## KEGG - sig enrich Pst complete by z score, sig in at least two of three replicates ################
abundance <- read_tsv("processed_data/pst_z_sigenrich_comp.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>% 
  filter(significant_comp == TRUE)

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
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

kterm_sigenrich_z_Pstcomp <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sigenrich_z_Pstcomp_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
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

pathA_sigenrich_z_Pstcomp <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sigenrich_z_Pstcomp_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
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

pathB_sigenrich_z_Pstcomp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sigenrich_z_Pstcomp_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_Pstcomp <- temp %>% 
  add_column(category = NA, .before = "B") %>% 
  replace_na(list(category = "Level B")) %>% 
  rename(path = B)


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
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

pathC_sigenrich_z_Pstcomp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sigenrich_z_Pstcomp_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15) %>% 
  rename(path = C)

kegg_Pstcomp <- add_row(kegg_Pstcomp, temp) %>% 
  replace_na(list(category = "Level C"))

# graph KEGG B and C in one plot
keggBC_sigenrich_z_Pstcomp <- kegg_Pstcomp %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(path, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Pathway Enrichment - Level B and C",
       subtitle = "Combined Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significantly Higher than the median",
       y= "KEGG Pathway")+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(fct_relevel(category, "Level B", "Level C") ~ ., scales = "free_y", space = "free")+
  #xlim(0,550)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich Pst minimal by z score, sig in at least two of three replicates ################
abundance <- read_tsv("processed_data/pst_z_sigenrich_min.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>%  
  filter(significant_min == TRUE)

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_min) & grepl(combo_term_count[i,1], kegg$ko)))
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

kterm_sigenrich_z_Pstmin <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sigenrich_z_Pstmin_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_path_count[i,1], kegg_with_path$A)))
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

pathA_sigenrich_z_Pstmin <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sigenrich_z_Pstmin_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
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

pathB_sigenrich_z_Pstmin <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sigenrich_z_Pstmin_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_Pstmin <- temp %>% 
  add_column(category = NA, .before = "term") %>% 
  replace_na(list(category = "Level C"))


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
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

pathC_sigenrich_z_Pstmin <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sigenrich_z_Pstmin_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_Pstmin <- add_row(kegg_Pstmin, temp) %>% 
  replace_na(list(category = "Level C"))

# graph KEGG B and C in one plot
keggBC_sigenrich_z_Pstmin <- kegg_Pstmin %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(term, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Pathway Enrichment - Level B and C",
       subtitle = "Combined Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significantly Higher than the median",
       y= "KEGG Pathway")+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(fct_relevel(category, "Level B", "Level C") ~ ., scales = "free_y", space = "free")+
  #xlim(0,550)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich and FC higher Pst complete by z score, sig in at least two of three replicates ################
abundance <- read_tsv("processed_data/pst_z_enrich_and_fc_comp.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>% 
  filter(significant_comp == TRUE)

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
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

kterm_fc_sigenrich_z_Pstcomp <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_fc_sigenrich_z_Pstcomp_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
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

pathA_fc_sigenrich_z_Pstcomp <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_fc_sigenrich_z_Pstcomp_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
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

pathB_fc_sigenrich_z_Pstcomp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_fc_sigenrich_z_Pstcomp_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_fcPstcomp <- temp %>% 
  add_column(category = NA, .before = "term") %>% 
  replace_na(list(category = "Level B"))


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_comp))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
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

pathC_fc_sigenrich_z_Pstcomp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_fc_sigenrich_z_Pstcomp_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_fcPstcomp <- add_row(kegg_fcPstcomp, temp) %>% 
  replace_na(list(category = "bp"))

# graph KEGG B and C in one plot
keggBC_fc_sigenrich_z_Pstcomp <- kegg_fcPstcomp %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(term, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Pathway Enrichment - Level B and C",
       subtitle = "Combined Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the median",
       y= "KEGG Pathway")+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(fct_relevel(category, "Level B", "Level C") ~ ., scales = "free_y", space = "free")+
  #xlim(0,550)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich and FC higher Pst minimal by z score, sig in at least two of three replicates ################
abundance <- read_tsv("processed_data/pst_z_enrich_and_fc_min.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", "imp_mednorm_log_pst_complete_1", 
                "imp_mednorm_log_pst_complete_2", "imp_mednorm_log_pst_complete_3",
                "imp_mednorm_log_pst_minimal_1", "imp_mednorm_log_pst_minimal_2",
                "imp_mednorm_log_pst_minimal_3", "pst_comp_sig", "pst_min_sig") %>% 
  mutate(significant_comp = if_else(pst_comp_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_comp = as.logical(significant_comp)) %>% 
  mutate(significant_min = if_else(pst_min_sig >=2, "TRUE", "FALSE")) %>% 
  mutate(significant_min = as.logical(significant_min)) %>% 
  inner_join(total_k_terms, by = "gi_number") %>%  
  filter(significant_min == TRUE)

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]

#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_ko <- double()
number_not_ko <- double()
number_diff_and_ko <- double()
contingency_table <- list()
end_of_loop <- nrow(combo_term_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_ko <- c(number_in_ko, combo_term_count[i,3])
  number_not_ko <- c(number_not_ko, sum(combo_term_count$N.y) - as.numeric(number_in_ko[i]))
  number_diff_and_ko <- c(number_diff_and_ko, sum(grepl(TRUE, kegg$significant_min) & grepl(combo_term_count[i,1], kegg$ko)))
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

kterm_fc_sigenrich_z_Pstmin <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_fc_sigenrich_z_Pstmin_sortbyabund <- k_term_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(ko, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level D, K-terms",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### assign k terms in set of interest to pathways ####
kegg_with_path <- kegg_db %>% 
  filter(ko %in% k_term_nona$ko) %>% 
  inner_join(kegg, kegg_db, by = "ko")


#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_path <- c(number_in_path, combo_path_count[i,3])
  number_not_path <- c(number_not_path, sum(combo_path_count$N.y) - as.numeric(number_in_path[i]))
  number_diff_and_path <- c(number_diff_and_path, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_path_count[i,1], kegg_with_path$A)))
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

pathA_fc_sigenrich_z_Pstmin <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_fc_sigenrich_z_Pstmin_sortbyabund <- k_path_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  #slice(1:10) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(A, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level A",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")


#### evaluate kegg level B ####
# make contingency_tables
k_pathB_count <- kegg_with_path %>% 
  group_by(B) %>% 
  summarize(N=n())

combo_pathB_count <- full_join(k_pathB_count, db_k_pathB_count, by = "B")
combo_pathB_count <- combo_pathB_count[order(combo_pathB_count$B),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i]))
  number_diff_and_pathB <- c(number_diff_and_pathB, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathB_count[i,1], kegg_with_path$B)))
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

pathB_fc_sigenrich_z_Pstmin <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_fc_sigenrich_z_Pstmin_sortbyabund <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(B, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level B",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathB_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_fcPstmin <- temp %>% 
  add_column(category = NA, .before = "term") %>% 
  replace_na(list(category = "Level B"))


#### evaluate kegg level C ####
# make contingency_tables
k_pathC_count <- kegg_with_path %>% 
  group_by(C) %>% 
  summarize(N=n())

combo_pathC_count <- full_join(k_pathC_count, db_k_pathC_count, by = "C")
combo_pathC_count <- combo_pathC_count[order(combo_pathC_count$C),]

diff_abund <- sum(grepl(TRUE, kegg$significant_min))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)#-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i]))
  number_diff_and_pathC <- c(number_diff_and_pathC, sum(grepl(TRUE, kegg_with_path$significant_min) & grepl(combo_pathC_count[i,1], kegg_with_path$C)))
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

pathC_fc_sigenrich_z_Pstmin <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
#  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_fc_sigenrich_z_Pstmin_sortbyabund <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(-actual_count) %>% 
  slice(1:15) %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Enrichment - level C",
       subtitle = "Pst Dataset
       Significant in at least Two of Three Minimal Replicates by Z-score
       Significant and FC Higher than the Median
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

# arrange for one graph
temp <- k_pathC_fisher %>% 
  mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
  rename(actual_count = N.x) %>% 
  arrange(adj_fisher_p) %>% 
  filter(adj_fisher_p <= 0.05) %>% 
  slice(1:15)

kegg_fcPstmin <- add_row(kegg_fcPstmin, temp) %>% 
  replace_na(list(category = "Level C"))

# graph KEGG B and C in one plot
keggBC_fc_sigenrich_z_Pstmin <- kegg_fcPstmin %>% 
  pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
  ggplot(aes(y = reorder(term, -adj_fisher_p), fill = protein_set, x = counts))+
  geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("gray26", "gray88"))+
  labs(title = "KEGG Pathway Enrichment - Level B and C",
       subtitle = "Combined Dataset
       Significant in at least Two of Three Complete Replicates by Z-score
       Significant and FC Higher than the median",
       y= "KEGG Pathway")+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(fct_relevel(category, "Level B", "Level C") ~ ., scales = "free_y", space = "free")+
  #xlim(0,550)+
  theme_classic()+
  theme(plot.title.position = "plot")


################ print figs ################
pdf("figures_kegg_Pst_z_Dataset.pdf")
kterm_sig_z_Pstcomp
pathA_sig_z_Pstcomp
pathB_sig_z_Pstcomp
pathC_sig_z_Pstcomp
kterm_sig_z_Pstmin
pathA_sig_z_Pstmin
pathB_sig_z_Pstmin
pathC_sig_z_Pstmin
kterm_sigenrich_z_Pstcomp
pathA_sigenrich_z_Pstcomp
pathB_sigenrich_z_Pstcomp
pathC_sigenrich_z_Pstcomp
kterm_sigenrich_z_Pstmin
pathA_sigenrich_z_Pstmin
pathB_sigenrich_z_Pstmin
pathC_sigenrich_z_Pstmin
kterm_fc_sigenrich_z_Pstcomp
pathA_fc_sigenrich_z_Pstcomp
pathB_fc_sigenrich_z_Pstcomp
pathC_fc_sigenrich_z_Pstcomp
kterm_fc_sigenrich_z_Pstmin
pathA_fc_sigenrich_z_Pstmin
pathB_fc_sigenrich_z_Pstmin
pathC_fc_sigenrich_z_Pstmin
kterm_sig_z_Pstcomp_sortbyabund
pathA_sig_z_Pstcomp_sortbyabund
pathB_sig_z_Pstcomp_sortbyabund
pathC_sig_z_Pstcomp_sortbyabund
kterm_sig_z_Pstmin_sortbyabund
pathA_sig_z_Pstmin_sortbyabund
pathB_sig_z_Pstmin_sortbyabund
pathC_sig_z_Pstmin_sortbyabund
kterm_sigenrich_z_Pstcomp_sortbyabund
pathA_sigenrich_z_Pstcomp_sortbyabund
pathB_sigenrich_z_Pstcomp_sortbyabund
pathC_sigenrich_z_Pstcomp_sortbyabund
kterm_sigenrich_z_Pstmin_sortbyabund
pathA_sigenrich_z_Pstmin_sortbyabund
pathB_sigenrich_z_Pstmin_sortbyabund
pathC_sigenrich_z_Pstmin_sortbyabund
kterm_fc_sigenrich_z_Pstcomp_sortbyabund
pathA_fc_sigenrich_z_Pstcomp_sortbyabund
pathB_fc_sigenrich_z_Pstcomp_sortbyabund
pathC_fc_sigenrich_z_Pstcomp_sortbyabund
kterm_fc_sigenrich_z_Pstmin_sortbyabund
pathA_fc_sigenrich_z_Pstmin_sortbyabund
pathB_fc_sigenrich_z_Pstmin_sortbyabund
pathC_fc_sigenrich_z_Pstmin_sortbyabund
dev.off()

svglite::svglite("figures/Pst_z_comp_KEGG_B.svg", height = 1.75, width = 5)
pathB_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_comp_KEGG_C.svg", height = 3, width = 5)
pathC_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_min_KEGG_B.svg", height = 1.9, width = 6)
pathB_sigenrich_z_Pstmin
dev.off()

svglite::svglite("figures/Pst_z_min_KEGG_C.svg", height = 3, width = 5)
pathC_sigenrich_z_Pstmin
dev.off()

svglite::svglite("figures/Pst_z_fccomp_KEGG_B.svg", height = 3, width = 5)
pathB_fc_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_fccomp_KEGG_C.svg", height = 3, width = 5)
pathC_fc_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_fcmin_KEGG_B.svg", height = 3, width = 5)
pathB_fc_sigenrich_z_Pstmin
dev.off()

svglite::svglite("figures/Pst_z_fcmin_KEGG_C.svg", height = 3, width = 5)
pathC_fc_sigenrich_z_Pstmin
dev.off()



# these would print the combined graphs, but there's a bug since most don't have significant pathways - or the term is just wrong in the code
svglite::svglite("figures/Pst_z_comp_KEGG.svg", height = 6, width = 6)
keggBC_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_min_KEGG.svg", height = 6, width = 6)
keggBC_sigenrich_z_Pstmin
dev.off()

svglite::svglite("figures/Pst_z_fccomp_KEGG.svg", height = 6, width = 6)
keggBC_fc_sigenrich_z_Pstcomp
dev.off()

svglite::svglite("figures/Pst_z_fcmin_KEGG.svg", height = 6, width = 6)
keggBC_fc_sigenrich_z_Pstmin
dev.off()