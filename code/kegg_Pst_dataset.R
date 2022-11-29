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


######## KEGG - sig enrich Pst complete by t test ################
abundance <- read_tsv("processed_data/pst_sigenrich_comp.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_p_value_comp_vs_min", "significant_comp_vs_min", "log2_fc_min_vs_comp") %>%
  rename(significant_comp = significant_comp_vs_min) %>% 
  inner_join(total_k_terms, by = "gi_number")

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

kterm_sigenrich_Pstcomp <- k_term_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sigenrich_Pstcomp_sortbyabund <- k_term_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media
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

pathA_sigenrich_Pstcomp <- k_path_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sigenrich_Pstcomp_sortbyabund <- k_path_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media
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

pathB_sigenrich_Pstcomp <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sigenrich_Pstcomp_sortbyabund <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media
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

pathC_sigenrich_Pstcomp <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sigenrich_Pstcomp_sortbyabund <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Complete Media
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich Pst minimal by t test ################
abundance <- read_tsv("processed_data/pst_sigenrich_min.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_p_value_comp_vs_min", "significant_comp_vs_min", "log2_fc_min_vs_comp") %>%
  rename(significant_comp = significant_comp_vs_min) %>% 
  inner_join(total_k_terms, by = "gi_number")

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

kterm_sigenrich_Pstmin <- k_term_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_sigenrich_Pstmin_sortbyabund <- k_term_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media
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

pathA_sigenrich_Pstmin <- k_path_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_sigenrich_Pstmin_sortbyabund <- k_path_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media
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

pathB_sigenrich_Pstmin <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_sigenrich_Pstmin_sortbyabund <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media
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

pathC_sigenrich_Pstmin <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_sigenrich_Pstmin_sortbyabund <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significantly Higher in Minimal Media
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich and FC higher Pst complete by t test ################
abundance <- read_tsv("processed_data/pst_enrich_and_fc_comp.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_p_value_comp_vs_min", "significant_comp_vs_min", "log2_fc_min_vs_comp") %>% 
  rename(significant_comp = significant_comp_vs_min) %>% 
  inner_join(total_k_terms, by = "gi_number")

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

kterm_fc_sigenrich_Pstcomp <- k_term_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_fc_sigenrich_Pstcomp_sortbyabund <- k_term_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media
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

pathA_fc_sigenrich_Pstcomp <- k_path_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_fc_sigenrich_Pstcomp_sortbyabund <- k_path_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media
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

pathB_fc_sigenrich_Pstcomp <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_fc_sigenrich_Pstcomp_sortbyabund <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media
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

pathC_fc_sigenrich_Pstcomp <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_fc_sigenrich_Pstcomp_sortbyabund <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Complete Media
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")


######## KEGG - sig enrich and FC higher Pst minimal by t test ################
abundance <- read_tsv("processed_data/pst_enrich_and_fc_min.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_p_value_comp_vs_min", "significant_comp_vs_min", "log2_fc_min_vs_comp") %>% 
  rename(significant_comp = significant_comp_vs_min) %>% 
  inner_join(total_k_terms, by = "gi_number")

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

kterm_fc_sigenrich_Pstmin <- k_term_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Term (level D)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,300)+
  theme_classic()+
  theme(plot.title.position = "plot")

kterm_fc_sigenrich_Pstmin_sortbyabund <- k_term_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media
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

pathA_fc_sigenrich_Pstmin <- k_path_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level A)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,2000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathA_fc_sigenrich_Pstmin_sortbyabund <- k_path_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media
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

pathB_fc_sigenrich_Pstmin <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level B)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,1000)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathB_fc_sigenrich_Pstmin_sortbyabund <- k_pathB_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media
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

pathC_fc_sigenrich_Pstmin <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")

pathC_fc_sigenrich_Pstmin_sortbyabund <- k_pathC_fisher %>% 
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
       Significant by t-test
       Significant and FC Higher in Minimal Media
       *Arranged by Actual Count",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")


################ print figs ################
pdf("figures_kegg_Pst_Dataset.pdf")
# kterm_sig_Pstcomp # have not yet done analysis with all significant, including both higher and lower values
# pathA_sig_Pstcomp
# pathB_sig_Pstcomp
# pathC_sig_Pstcomp
# kterm_sig_Pstmin
# pathA_sig_Pstmin
# pathB_sig_Pstmin
# pathC_sig_Pstmin
kterm_sigenrich_Pstcomp
pathA_sigenrich_Pstcomp
pathB_sigenrich_Pstcomp
pathC_sigenrich_Pstcomp
kterm_sigenrich_Pstmin
pathA_sigenrich_Pstmin
pathB_sigenrich_Pstmin
pathC_sigenrich_Pstmin
# kterm_fc_sigenrich_Pstcomp # there's only one protein in this category and it doesn't map to a k-term
# pathA_fc_sigenrich_Pstcomp
# pathB_fc_sigenrich_Pstcomp
# pathC_fc_sigenrich_Pstcomp
kterm_fc_sigenrich_Pstmin
pathA_fc_sigenrich_Pstmin
pathB_fc_sigenrich_Pstmin
pathC_fc_sigenrich_Pstmin
# kterm_sig_Pstcomp_sortbyabund # these would be the sorted ones for all significance grouping
# pathA_sig_Pstcomp_sortbyabund
# pathB_sig_Pstcomp_sortbyabund
# pathC_sig_Pstcomp_sortbyabund
# kterm_sig_Pstmin_sortbyabund
# pathA_sig_Pstmin_sortbyabund
# pathB_sig_Pstmin_sortbyabund
# pathC_sig_Pstmin_sortbyabund
kterm_sigenrich_Pstcomp_sortbyabund
pathA_sigenrich_Pstcomp_sortbyabund
pathB_sigenrich_Pstcomp_sortbyabund
pathC_sigenrich_Pstcomp_sortbyabund
kterm_sigenrich_Pstmin_sortbyabund
pathA_sigenrich_Pstmin_sortbyabund
pathB_sigenrich_Pstmin_sortbyabund
pathC_sigenrich_Pstmin_sortbyabund
# kterm_fc_sigenrich_Pstcomp_sortbyabund # there's only one protein in this category and it doesn't map to a k-term
# pathA_fc_sigenrich_Pstcomp_sortbyabund
# pathB_fc_sigenrich_Pstcomp_sortbyabund
# pathC_fc_sigenrich_Pstcomp_sortbyabund
kterm_fc_sigenrich_Pstmin_sortbyabund
pathA_fc_sigenrich_Pstmin_sortbyabund
pathB_fc_sigenrich_Pstmin_sortbyabund
pathC_fc_sigenrich_Pstmin_sortbyabund
dev.off()

svglite::svglite("figures/Pst_comp_KEGG_B.svg", height = 3, width = 5)
pathB_sigenrich_Pstcomp
dev.off()

svglite::svglite("figures/Pst_comp_KEGG_C.svg", height = 3, width = 5)
pathC_sigenrich_Pstcomp
dev.off()

svglite::svglite("figures/Pst_min_KEGG_B.svg", height = 2.5, width = 6)
pathB_sigenrich_Pstmin
dev.off()

svglite::svglite("figures/Pst_min_KEGG_C.svg", height = 5, width = 7)
pathC_sigenrich_Pstmin
dev.off()

svglite::svglite("figures/Pst_fccomp_KEGG_B.svg", height = 3, width = 5)
pathB_fc_sigenrich_Pstcomp
dev.off()

svglite::svglite("figures/Pst_fccomp_KEGG_C.svg", height = 3, width = 5)
pathC_fc_sigenrich_Pstcomp
dev.off()

svglite::svglite("figures/Pst_fcmin_KEGG_B.svg", height = 3, width = 5)
pathB_fc_sigenrich_Pstmin
dev.off()

svglite::svglite("figures/Pst_fcmin_KEGG_C.svg", height = 3, width = 5)
pathC_fc_sigenrich_Pstmin
dev.off()