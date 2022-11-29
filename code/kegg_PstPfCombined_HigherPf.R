library(tidyverse)
library(dplyr)
library(readxl)

######## needed for all comparisons ########
db_total_k_terms <- read_tsv("raw_data/kegg_combined_DB_hits.txt") %>% 
  dplyr::select("accession", "ko")

db_k_term_count <- db_total_k_terms %>% 
  group_by(ko) %>% 
  summarize(N=n())

total_k_terms <- read_excel("raw_data/kegg_combined_total.xlsx") %>% 
  dplyr::select("gi_number", "accession", "ko") %>% 
  mutate(gi_number = as.double(gi_number)) 

# there's some homology between Pst and Pf that is accounted for in the DB but not in the "combined total"
# innerjoin to eliminate homologous entries that will artifactually inflate k terms
total_k_terms <- inner_join(total_k_terms, db_total_k_terms, by = "accession") %>% 
  select("gi_number", "ko.x") %>% 
  rename(ko = ko.x)

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


######## KEGG - sig species effect, higher in Pf comp and min than in Pst comp and min ################
abundance <- read_tsv("processed_data/higher_Pf_comp_and_min_vs_Pst.tsv")

# assign k terms to proteins in the set of interest
kegg <- abundance %>% 
  dplyr::select("gi_number", "description", contains("imp_mc_log"), 
                "adj_complete_p", "significant_comp", "comp_fc", "log2_comp_fc", "significant_species") %>% 
  inner_join(total_k_terms, by = "gi_number")

k_term_count <- kegg %>% 
  group_by(ko) %>% 
  summarize(N=n())

# make list with counts from db and set of interest
combo_term_count <- full_join(k_term_count, db_k_term_count, by = "ko")
combo_term_count <- combo_term_count[order(combo_term_count$ko),]


#### evaluate lowest level of kegg - the k terms ####
# make contingency_tables
diff_abund <- sum(grepl(TRUE, kegg$significant_species))
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

kterm_highPf <- k_term_fisher %>% 
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
       subtitle = "Combined Dataset
       Significant Species Effect 
       Significantly Higher Pf Complete and Minimal vs Pst Complete and Minimal",
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

write_tsv(kegg_with_path, "processed_data/kegg_higher_Pf_comp_and_min_vs_Pst.tsv")

#### evaluate highest kegg level A ####
# make contingency_tables
k_path_count <- kegg_with_path %>% 
  group_by(A) %>% 
  summarize(N=n())

combo_path_count <- full_join(k_path_count, db_k_path_count, by = "A")
combo_path_count <- combo_path_count[order(combo_path_count$A),]

diff_abund <- sum(grepl(TRUE, kegg$significant_species))
number_in_path <- double()
number_not_path <- double()
number_diff_and_path <- double()
contingency_table_path <- list()
end_of_loop <- nrow(combo_path_count)-1
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

pathA_highPf <- k_path_fisher %>% 
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
       subtitle = "Combined Dataset
       Significant Species Effect 
       Significantly Higher Pf Complete and Minimal vs Pst Complete and Minimal",
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

diff_abund <- sum(grepl(TRUE, kegg$significant_species))
number_in_pathB <- double()
number_not_pathB <- double()
number_diff_and_pathB <- double()
contingency_table_pathB <- list()
end_of_loop <- nrow(combo_pathB_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathB <- c(number_in_pathB, combo_pathB_count[i,3])
  number_not_pathB <- c(number_not_pathB, sum(combo_pathB_count$N.y) - as.numeric(number_in_pathB[i])) # grepl(TRUE, !is.na(db_kegg_with_path$B))
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

pathB_highPf <- k_pathB_fisher %>% 
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
       subtitle = "Combined Dataset
       Significant Species Effect 
       Significantly Higher Pf Complete and Minimal vs Pst Complete and Minimal",
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

diff_abund <- sum(grepl(TRUE, kegg$significant_species))
number_in_pathC <- double()
number_not_pathC <- double()
number_diff_and_pathC <- double()
contingency_table_pathC <- list()
end_of_loop <- nrow(combo_pathC_count)-1
for(i in 1:end_of_loop){ #get rid of NA entry, many don't map to a kegg term
  number_in_pathC <- c(number_in_pathC, combo_pathC_count[i,3])
  number_not_pathC <- c(number_not_pathC, sum(combo_pathC_count$N.y) - as.numeric(number_in_pathC[i])) # grepl(TRUE, !is.na(db_kegg_with_path$C))
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

pathC_highPf <- k_pathC_fisher %>% 
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
       subtitle = "Combined Dataset
       Significant Species Effect
       Significantly Higher Pf Complete and Minimal vs Pst Complete and Minimal",
       x= "Count (Adj. P-value)",
       y= "KEGG Pathway (level C)")+
  guides(fill = guide_legend(reverse = TRUE))+
  #xlim(0,800)+
  theme_classic()+
  theme(plot.title.position = "plot")




# nothing enriched? plot just the terms that are present in this set based on abundance
# pathC_sigint_fcenrichPst_sighighPf_virpred_sortbyabund <- k_pathC_fisher %>% 
#   mutate(expected_count = (N.y/sum(N.y))*sum(N.x, na.rm = TRUE)) %>% 
#   rename(actual_count = N.x) %>% 
#   arrange(-actual_count) %>% 
#   slice(1:15) %>% 
#   pivot_longer(cols = c(actual_count, expected_count), names_to = "protein_set", values_to = "counts") %>% 
#   ggplot(aes(y = reorder(C, -adj_fisher_p), fill = protein_set, x = counts))+
#   geom_text(aes(label = signif(adj_fisher_p, 3)), position = position_dodge(1), hjust = -0.1, size = 3)+
#   geom_col(position = "dodge")+
#   scale_fill_manual(values = c("gray26", "gray88"))+
#   labs(title = "KEGG Enrichment - level C",
#        subtitle = "Combined Dataset
#        Significant Interaction Effect Species x Media
#        Significantly Enriched and FC Higher Pst Minimal vs Pst Complete
#        Significantly Higher Pf Complete vs Pst Complete
#        Predicted Virulence Factor
#        *Arranged by Actual Count",
#        x= "Count (Adj. P-value)",
#        y= "KEGG Pathway (level C)")+
#   guides(fill = guide_legend(reverse = TRUE))+
#   #xlim(0,800)+
#   theme_classic()+
#   theme(plot.title.position = "plot")


################ print figs ################
pdf("figures_kegg_PstPfCombined_HighPf.pdf")
kterm_highPf
pathA_highPf
pathB_highPf
pathC_highPf
dev.off()

svglite::svglite("figures/higher_Pf_vs_Pst_KEGG_B.svg", height = 2, width = 6)
pathB_highPf
dev.off()

svglite::svglite("figures/higher_Pf_vs_Pst_KEGG_C.svg", height = 2, width = 6)
pathC_highPf
dev.off()

