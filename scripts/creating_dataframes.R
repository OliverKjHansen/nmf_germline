library("callr")
library("tidyverse")


args<-commandArgs(trailingOnly = TRUE) 


background_counts <- args[1]
insertions_counts  <- args[2]
deletions_counts <- args[3]

path_to_output_ins <- args[4]
# path_to_input <- args[4]
path_to_output_del <- args[5]

path_to_output_merged <- args[6]



background_matrix <- read.delim(background_counts, sep="\t", header = FALSE)
colnames(background_matrix) <- c("region","kmer","counts")
insertions_matrix <- read.delim(insertions_counts, sep="\t", header = FALSE)
colnames(insertions_matrix) <- c("region","kmer","counts")
deletions_matrix<- read.delim(deletions_counts,  sep="\t", header = FALSE)
colnames(deletions_matrix) <- c("region","kmer","counts")

ms_bck_matrix <- pivot_wider(background_matrix, names_from = kmer, values_from = counts, values_fill = 0, names_sort = TRUE)
ms_bck_matrix <- ms_bck_matrix[str_order(ms_bck_matrix$region, numeric = TRUE),]

ms_ins_matrix <- pivot_wider(insertions_matrix, names_from = kmer, values_from = counts, values_fill = 0, names_sort = TRUE)
ms_ins_matrix <- ms_ins_matrix[str_order(ms_ins_matrix$region, numeric = TRUE),]

ms_del_matrix <- pivot_wider(deletions_matrix, names_from = kmer, values_from = counts, values_fill = 0, names_sort = TRUE)
ms_del_matrix <- ms_del_matrix[str_order(ms_del_matrix$region, numeric = TRUE),]

ms_ins_matrix <- ms_ins_matrix %>% remove_rownames %>% column_to_rownames(var="region")
ms_del_matrix <- ms_del_matrix %>% remove_rownames %>% column_to_rownames(var="region")
ms_bck_matrix <- ms_bck_matrix %>% remove_rownames %>% column_to_rownames(var="region")

ms_ins_df <- list("count_matrix" = ms_ins_matrix,
             "background_matrix" = ms_bck_matrix)
  

ms_del_df <- list("count_matrix" = ms_del_matrix,
             "background_matrix" = ms_bck_matrix)

ms_merged_matrix <- merge(x= ms_ins_matrix, y = ms_del_matrix by = 0, all = TRUE, suffixes = c("_ins","_del"))
ms_merged_matrix <- ms_merged_matrix[str_order(ms_merged_matrix$Row.names, numeric = TRUE),] %>% remove_rownames %>% column_to_rownames(var="Row.names")

ms_bck_merged_matrix <- merge(x= ms_bck_matrix, y = ms_bck_matrix, by = 0, all = TRUE, suffixes = c("_ins","_del"))
ms_bck_merged_matrix <- ms_bck_merged_matrix[str_order(ms_bck_merged_matrix$Row.names, numeric = TRUE),] %>% remove_rownames %>% column_to_rownames(var="Row.names")

ms_merged_df <- list("count_matrix" = ms_merged_matrix,
             "background_matrix" = ms_bck_merged_matrix)

saveRDS(ms_del_df, path_to_output_del)
saveRDS(ms_ins_df, path_to_output_ins)
saveRDS(ms_merged_df, path_to_output_merged)