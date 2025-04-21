library(SCopeLoomR)
library(SCENIC)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggplot2)

loom <- open_loom('SCENIC+_gene_based_filter.loom')
cells_AUC_gene <- AUCell::getAUC(get_regulons_AUC(loom, column.attr.name='RegulonsAUC'))
rownames(cells_AUC_gene) <- gsub('_extended', '', rownames(cells_AUC_gene) )
cell_data <- get_cell_annotation(loom)
dgem <- get_dgem(loom)

All_selected_eRegulon = scan('selected_eRegulons_All_0.6.txt', what="", sep="\n")
Myeloid_selected_eRegulon = scan('selected_eRegulons_Myeloid_0.8.txt', what="", sep="\n")
B_selected_eRegulon = scan('selected_eRegulons_B_0.8.txt', what="", sep="\n")
NK_selected_eRegulon = scan('selected_eRegulons_NK_0.8.txt', what="", sep="\n")
CD8T_selected_eRegulon = scan('selected_eRegulons_CD8T_0.8.txt', what="", sep="\n")
CD4T_selected_eRegulon = scan('selected_eRegulons_CD4T_0.8.txt', what="", sep="\n")

selected_eRegulon_names = unique(c(All_selected_eRegulon,Myeloid_selected_eRegulon,B_selected_eRegulon,NK_selected_eRegulon,CD8T_selected_eRegulon,CD4T_selected_eRegulon))

regulons_gene <- get_regulons(loom, column.attr.name='Regulons')
regulon_gene_list <- list()
for (row in rownames(regulons_gene)){
  regulon_gene_list[[row]] <- colnames(regulons_gene)[which(regulons_gene[row,] == 1)]
}
regulon_gene_list <- regulon_gene_list[selected_eRegulon_names]
names(regulon_gene_list) <- gsub('_extended', '', names(regulon_gene_list))

regulon_gene_list <- tapply(unlist(regulon_gene_list , use.names = FALSE), rep(names(regulon_gene_list), lengths(regulon_gene_list)), FUN = c)

positive <- names(regulon_gene_list)[grep('_+', names(regulon_gene_list), fixed=TRUE)]
negative <- names(regulon_gene_list)[grep('_-', names(regulon_gene_list), fixed=TRUE)]


rss_values <- calcRSS(cells_AUC_gene, cell_data$GEX_final_annotation)
rss_values <- sweep(rss_values,2,colSums(rss_values),`/`)*100
rss_values <- rss_values[,sort(colnames(rss_values))]


cells_AUC_gene_df <- as.data.frame(cells_AUC_gene)
cells_AUC_gene_df$Regulon <- rownames(cells_AUC_gene)

df_long <- cells_AUC_gene_df %>%
  pivot_longer(
    cols = -Regulon,
    names_to = "Cell",
    values_to = "AUC"
  )

df_long <- df_long %>%
  mutate(Cell_type = str_extract(Cell, "(?<=_)[^_]+(?:_[^_]+)*?(?=_[^_]+$)"))

df_long <- df_long %>%
  select(Regulon, Cell_type, AUC)
long_cells_AUC_gene_df = df_long

mean_AUC_df <- long_cells_AUC_gene_df %>%
  group_by(Regulon, Cell_type) %>%
  summarise(mean_AUC = mean(AUC, na.rm = TRUE))

mean_AUC_mat <- mean_AUC_df %>%
  pivot_wider(names_from = Cell_type, values_from = mean_AUC) %>%
  column_to_rownames("Regulon") %>%
  as.matrix()


# Prepare data for plotting
expression_list <- list()
for (x in unique(cell_data$GEX_final_annotation)){
  print(x)
  expression_list[[x]] <- as.data.frame(t(log(rowSums(dgem[,grep(x, cell_data$GEX_final_annotation, fixed = TRUE)])/sum(rowSums(dgem[,grep(x, cell_data$GEX_final_annotation, fixed = TRUE)]))*10^6+1)))
}
exp_mat <- t(data.table::rbindlist(expression_list))
colnames(exp_mat) <- names(expression_list)
exp_mat <- exp_mat



line_order <- c(
              'Bn_TCL1A','Bn_IFIT3','Transitional_B_SOX4','Unswitched_Bm_CD1C','pre-Switched_Bm_JAM3','Atypical_Bm_ITGAX',
                'Switched_Bm_IGHDneg','Switched_Bm_IGHE','Switched_activated_Bm_CD86','Plasma_IGHA1','Plasma_IGHG1','Plasmablast_MKI67',
                
                'pDC_IRF4','AS_DC','ncMono_C1QA','ncMono_FCGR3A','ncMono_IFIT1','intMono_GFRA2','cMono_CD14','cMono_CXCL10','cMono_IFI44L','cMono_IL1B',
                'cDC1_BATF3','cDC2_CD1C','cDC_CSF2RA',
                'HSPC_CD34',
                'ILC2_IL2RA','NK_bright_XCL1','Transitional_NK_GZMK','Mature_NK_dim_FCGR3A','Terminal_NK_dim_CD160neg','Cycling_NK_MKI67',
                
                'MAIT_SLC4A10','gdT1_TRDV1','gdT2_GZMH','gdT2_GZMK','gdT2_IL12RB2','Cycling_T_MKI67','NKT_NCR1',
                'CD8_Tn_CCR7','CD8_Tcm_IFI44L','CD8_Tem_CCR7neg','CD8_CTL_GZMB','CD8_CTL_GZMK','CD8_CTL_IFI44L' ,
                'CD4_Tn_CCR7','CD4_Tcm_CXCR5','CD4_Tcm_IFI44L','CD4_Tem_CCR5','CD4_Tem_CCR7neg','CD4_Tfh-like_CXCR5',
                'CD4_Th1-like_GZMK','CD4_Th17-like_RORC','CD4_Th22-like_CCR10','CD4_Th_CCR4','CD4_Th_LMNA','CD4_Th_TNFRSF11A','CD4_Tr1-like_IL10',
                'CD4_Treg_FCRL3','CD4_Treg_FOXP3','CD4_CTL_GZMH'             
               )

sel_rel <- c(positive, negative)
order_list <- list()
for (i in 1:ncol(rss_values)){
  order_list[[i]] <- vector()
}
for (x in sel_rel){
  i <- which.max(rss_values[x,line_order])
  order_list[[i]] <- c(order_list[[i]], x)
}
sel_rel <- unlist(order_list)

rss_values <- t(apply(rss_values, 1, function(x)(x-min(x))/(max(x)-min(x))))
rel_list <- sel_rel
rel_data <- data.frame()
for (rel in rel_list){
  for (name in colnames(rss_values)){
      tf <- strsplit(rel, split = "_")[[1]][1]
      row_to_add <- t(c(rel, name, exp_mat[tf, name], rss_values[rel,name],mean_AUC_mat[rel,name]))
      rel_data <- rbind(rel_data, row_to_add)
   }
}

colnames(rel_data) <- c('Regulon', 'Cell_type', 'Expression', 'RSS','mean_AUC')
sel_rel_color <- rep('grey',length(sel_rel))
sel_rel_color[which(sel_rel %in% positive)] <- 'forestgreen'
sel_rel_color[which(sel_rel %in% negative)] <- 'red'


rel_data$Regulon <- factor(rel_data$Regulon, levels=sel_rel)
rel_data$Cell_type <- factor(rel_data$Cell_type, levels= line_order)
rel_data$RSS <- as.numeric(rel_data$RSS)
rel_data$mean_AUC <- as.numeric(rel_data$mean_AUC)
rel_data$AUC_scale <- ave(rel_data$mean_AUC, rel_data$Regulon, FUN = scale)
rel_data$AUC_scale[which(rel_data$AUC_scale < -2.5)] <- -2.5
rel_data$AUC_scale[which(rel_data$AUC_scale > 2.5)] <- 2.5
rel_data$Expression <- as.numeric(rel_data$Expression)
rel_data$Expression_scale <- ave(rel_data$Expression, rel_data$Regulon, FUN = scale)
rel_data$Expression_scale[which(rel_data$Expression_scale < -2.5)] <- -2.5
rel_data$Expression_scale[which(rel_data$Expression_scale > 2.5)] <- 2.5
rel_data$Repressor_activator <- unlist(lapply(grepl("+", rel_data$Regulon, fixed=TRUE), function(x) if (x) 'Activators' else 'Repressors'))


# Custom function to draw split tiles
split_tile <- function(x, y, fill1, fill2) {
  data.frame(
    x = c(x - 0.5, x + 0.5, x + 0.5, NA, x + 0.5, x - 0.5, x - 0.5, NA),
    y = c(y - 0.5, y - 0.5, y + 0.5, NA, y + 0.5, y + 0.5, y - 0.5, NA),
    fill = c(fill1, fill1, fill1, NA, fill2, fill2, fill2, NA),
    part = c(rep(1, 4), rep(2, 4))
  )
}

rel_data_Activators = rel_data[rel_data$Repressor_activator == 'Activators', ]

# Store the original factor levels
cell_type_levels <- levels(factor(rel_data_Activators$Cell_type))
regulon_levels <- levels(factor(rel_data_Activators$Regulon))

# Convert Cell_type and Regulon to numeric for plotting
rel_data_Activators <- rel_data_Activators %>%
  mutate(Cell_type_numeric = as.numeric(factor(Cell_type, levels = cell_type_levels)),
         Regulon_numeric = as.numeric(factor(Regulon, levels = regulon_levels)))

# Applying the function to the data
split_tiles <- do.call(rbind, lapply(1:nrow(rel_data_Activators), function(i) {
  split_tile(rel_data_Activators$Cell_type_numeric[i], rel_data_Activators$Regulon_numeric[i], 
             rel_data_Activators$AUC_scale[i], rel_data_Activators$Expression_scale[i])
}))

# Merging the split tile data with the original data
split_tiles$Cell_type <- rel_data_Activators$Cell_type[rep(1:nrow(rel_data_Activators), each = 8)]
split_tiles$Regulon <- rel_data_Activators$Regulon[rep(1:nrow(rel_data_Activators), each = 8)]

# Reverse the color palette
reversed_palette <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

# Plotting
g <- ggplot() +
  geom_polygon(data = split_tiles, aes(x = x, y = y, group = interaction(Cell_type, Regulon, part), fill = fill)) +
  geom_point(data = subset(rel_data_Activators, Top_RSS == TRUE), 
             aes(x = Cell_type_numeric, y = Regulon_numeric, size = RSS), color = 'white', shape = 18) +
  scale_fill_gradientn(colors = reversed_palette, limits = c(-2.5, 2.5)) +
  scale_radius(range = c(0.8, 3), limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:length(cell_type_levels), labels = cell_type_levels) +
  scale_y_continuous(breaks = 1:length(regulon_levels), labels = regulon_levels) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

ggsave('R_Dotplot_RSS_hm_activators_0.6All&0.8Sub_AUC_Expression_top15_Max2.5.pdf', g, width = 35, height = 10)
g


rel_data_Repressors = rel_data[rel_data$Repressor_activator == 'Repressors', ]
rel_data_Repressors <- rel_data_Repressors %>%
  group_by(Cell_type) %>%
  mutate(Rank = rank(-RSS, ties.method = "first")) %>%
  mutate(Top_RSS = ifelse(Rank <= 10, TRUE, FALSE)) %>%
  mutate(Shape = ifelse(Top_RSS, '19', '17')) %>%
  ungroup()

# Store the original factor levels
cell_type_levels <- levels(factor(rel_data_Repressors$Cell_type))
regulon_levels <- levels(factor(rel_data_Repressors$Regulon))

# Convert Cell_type and Regulon to numeric for plotting
rel_data_Repressors <- rel_data_Repressors %>%
  mutate(Cell_type_numeric = as.numeric(factor(Cell_type, levels = cell_type_levels)),
         Regulon_numeric = as.numeric(factor(Regulon, levels = regulon_levels)))

# Applying the function to the data
split_tiles <- do.call(rbind, lapply(1:nrow(rel_data_Repressors), function(i) {
  split_tile(rel_data_Repressors$Cell_type_numeric[i], rel_data_Repressors$Regulon_numeric[i], 
             rel_data_Repressors$AUC_scale[i], rel_data_Repressors$Expression_scale[i])
}))

# Merging the split tile data with the original data
split_tiles$Cell_type <- rel_data_Repressors$Cell_type[rep(1:nrow(rel_data_Repressors), each = 8)]
split_tiles$Regulon <- rel_data_Repressors$Regulon[rep(1:nrow(rel_data_Repressors), each = 8)]

# Reverse the color palette
reversed_palette <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

# Plotting
g <- ggplot() +
  geom_polygon(data = split_tiles, aes(x = x, y = y, group = interaction(Cell_type, Regulon, part), fill = fill)) +
  #geom_point(data = subset(rel_data_Repressors, Top_RSS == TRUE), 
  #           aes(x = Cell_type_numeric, y = Regulon_numeric, size = RSS), color = 'white', shape = 18) +
  scale_fill_gradientn(colors = reversed_palette, limits = c(-2.5, 2.5)) +
  scale_radius(range = c(0.8, 3), limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:length(cell_type_levels), labels = cell_type_levels) +
  scale_y_continuous(breaks = 1:length(regulon_levels), labels = regulon_levels) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

ggsave('R_Dotplot_RSS_hm_Repressor_0.6All&0.8Sub_AUC_Expression_top15_Max2.5.pdf', g, width = 8, height = 10)
g


rel_data_AUC = rel_data[rel_data$AUC_scale >= 1, ]

celltype_regulon_count <- rel_data_AUC %>%
  group_by(Cell_type,Repressor_activator) %>%
  summarise(Regulon_count = n_distinct(Regulon),
  Regulon_list = paste(unique(Regulon), collapse = ", "),
           .groups = 'drop') 

write.csv(celltype_regulon_count,'./eRegulon_0.6All&0.8Sub_count_AUC1.csv')

celltype_regulon_count$Cell_type <- factor(celltype_regulon_count$Cell_type, levels = line_order)
repressor_activator_order <- c("Repressors", "Activators")
celltype_regulon_count$Repressor_activator <- factor(celltype_regulon_count$Repressor_activator, levels = repressor_activator_order)
p = ggplot(celltype_regulon_count, aes(y = Cell_type, x = Regulon_count, fill = Repressor_activator)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Regulon Count by Cell Type and Repressor/Activator",
       y = "Cell Type",
       x = "Regulon Count",
       fill = "Repressor/Activator") +
  theme_minimal()
p
ggsave('eRegulon_0.6All&0.8Sub_count_AUC1.pdf', p, width=7, height=10)