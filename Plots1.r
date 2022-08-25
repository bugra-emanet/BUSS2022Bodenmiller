library(tictoc)
library(ggplot2)
library(SingleCellExperiment)
library(scales)
library(ggpubr)
library(viridis)  

### Read data
sce <- readRDS("/home/bugra/data/SCE/sce_t1d.rds")

colData(sce)$distance_to_islet <- -1 * colData(sce)$distance_to_islet 

colData(sce)$sample_id <- NULL
colData(sce)$image_area <- colData(sce)$image_width * colData(sce)$image_height
colnames(colData(sce))[colnames(colData(sce))=="CellType"] <- "cell_type"
unique(colData(sce)$cell_type)

immune_cell_type <- c("Macrophage_1","Macrophage_2","Neutrophil","T_CD8","T_CD4","B")

islet_cell_type <- c("Delta","Alpha","Beta")

exocrine_cell_type <- c("Acinar","Ductal")

other_cell_type <- c("Endothelial","Other")

a <- colData(sce)$cell_type %in% immune_cell_type
b <- colData(sce)$cell_type %in% islet_cell_type
c <- colData(sce)$cell_type %in% exocrine_cell_type
d <- colData(sce)$cell_type %in% other_cell_type

colData(sce)[a,"cell_category"] <- "immune"
colData(sce)[b,"cell_category"] <- "islet"
colData(sce)[c,"cell_category"] <- "exocrine"
colData(sce)[d,"cell_category"] <- "other"
source("helpers.R", echo=F)


### UMAP
library(scater)
umap_plot <- plotReducedDim(sce,"UMAP_scaled",colour_by = "cell_type") + xlab("UMAP1") + ylab("UMAP2")


### Get Global Density

tic("Global Intensity for immune")
colData(sce)[,"immune_global_intensity"] <- get_global_intensity(sce,
                                                                 cell_type = "immune" ,
                                                                 cell_type_column = "cell_category", 
                                                                 distance_column = "distance_to_islet", 
                                                                 sample_id_column = "image_id")

toc()


### Get global density for immune cell types
tic("Global Intensity for every immune cell type")
for (cell_type in immune_cell_type){
  name <- paste(cell_type,"global_intensity",sep = "_")
  colData(sce)[,name] <-  get_global_intensity(sce,
                                               cell_type = cell_type ,
                                               cell_type_column = "cell_type", 
                                               distance_column = "distance_to_islet", 
                                               sample_id_column = "image_id")
}

 ### Get Data in long format

scores <- as.data.frame(colData(sce)) %>%select("image_id","donor_type","case_id",ends_with("global_intensity"))

scores <- scores[!duplicated(scores$image_id),]

scores[is.na(scores)] <- 0

df_long <- scores %>% pivot_longer(cols = ends_with("global_intensity"),names_to = "score_type",values_to = "score")

t1d_stages <- c("NoDiabetes", "mAAb+", "RecentOnset", "LongDuration")
order_samples <- unique(colData(sce)[, c("donor_type", "case_id")])
order_samples <- order_samples[order(match(order_samples$donor_type, t1d_stages)), ]

df_long <-df_long %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)

immune_long<-  df_long[df_long$score_type == "immune_global_intensity",]

type_long <-  df_long[df_long$score_type != "immune_global_intensity",]


toc()


### Get Local Density

ROI <- 30


tic("Local Intensity Immune")
sce[["immune_local_intensity_30"]] <- get_local_intensity(sce, 
                                                          ROI_size = ROI, 
                                                          cell_type = "immune", 
                                                          cell_type_column = "cell_category", 
                                                          distance_column = "distance_to_islet", 
                                                          sample_id_column = "image_id")
toc()


tic("Local Intensity Every Immune Cell Type")
for (cell_type in immune_cell_type){
  name <- paste(cell_type,"local_intensity_30",sep = "_")
  colData(sce)[,name] <- get_local_intensity(sce, ROI, cell_type = cell_type , 
                                             cell_type_column = "cell_type", 
                                             distance_column = "distance_to_islet", 
                                             sample_id_column = "image_id")
}
toc()


### Get Data in long format
scores_local <- as.data.frame(colData(sce)) %>% select("image_id","donor_type","case_id",ends_with("local_intensity_30"))

scores_local <- scores_local[!duplicated(scores_local$image_id),]

df_long_local <- scores_local %>% pivot_longer(cols = ends_with("local_intensity_30"),names_to = "score_type",values_to = "score")


df_long_local <-df_long_local %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)

immune_long_local<-  df_long_local[df_long_local$score_type == "immune_local_intensity_30",]

type_long_local <-  df_long_local[df_long_local$score_type != "immune_local_intensity_30",]


immune_data_long <- rbind(immune_long,immune_long_local)



### Get distance score immune

tic()
colData(sce)[["immune_distance_score_R"]] <- get_distance_score_reciprocal(sce, cell_type = "immune",cell_type_column = "cell_category",distance_column = "distance_to_islet",sample_id_column = "image_id")
toc()

### Get distance score for every cell type
tic()
for (cell_type in immune_cell_type){
  name_R <- paste(cell_type,"distance_score_R",sep = "_")
  colData(sce)[,name_R] <- get_distance_score_reciprocal(sce, cell_type = cell_type , 
                                                         cell_type_column = "cell_type", 
                                                         distance_column = "distance_to_islet", 
                                                         sample_id_column = "image_id")
  
  
}
toc()



scores_distance <- as.data.frame(colData(sce)) %>% select("image_id","donor_type","case_id",ends_with("distance_score_R"))

scores_distance <- scores_distance[!duplicated(scores_distance$image_id),]

scores_distance[is.na(scores_distance)] <- 0


df_long_distance <- scores_distance %>% pivot_longer(cols = ends_with("distance_score_R"),names_to = "score_type",values_to = "score")

df_long_distance<-df_long_distance %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)

immune_long_distance <-  df_long_distance[df_long_distance$score_type == "immune_distance_score_R",]

type_long_distance <-  df_long_distance[df_long_distance$score_type != "immune_distance_score_R",]

### Get plots
my_labeller<- function(x){
  x <- NULL
  return(data.frame(labels = c("B","Macrophage 1","Macrophage 2","Neutrophil","T CD4","T CD8")))
}


### Get mean global & local densities per patient


scores <-scores %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)


mean_prop<- scores%>%
  group_by(case_id,donor_type) %>%
  summarise_at(vars(colnames(scores)[5:10]), list(name = mean))


mean_prop <- mean_prop %>% pivot_longer(cols = ends_with("intensity_name"),names_to = "score_type",values_to = "score")
my_theme_2  <- theme(panel.grid.major.x = element_blank(),
                     axis.title.x = element_text(size = 10),
                     axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.text.y = element_text(size = 6),
                     axis.title.y = element_text(size = 10),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 8),
                     legend.key.size = unit(0.3,"cm"),
                     strip.text.x = element_text(size = 10),
                     axis.ticks.length=unit(.05, "cm"))


scores_local <-scores_local %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)


mean_prop_local<- scores_local%>%
  group_by(case_id,donor_type) %>%
  summarise_at(vars(colnames(scores_local)[5:10]), list(name = mean))


mean_prop_local <- mean_prop_local %>% pivot_longer(cols = ends_with("intensity_30_name"),names_to = "score_type",values_to = "score")








my_theme <- theme(
  axis.text =element_text(size = 6) ,
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 10),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.3,"cm"),
  strip.text.x = element_text(size = 10),
  axis.ticks.length=unit(.05, "cm"))

my_theme_bar <- theme(panel.grid.major.x = element_blank(),
                                      axis.title.x = element_text(size = 10),
                                      axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
                                      axis.text.y = element_text(size = 6),
                                      axis.title.y = element_text(size = 10),
                                      legend.title = element_text(size = 8),
                                      legend.text = element_text(size = 8),
                                      legend.key.size = unit(0.3,"cm"),
                                      strip.text.x = element_text(size = 10),
                                      axis.ticks.length=unit(.05, "cm"))

point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)
my_comparisons = list( c("mAAb+", "RecentOnset"),c("NoDiabetes", "RecentOnset"),c("LongDuration", "RecentOnset") )


p1 <- ggplot(data = immune_data_long, mapping = aes(y = score, x = donor_type, col = donor_type)) + geom_violin() + facet_wrap(~score_type,scales = "fixed",labeller =  function(x){
  x <- NULL
  return(data.frame(labels = c("Immune Global", "Immune Local")))
}) + my_theme + stat_summary(fun = "mean", mapping = aes(shape="mean"), geom = "point",color = "red",size = 1.2) + geom_jitter(size = 0.4)+scale_shape_manual(values=c(16)) + labs(shape = NULL, x = "Case ID", y = "Relative Density", col = "Disease Stage") +scale_y_continuous(labels = point) +
  stat_compare_means(comparisons = my_comparisons, show.legend = F,label.y = c(0.20,0.15,0.10),paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),hide.ns = T))+
  stat_compare_means(label.y = 0.185,show.legend = F, size = 3)

p2 <- ggplot(type_long,aes(y= score,x = donor_type,col=donor_type)) + geom_violin() + geom_jitter(size = 0.2)+ facet_wrap(~score_type,labeller = my_labeller,scale = "fixed") + my_theme  + stat_summary(fun = "mean", mapping = aes(shape="mean"), geom = "point",color = "red",size = 1.2) + scale_shape_manual(values=c(16))+ labs(shape = NULL, x = "Case ID", y = "Global Relative Density", col = "Disease Stage") + scale_y_continuous(labels = point) +
  stat_compare_means(comparisons = my_comparisons, show.legend = F,label.y = c(0.1,0.05,0.030),paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),hide.ns = T))+
  stat_compare_means(show.legend = F, size = 2)

p3 <- ggplot(type_long_local,aes(y= score,x = donor_type,col= donor_type)) + geom_violin( )  + geom_jitter(size = 0.2)+ facet_wrap(~score_type,scales = "fixed",labeller = my_labeller) + my_theme +  stat_summary(fun = "mean", mapping = aes(shape="mean"), geom = "point",color = "red",size = 1.2) + scale_shape_manual(values=c(16)) + labs(shape = NULL, x = "Case ID", y = "Local Relative Density", col = "Disease Stage") + scale_y_continuous(labels = point)+
  stat_compare_means(comparisons = my_comparisons, show.legend = F,label.y = c(0.1,0.05,0.03),paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),hide.ns = T))+
  stat_compare_means(show.legend = F, size = 2)
      
p4 <- ggplot(type_long_distance,aes(y= score,x = donor_type,col= donor_type)) + geom_violin() + geom_jitter(size = 0.2)+ facet_wrap(~score_type,scales = "fixed",labeller = my_labeller) + my_theme +  stat_summary(fun = "mean", mapping = aes(shape="mean"), geom = "point",color = "red",size = 1.2) + scale_shape_manual(values=c(16)) + labs(shape = NULL, x = "Case ID", y = "Distance Score", col = "Disease Stage") + scale_y_continuous(labels = point)   +stat_compare_means(label.y = 1,show.legend = F, size =2)+
  stat_compare_means(comparisons = my_comparisons, show.legend = F,hide.ns = F,label.y = c(0.40,0.25,0.10),paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),hide.ns = T))

p5 <- ggplot(immune_long_distance,aes(y= score,x = donor_type,col= donor_type)) + geom_violin()+ geom_jitter(size = 0.5) + facet_wrap(~score_type,scales = "free_y",labeller =  function(x){
  x <- NULL
  return(data.frame(labels = c("Immune Distance Score")))}) + my_theme +  stat_summary(fun = "mean", mapping = aes(shape="mean"), geom = "point",color = "red",size = 1.2) + scale_shape_manual(values=c(16)) + labs(shape = NULL, x = "Case ID", y = "Distance Score", col = "Disease Stage") + scale_y_continuous(labels = point)+
  stat_compare_means(comparisons = my_comparisons, show.legend = F,label.y = c(0.20,0.15,0.10),paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"),hide.ns = T))+
  stat_compare_means(label.y = 0.55,show.legend = F, size = 3)

bar_chart_local <-ggplot(data = mean_prop_local,mapping = aes(x= case_id,y = score,fill = score_type)) + geom_bar(stat= "identity") + my_theme_bar +
  xlab("Case ID") + ylab("Local Relative Density") + labs(fill = "Cell Type") + geom_vline(xintercept = seq(0.5, nrow(mean_prop), by = 5), color="red", size=1,, alpha=0.2,linetype = "dashed")+
  annotate(geom="text", x=3, y=0.14, label="NoDiabetes",
           color="black")+annotate(geom="text", x=8, y=0.14, label="mAAb+",
                                   color="black")+annotate(geom="text", x=13, y=0.14, label="recentOnset",
                                                           color="black")+annotate(geom="text", x=18, y=0.14, label="LongDuration",
                                                                                   color="black")+ scale_fill_manual(labels = c("B","Macrophage 1","Macrophage 2","Neutrophil","T CD4","T CD8"),values = viridis(6))



bar_chart_global<- ggplot(data = mean_prop,mapping = aes(x= case_id,y = score,fill = score_type)) + geom_bar(stat= "identity") + my_theme_bar +
  xlab("Case ID") + ylab("Global Relative Density") + labs(fill = "Cell Type") + geom_vline(xintercept = seq(0.5, nrow(mean_prop), by = 5), color="red", size=1,, alpha=0.2,linetype = "dashed")+
  annotate(geom="text", x=3, y=0.14, label="NoDiabetes",
           color="black")+annotate(geom="text", x=8, y=0.14, label="mAAb+",
                                   color="black")+annotate(geom="text", x=13, y=0.14, label="recentOnset",
                                                           color="black")+annotate(geom="text", x=18, y=0.14, label="LongDuration",
                                                                                   color="black")+ scale_fill_manual(labels = c("B","Macrophage 1","Macrophage 2","Neutrophil","T CD4","T CD8"),values = viridis(6))
### Save plots

ggsave(filename ="immune_densities.pdf",plot = p1,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="global_densities.pdf", plot = p2,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="local_densities.pdf", plot = p3,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="distance_type.pdf", plot = p4,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="distance_immune.pdf", plot = p5,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="bar_chart_global.pdf", plot = bar_chart_global,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="bar_chart_local.pdf", plot = bar_chart_local,device= "pdf", width = 18,height = 12, units = "cm")

ggsave(filename ="umap_plot.pdf", plot = umap_plot, device= "pdf", width = 18,height = 12, units = "cm")


