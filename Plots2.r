library(tictoc)
library(tidyverse)
library(spatstat)
library(scales)
library(SingleCellExperiment)
library(ggplot2)
library(ggpubr)


### Get data and helper functions
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

data <- hyperframe()
## Get patterns

tic("Get patterns for immune")
patterns <- as.solist(lapply(unique(colData(sce)$image_id),function(sample_id){get_ppp_for_image(sce,
                                                                                                 sample_id = sample_id, 
                                                                                                 sample_id_column = "image_id",
                                                                                                 cell_type = immune_cell_type,
                                                                                                 cell_type_column = "cell_type",
                                                                                                 cell_x_column = "cell_x",
                                                                                                 cell_y_column = "cell_y",
                                                                                                 image_width_column = "image_width", 
                                                                                                 image_height_column = "image_height",
                                                                                                 distance_column = "distance_to_islet"
)}))
toc()

tic("Get patterns for all cells")
patterns_all <- as.solist(lapply(unique(colData(sce)$image_id),function(sample_id){get_ppp_for_image(sce,
                                                                                                     sample_id = sample_id, 
                                                                                                     sample_id_column = "image_id",
                                                                                                     cell_type = unique(colData(sce)$cell_type),
                                                                                                     cell_type_column = "cell_type",
                                                                                                     cell_x_column = "cell_x",
                                                                                                     cell_y_column = "cell_y",
                                                                                                     image_width_column = "image_width", 
                                                                                                     image_height_column = "image_height",
                                                                                                     distance_column = "distance_to_islet"
)}))
toc()
data$immune_pattern <-  patterns
data$all_pattern <-  patterns_all
rownames(data)<- unique(colData(sce)$image_id)

data$case_id <- unlist(strsplit(rownames(data),"[_]"))[c(T,F,F)]

by_donor<- split(colData(sce),colData(sce)$donor_type)
case_id_by_donor<- lapply(by_donor,function(x){unique(x$case_id)})

data$donor_type <- vector("list",nrow(data))
for(i in 1:length(case_id_by_donor)){
  data[data$case_id %in% case_id_by_donor[[i]],"donor_type"]<- names(case_id_by_donor)[[i]]
}
data$donor_type <- unlist(data$donor_type,recursive = F)


tic("Get quadrats")
quadrats <- lapply(data$all_pattern,function(x){
  
  marks(x)[marks(x)$distance <= 0,"distance"] <- 0
  
  return(x)
  
})

Z<- lapply(quadrats,function(x){
  marks(x)$cell_type <- NULL
  Z <- Smooth(x, sigma = bw.CvL)
  return(Z)
})

data$cov <- Z


data$tess <- lapply(data$cov,function(Z){
  b <- quantile(Z,seq(0,1,1/3))
  Zcut <- cut(Z,breaks = b,labels = c("islet","near_islet","far"))
  V <- tess(image = Zcut)
  return(V)
})

toc()  
### Get rid of distance marks
data$immune_pattern<- lapply(data$immune_pattern,function(x){marks(x)$distance <- NULL
levels(marks(x)) <- immune_cell_type
return(x)})


### Get intensities for quadrats/cell types
quadrat_intensities <- vector("list", nrow(data))
for (i in 1:nrow(data)){
  x <- split(data$immune_pattern[[i]],marks(data$immune_pattern[[i]]))
  
  result <- lapply(quadratcount(x,tess = data$tess[[i]]),intensity)
  quadrat_intensities[[i]] <- result 
  
}


names(quadrat_intensities) <- rownames(data)

### Fix data assign disease stage cell id tile image id and patient id 
a <-  unlist(quadrat_intensities,recursive = T)
id<- unlist(strsplit(names(a),"[.]"))[c(T,F,F)]
cell <- unlist(strsplit(names(a),"[.]"))[c(F,T,F)]
tile <- unlist(strsplit(names(a),"[.]"))[c(F,F,T)]
d<- data.frame(intensity = a,image_id = id ,cell_type = cell ,tile = tile)
d$case_id <- unlist(strsplit(d$image_id,"[_]"))[c(T,F,F)]
by_donor<- split(colData(sce),colData(sce)$donor_type)
case_id_by_donor<- lapply(by_donor,function(x){unique(x$case_id)})

d$donor_type <- vector("list",nrow(d))
for(i in 1:length(case_id_by_donor)){
  d[d$case_id %in% case_id_by_donor[[i]],"donor_type"]<- unlist(names(case_id_by_donor)[[i]])
}
d$donor_type <- unlist(d$donor_type)

t1d_stages <- c("NoDiabetes", "mAAb+", "RecentOnset", "LongDuration")
order_samples <- unique(colData(sce)[, c("donor_type", "case_id")])
order_samples <- order_samples[order(match(order_samples$donor_type, t1d_stages)), ]

d <-d %>%
  # as_tibble() %>%
  group_by(case_id, donor_type) %>%
  mutate(case_id = factor(case_id, levels = order_samples$case_id),
         donor_type = factor(donor_type, levels = t1d_stages))%>%
  arrange(case_id, donor_type)

### Split by cell type 

d_split<- split(d,d$cell_type)

### Get plots

my_theme <- theme(
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
my_comparisons = list( c("mAAb+", "RecentOnset"),c("NoDiabetes", "RecentOnset"),c("LongDuration", "RecentOnset") )
### Get plots
plots <- lapply(d_split,function(x){
  x  <- x %>%
    mutate(across(tile, factor, levels=c("islet","near_islet","far"))) 
  
  point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)
  
  p <- ggplot(as.data.frame(x), aes(x = donor_type, y = intensity,col = donor_type )) + 
    geom_violin() +geom_jitter(size = 0.2) + my_theme + 
    facet_wrap(~tile,labeller =  function(x){
    x <- NULL
    return(data.frame(labels = c("Islet", "Near Islet","Far")))
  })  + stat_summary(fun = "mean", geom = "point",mapping = aes(shape="mean"),color = "red",size = 1.2) +  scale_shape_manual(values=c(16)) +  scale_y_continuous(labels = point) + 
    labs(shape = NULL, x = "Case ID", y = "Absolute Density", col = "Disease Stage") + stat_compare_means(comparisons = my_comparisons, label.y = c(0.0001,0.0002,0.0003),show.legend = F,paired = F,method = "wilcox",symnum.args =   list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) + 
    stat_compare_means(show.legend = F,label.y = c(0.00075),size = 2.5) 
  
                                                                                                                                                                          
  return(p)
  
})

### Save plots
i <- 1
for(my_plot in plots){
  file <- paste("quadrat_densities(",as.character(i),").pdf")
  ggsave(filename =file,plot = my_plot,device= "pdf", width = 18,height = 12, units = "cm")
  i <- i + 1
}



#### Clark evans

data$clarkevans <- sapply(data$immune_pattern,function(x){return(clarkevans(x,correction = "Donnelly"))})


a <- as.data.frame(data[,c("case_id","clarkevans","donor_type")])


a <- a[is.finite(a$clarkevans),]

                      
p6 <- ggplot(data = a,mapping= aes(x = donor_type,y = clarkevans ,col = donor_type),show.legend = F)+ geom_boxplot(outlier.shape = NA,show.legend = T) + my_theme + geom_jitter(size=0.5,show.legend = T) + ylab("Clark Evans Index") + xlab("Disease Stage") +labs(col = "Disease Stage")+
  stat_compare_means(method = "anova",show.legend = FALSE,size = 3) + 
  geom_abline(aes(intercept=1, slope=0), linetype = "dashed") 


ggsave(filename ="clarkevans_bystage.pdf",plot = p6,device= "pdf", width = 18,height = 12, units = "cm")


### L functions

data$Limmune <- lapply(data$immune_pattern,function(x){return(Lest(x,correction = "isotropic",ratio = T))})


mask <- sapply(data$immune_pattern, npoints) >= 3



LP <- split(data[mask], f=data$donor_type[mask])


pool_LP <- lapply(LP,function(x){return(pool.anylist(x$Limmune,relabel=T))})

av_mean_rad <- mean(sqrt(((colData(sce)$cell_minor_axis_length/2)^2+(colData(sce)$cell_major_axis_length/2)^2)/2))*2

pdf(file = "L-function.pdf",height = 4.72441,width = 7.08661)
plot(as.anylist(pool_LP), ylim = c(0,20),cbind(pooliso,pooltheo,hiiso,loiso)-r~r,legend = F,shade=c("hiiso","loiso"),xlim =c(av_mean_rad,60),main = "") 
dev.off()

### Lcross function 




Lcross_beta_tc8 <-get_test_result(sce,image_number = "image_id",correction = "all",from = "Beta",to = "T_CD8",use.weights = T)

Lcross_beta_tc8<- Lcross_beta_tc8[!(Lcross_beta_tc8$case_id %in% c("6063","6135","6285","6418","6422","6289")),] # Ld remove 6289 removed because it had only one image that had more than 1 tc8 cell

pool_Lc_tc8 <- lapply(split(Lcross_beta_tc8$test,f = Lcross_beta_tc8$donor_type),function(x){pool(x)})

Lcross_beta_tc4 <-get_test_result(sce,image_number = "image_id",correction = "all",from = "Beta",to = "T_CD4",use.weights = T)

Lcross_beta_tc4<- Lcross_beta_tc4[!(Lcross_beta_tc4$case_id %in% c("6063","6135","6285","6418","6422","6289")),]

pool_Lc_tc4 <- lapply(split(Lcross_beta_tc4$test,f = Lcross_beta_tc4$donor_type),pool)



pdf("cross_tcd4.pdf",height = 4.72441,width = 7.08661)
plot.anylist(pool_Lc_tc4 ,cbind(pooliso,pooltheo,hiiso,loiso)-r~r,legend = F,shade=c("hiiso","loiso"),xlim =c(av_mean_rad,60), main = "",equal.scales = T)
dev.off()
pdf("cross_tcd8.pdf",height = 4.72441,width = 7.08661)
plot.anylist(pool_Lc_tc8,cbind(pooliso,pooltheo,hiiso,loiso)-r~r,legend = F,shade=c("hiiso","loiso"),xlim =c(av_mean_rad,60), main = "",equal.scales = T)
dev.off()


### Model


data_6534 <- data[data$case_id == "6534", ]
data_6534 <-data_6534[rownames(data_6534) != "6534_Immune_082",] # remove the corrupted image

my_model_6534 <- mppm(immune_pattern ~ marks*tess ,random = ~1 | id ,data = data_6534)



data_6289 <- data[data$case_id == "6289", ]


my_model_6289 <- mppm(immune_pattern ~ marks*tess ,random = ~1 | id ,data = data_6289)


preds <- predict(my_model_6534)
preds_ND <- predict(my_model_6289)
pdf("6534_density_by_type_im_059.pdf",height = 4.72441,width = 7.08661)
plot(density(split(data_6534$immune_pattern$`6534_Immune_059`,marks(data_6534$immune_pattern$`6534_Immune_059`)))[c("Macrophage_1","Macrophage_2","T_CD8","T_CD4")],main = "",ribbon = T)
dev.off()

my_ims <- anylist(Macrophage_1 = preds$cifMacrophage_1$`6534_Immune_059`,Macrophage_2 = preds$cifMacrophage_2$`6534_Immune_059`,T_CD8 = preds$cifT_CD8$`6534_Immune_059`,T_CD4 = preds$cifT_CD4$`6534_Immune_059`)
pdf("6534_trend_by_type_im_059.pdf",height = 4.72441,width = 7.08661)
plot(my_ims,main = "",ribbon = T)
dev.off()

pdf("6289_density_by_type_im_072.pdf",height = 4.72441,width = 7.08661)
plot.anylist(anylist(density(split(data$immune_pattern$`6289_Immune_072`,marks(data$immune_pattern$`6289_Immune_072`))[c("Macrophage_1","Macrophage_2","T_CD8","T_CD4")])),main = "",ribbon = T,mar.panel = c(1,1,1,1),main.panel = c("Macrophage1","Macrophage2","T CD4","T CD8"))
dev.off()


my_ims_ND <- anylist(Macrophage_1 = preds_ND$cifMacrophage_1$`6289_Immune_072`,Macrophage_2 = preds_ND$cifMacrophage_2$`6289_Immune_072`,T_CD8 = preds_ND$cifT_CD8$`6289_Immune_072`,T_CD4 = preds_ND$cifT_CD4$`6289_Immune_072`)
pdf("6289_trend_by_type_im_072.pdf",height = 4.72441,width = 7.08661)
plot(my_ims_ND,main = "",ribbon = T)
dev.off()