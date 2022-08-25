library(spatstat)
library(S4Vectors)
library(SingleCellExperiment)
get_local_intensity <- function(spe, ROI_size, cell_type , cell_type_column, distance_column, sample_id_column ){
  
  mask <-  colData(spe)[,cell_type_column] == cell_type
  
  outside_mask <-  (colData(spe)[,distance_column] > 0 & 
                      colData(spe)[,distance_column] <= ROI_size)
  inside_mask <-  colData(spe)[,distance_column] <= 0
  
  In_ROI <- outside_mask | inside_mask
  

  for (sample_id in levels(as.factor(as.character(colData(spe)[,sample_id_column])))) {
    in_image_mask <- colData(spe)[,sample_id_column] == sample_id
    colData(spe)[in_image_mask, "Proportion_Score"] <- sum((mask & in_image_mask) & 
                                                             In_ROI)/sum(in_image_mask & In_ROI) 
    
  }
  
  
  return(spe[["Proportion_Score"]])
  
}


get_global_intensity <- function(spe, cell_type , cell_type_column, distance_column, sample_id_column){
  
  mask <-  colData(spe)[,cell_type_column] == cell_type
  
  for (sample_id in levels(as.factor(as.character(colData(spe)[,sample_id_column])))) {
    in_image_mask <- colData(spe)[,sample_id_column] == sample_id
    colData(spe)[in_image_mask, "Global_Intensity"] <- sum((mask & in_image_mask)) / 
     sum(in_image_mask)
  }
  
  
  return(spe[["Global_Intensity"]])
  
}



get_distance_score_reciprocal<- function(spe,cell_type , cell_type_column, distance_column, sample_id_column ){
  
  type_mask <- colData(spe)[,cell_type_column] == cell_type
  
  outside_mask <-  colData(spe)[,distance_column] > 0 & type_mask
  inside_mask <-  colData(spe)[ ,distance_column] <= 0 & type_mask
  
  colData(spe)[,"scores"] <- NaN 
  
  colData(spe)[outside_mask,"scores"] <- 1/colData(spe)[outside_mask,distance_column]
  colData(spe)[inside_mask, "scores"] <- 1
  
  
  dataframe_list <- as.list(split(colData(spe), f = colData(spe)[[sample_id_column]]))
  for(image_data in dataframe_list){
    
    colData(spe)[colData(spe)[,sample_id_column] == unique(image_data[,sample_id_column]),"mean_scores"] <- mean(image_data[,"scores"] , na.rm=T)
  }
  
  
  return(colData(spe)[["mean_scores"]])
  
}


# Distance scoring using normal distrubution (score = dnorm(distance/max(distance),mean = 0, sd = 0.4))

get_distance_score_normal<- function(spe,cell_type , cell_type_column, distance_column, sample_id_column ){
  
  type_mask <- colData(spe)[,cell_type_column] == cell_type
  
  
  outside_mask <-  colData(spe)[,distance_column] > 0 & type_mask
  inside_mask <-  colData(spe)[ ,distance_column] <= 0 & type_mask
  
  colData(spe)[,"scores"] <- NaN 
  normalized_distance <- (colData(spe)[outside_mask, distance_column] - min(colData(spe)[outside_mask, distance_column]))/( max(colData(spe)[outside_mask, distance_column])- min(colData(spe)[outside_mask, distance_column]))

  colData(spe)[outside_mask,"scores"] <- dnorm(normalized_distance, mean = 0, sd = 0.4)
  colData(spe)[inside_mask, "scores"] <- dnorm(0,mean= 0,sd=0.4)
  
  
  
  dataframe_list <- as.list(split(colData(spe), f = colData(spe)[[sample_id_column]]))
  for(image_data in dataframe_list){
    
    colData(spe)[colData(spe)[,sample_id_column] == unique(image_data[,sample_id_column]),"mean_scores"] <- mean(image_data[,"scores"] , na.rm=T)
  }
  
  
  return(colData(spe)[["mean_scores"]])
  
}



get_ppp_for_image <- function(spe, 
                              sample_id, 
                              sample_id_column,
                              cell_x_column,
                              cell_y_column,
                              cell_type,
                              cell_type_column,
                              image_width_column,
                              image_height_column,
                              distance_column
)
{
  
  image_number_mask <- colData(spe)[[sample_id_column]] == sample_id
  
  cell_type_mask <- (colData(spe)[[cell_type_column]] %in% cell_type)
  
  x_range  <- c(0,unique(colData(spe)[image_number_mask,image_width_column]))
  
  
  
  y_range <- c(0,unique(colData(spe)[image_number_mask ,image_height_column]))
  
  markings_type <-factor(colData(spe)[image_number_mask & cell_type_mask, cell_type_column])
  markings_distance <- colData(spe)[image_number_mask & cell_type_mask, distance_column]
  
  
  
  pattern <- ppp(x = colData(spe)[image_number_mask & cell_type_mask,cell_x_column],
                 y = colData(spe)[image_number_mask & cell_type_mask,cell_y_column], 
                 xrange = x_range,
                 yrange = y_range,
                 marks = data.frame(cell_type = markings_type, distance = markings_distance))
  
  
  return(pattern)
  
  
  
}


get_test_result <- function(data,
                            cell_type="cell_type",
                            cell_x="cell_x",
                            cell_y="cell_y",
                            image_width="image_width",
                            image_height="image_height",
                            image_number="image_number",group_by = "donor_type",
                            from,
                            to,
                            correction = "all",
                            test = spatstat.core::Lcross,selection = 1,use.weights = T){
  
  ## Getting point patterns form SpatialExperiment/SingleCellExperiment
  grouped_colData_sce <- split(x= colData(data),f = colData(data)[[image_number]])
  
  
  patterns_all <- lapply(grouped_colData_sce,function(data){
    return(ppp(x = data[[cell_x]],y = data[[cell_y]],xrange = c(0,unique(data[[image_width]])),yrange = c(0,unique(data[[image_height]])),marks = factor(data[[cell_type]])))})
  
  names(patterns_all) <- unique(colData(data)[[image_number]])
  
  
  
  has_from<- sapply(patterns_all,function(x){return(ifelse(sum(from  == marks(x)) > selection,T,F))})
 
  
  has_to<- sapply(patterns_all,function(x){return(ifelse(sum(to == marks(x))>selection,T,F))})

  
  
 
  
  result_test <- lapply(patterns_all[has_to & has_from],function(x){return(test(X = x,i=  from,j= to,correction = correction))})
  names(result_test) <- names(patterns_all)[has_to & has_from]
  
  
  if (use.weights){
    weights <- sapply(patterns_all, npoints)
    weights <- weights[has_to & has_from]
    result <- hyperframe(test = as.anylist(result_test),weights = weights)
    result$case_id <- unlist(strsplit(rownames(result),"[_]"))[c(T,F,F)]
    by_donor<- split(colData(sce),colData(sce)[group_by])
    case_id_by_donor<- lapply(by_donor,function(x){unique(x$case_id)})
    result$donor_type <- vector("list",nrow(result))
    for(i in 1:length(case_id_by_donor)){
      result[result$case_id %in% case_id_by_donor[[i]],"donor_type"]<- unlist(names(case_id_by_donor)[[i]])
    }
    result$donor_type <- unlist(result$donor_type)
    return(result)
  }else{
    weights <- NULL
    result <- hyperframe(test = as.anylist(result_test))
    result$case_id <- unlist(strsplit(rownames(result),"[_]"))[c(T,F,F)]
    by_donor<- split(colData(sce),colData(sce)[group_by])
    case_id_by_donor<- lapply(by_donor,function(x){unique(x$case_id)})
    result$donor_type <- vector("list",nrow(result))
    for(i in 1:length(case_id_by_donor)){
      result[result$case_id %in% case_id_by_donor[[i]],"donor_type"]<- unlist(names(case_id_by_donor)[[i]])
    }
    result$donor_type <- unlist(result$donor_type)
    return(result)
    
  }
  
}


get_covariate <-      function(spe,
                               sample_id_column,
                               cell_type,
                               cell_type_column, 
                               cell_x_column, 
                               cell_y_column,
                               image_width_column, 
                               image_height_column, 
                               distance_column)
{
  
  cell_type_all <- unique(colData(spe)[[cell_type_column]])
  
  patterns_all <- lapply(unique(colData(spe)[[sample_id_column]]),
                         function(x){get_ppp_for_image(spe,
                                                       sample_id = x, 
                                                       sample_id_column = sample_id_column,
                                                       cell_type = cell_type_all,
                                                       cell_type_column = cell_type_column,
                                                       cell_x_column = cell_x_column,
                                                       cell_y_column = cell_y_column,
                                                       image_width_column = image_width_column, 
                                                       image_height_column = image_height_column,
                                                       distance_column = distance_column  )
                       })
  
  patterns <- lapply(patterns_all ,
                     function(x){return(split(x , f =  "cell_type",drop = T)[cell_type])})
  
  hpf <- hyperframe(all = patterns_all,type = patterns )
  
  quadrats <- lapply(hpf$all,function(x){
 
    
    
    marks(x)$cell_type <- NULL
    Z <- Smooth(x, sigma = bw.CvL)
    b <- quantile(Z,seq(0,1,0.1))
    Zcut <- cut(Z,breaks = b,labels = 1:10)
    V <- tess(image = Zcut)
    return(V)
  })
  hpf$tess <- quadrats
  

  
    
  return(hpf)
}


### Unnecessary Functions

min_max_norm <- function(x){
 
     normalized <- (x - min(x))/(max(x)- min(x))

  return(normalized)
}

pool_by_patient_id <- function(result){
  
  z <- unlist(result,recursive = FALSE)
  b<- unlist(sapply(names(z),function(x){return(strsplit(x,"[.]"))}))[c(T,F)]
  c <- hyperframe(test = unlist(result,recursive = FALSE) ,patient_id = b)
  rownames(c) <-names(z)
  hpf<- split(c,c$patient_id)
  

  
  return(lapply(hpf,function(x){pool.anylist(x$test)}))

}


## Hard coded function for my use 
plot_score<- function(sce,score_columns,scales = "fixed",normalize = TRUE,theme =theme_bw(),log_transform = F){
  scores <- as.data.frame(colData(sce)) %>% select("image_id","donor_type","case_id",ends_with(score_columns))
  
  scores <- scores[!duplicated(scores$image_id),]
  
  scores[is.na(scores)] <- 0
  if (normalize){
  scores<- data.frame(scores[,1:3],lapply(scores[,4:ncol(scores)],min_max_norm))
    
    
  }
  
  df_long <- scores %>% pivot_longer(cols = ends_with(score_columns),names_to = "score_type",values_to = "score")
  if(log_transform){
  df_long$score <- -log2(df_long$score) }
  p1 <- ggplot(df_long,aes(y= score,x = donor_type,col=donor_type)) + geom_violin() + facet_wrap(~score_type,scales = scales) + theme + geom_jitter(size=0.5) + stat_summary(fun = "mean",
               geom = "point",position=position_nudge(x=0.08),
               color = "red",show.legend = F)+ stat_summary(fun = "median",
               geom = "point",
               color = "green",show.legend = F)
  return(p1)
}

