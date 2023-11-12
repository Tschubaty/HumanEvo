



df_empirical <-
  setNames(
    data.frame(matrix(ncol = 6, nrow = 0)),
    c(
      "permuation_order",
      "min_delta",
      "min_log_p",
      "count",
      "data_type",
      "hemming.distance_permuation"
    )
  )


minus_log_alpha <- 6.5

for (min_delta in seq(from = 0.05, to = 0.2, by = 0.01)) {
  start_time <- Sys.time()
  print(min_delta)
  
  p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
  #print(p_expr)
  alpha <- eval(parse(text = p_expr))
  
  real_peak_results_filename <- file.path(
    OUTPUT_FOLDER,
    paste(
      "real_peak_results",
      ANNOTATION,
      N_SAMPLES,
      "min_type_delta",
      min_delta,
      "rds",
      sep = "."
    )
  )
  
  
  
  if (file.exists(real_peak_results_filename)) {
    real_peak_results <- readRDS(file = real_peak_results_filename)
  } else{
    real_peak_results <-
      find_best_CpG_for_peak(
        peak_coord = peak_coordinates,
        KW.test.stat = df_type_real_results,
        min_delta = min_delta
      )
    saveRDS(object = real_peak_results,
            file = file.path(
              OUTPUT_FOLDER,
              paste(
                "real_peak_results",
                ANNOTATION,
                N_SAMPLES,
                "min_type_delta",
                min_delta,
                "rds",
                sep = "."
              )
            ))
  }
  
  
  real_peak_results_complete <-
    real_peak_results[real_peak_results$n_complete > 0,]
  
  # real data
  df_empirical <- rbind(
    df_empirical,
    data.frame(
      permuation_order = paste0(1:37, collapse = "."),
      min_delta = min_delta,
      min_log_p = minus_log_alpha,
      count = sum(
        !is.na(real_peak_results_complete$best_p_val) &
          real_peak_results_complete$best_p_val < alpha
      ),
      data_type = "real",
      hemming.distance_permuation = 0
    )
  )
  
  #########################################################################
  
  
  # peak_coord <- peak_coordinates
  # KW.test.stat <- df_type_real_results
  # min_d <- 0.1
  #
  
  find_best_CpG_for_peak <-
    function(peak_coord,
             KW.test.stat,
             min_d) {
      df_best_CpG <-
        as.data.frame(matrix(
          nrow = nrow(peak_coord),
          ncol = ncol(KW.test.stat) + 1
        ))
      
      for (chr in CHR_NAMES) {
        #print(chr)
        #chr <- CHR_NAMES[1]
        index_chr <- chr == df_bed6$chrom
        df_bed6_chr <- df_bed6[index_chr,]
        KW.test.stat_chr <- KW.test.stat[index_chr,]
        index_peak_chr <- chr == peak_coord$chr
        peak_coord_chr <- peak_coord[index_peak_chr,]
        df_best_CpG_chr <-
          as.data.frame(matrix(
            nrow = nrow(peak_coord_chr),
            ncol = ncol(KW.test.stat_chr) + 1
          ))
        
        
        for (r in 1:nrow(peak_coord_chr)) {
          # r <- 2000
          peak <- peak_coord_chr[r,]
          if (peak$n_complete > 0) {
            indices_min_delta <-
              peak$start <= df_bed6_chr$start &
              df_bed6_chr$end <= peak$end &
              !is.na(KW.test.stat_chr$type_mean_delta) &
              KW.test.stat_chr$type_mean_delta >= min_d
            
            n_min_delta <- sum(indices_min_delta)
            
            if (n_min_delta > 0) {
              # CpG postiions satisfy conditions
              CpG_positions <- df_bed6_chr$start[indices_min_delta]
              
              p_vals <- KW.test.stat_chr$p_val[indices_min_delta]
              statistics <-
                KW.test.stat_chr$statistic[indices_min_delta]
              type_mean_deltas <-
                KW.test.stat_chr$type_mean_delta[indices_min_delta]
              
              # best position of those
              best_index <- which.min(p_vals)
              # best index position
              best_CpG_position <-
                CpG_positions[best_index]
              
              df_best_CpG_chr[r,] <-
                as.numeric(c(
                  best_CpG_position,
                  p_vals[best_index],
                  statistics[best_index],
                  type_mean_deltas[best_index]
                ))
              
            }
          }
        }
        df_best_CpG[index_peak_chr,] <- df_best_CpG_chr
      }
      colnames(df_best_CpG) <-
        c("best_position",
          "best_p_val",
          "best_statistic",
          "best_type_mean_delta")
      return(cbind(peak_coord, df_best_CpG))
    }
  
  
  
  
  # start_time <- Sys.time()
  # real_peak_results <-
  #   find_best_CpG_for_peak(
  #     peak_coord = peak_coordinates,
  #     KW.test.stat = df_type_real_results,
  #     min_d = min_delta
  #   )
  #   end_time <- Sys.time()
  #   print(end_time - start_time)
  

  # p <- ggplot(data = df_empirical,
  #             mapping = aes(x = min_log_p, y = count, color = data_type)) + geom_point() +
  #   ggtitle(min_delta, "min_delta") +
  #   theme_minimal()
  
  # ggsave(plot = p,
  #        filename =   file.path(
  #          OUTPUT_FOLDER,
  #          "plot",
  #          paste("min_delta",
  #                min_delta,
  #                "n_significant_peaks.png")
  #        ))
  
  # df_empirical_means <-
  #   aggregate(count ~ min_log_p, df_empirical[df_empirical$data_type == "sim", ], FUN =
  #               mean)
  # df_empirical_means$real_count <-
  #   df_empirical$count[df_empirical$data_type == "real"]
  #
  # p <- ggplot(data = df_empirical_means,
  #             mapping = aes(x = min_log_p, y = real_count / count)) + geom_point() + theme_minimal() +
  #   ggtitle(min_delta, "min_delta")
  
  # ggsave(plot = p,
  #        filename =   file.path(
  #          OUTPUT_FOLDER,
  #          "plot",
  #          paste(
  #            "min_delta",
  #            min_delta,
  #            "enrichement_significant_peaks.png"
  #          )
  #        ))
}

###################################### all data 
df_empirical_big <- data.frame()

for (min_delta in seq(from = 0.05, to = 0.2, by = 0.01)) {
  print(min_delta)
  df_empirical <- readRDS(file =  file.path(
    OUTPUT_FOLDER,
    "simulation",
    "min_delta",
    min_delta,
    paste(
      "auumulated_sim_peak_results.hg19" ,
      min_delta,
      "mean_delta.rds",
      sep = "."
    )
  ))
  df_empirical_big <- rbind(df_empirical_big, df_empirical)
}

ggplot(data = df_empirical_big[df_empirical_big$min_log_p == 6.5, ],
       mapping = aes(x = min_delta, y = count, color = data_type)) + 
  geom_point() +
  theme_minimal()

points ~ team + conf

df_empirical_means <-
  aggregate(count ~ min_delta + min_log_p, df_empirical_big[df_empirical_big$data_type == "sim",], FUN =
              mean)


# tail()
# 
# df_empirical_means$real_count <-
#   df_empirical$count[df_empirical$data_type == "real"]

p <- ggplot(data = df_empirical_means,
            mapping = aes(x = min_log_p, y = real_count / count)) + geom_point() + theme_minimal() +
  ggtitle(min_delta, "min_delta")


df_empirical_means <- df_empirical_means[order(df_empirical_means$min_delta,df_empirical_means$min_log_p),]

colnames(df_empirical_means)[colnames(df_empirical_means) == 'count'] <- 'mean_sim_count'

df_empirical_real <- df_empirical_big[df_empirical_big$data_type == "real",]
df_empirical_real <- df_empirical_real[order(df_empirical_real$min_delta,df_empirical_real$min_log_p),]

df_empirical_means$real_count <- df_empirical_real$count

df_empirical_means$enrichment <- df_empirical_means$real_count/df_empirical_means$mean_sim_count

library(plotly)

 plot_ly(df_empirical_means, x = ~min_delta , y = ~min_log_p, z = ~enrichment )

 search_df <- df_empirical_means[df_empirical_means$min_log_p < 7.6,]

 plot_ly(search_df , x = ~min_delta , y = ~min_log_p, z = ~enrichment ,marker = list(color = ~enrichment, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
 
 search_df[which.max(search_df$enrichment),]
 
# ggsave(plot = p,
#        filename =   file.path(
#          OUTPUT_FOLDER,
#          "plot",
#          paste(
#            "min_delta",
#            min_delta,
#            "enrichement_significant_peaks.png"
#          )
#        ))
# end_time <- Sys.time()
# print(end_time - start_time)
 
 meta39 <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/gatherData/meta.rds")
 
 meta39$sample <- factor(x =  meta39$sample,levels = meta39$sample[order(meta39$age_mean_BP)]) 
 
 
 p <-  ggplot(data = meta39,
              mapping = aes(x = sample, y = age_mean_BP, color = Type)) + geom_point() +
   theme_minimal()+
   theme(axis.text.x = element_text(
     angle = 90,
     vjust = 0.5,
     hjust = 1,
     color = sample
   ))
 
 
  g <- ggplot_build(p)
 g$data[[1]]["colour"]
 
 
 ggplot(data = meta39,
        mapping = aes(x = sample, y = age_mean_BP, color = Type)) + geom_point() +
   theme_minimal()+
   theme(axis.text.x = element_text(
     angle = 90,
     vjust = 0.5,
     hjust = 1#,
     #color = g$data[[1]]["colour"]
   ))
 

 
