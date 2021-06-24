volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "red",
                                     "royalblue3")) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}

runProgenyFast <- function(df,weight_matrix,k = 10000, z_scores = T, get_nulldist = F)
{
  resList <- list()
  if(get_nulldist)
  {
    nullDist_list <- list()
  }
  
  for(i in 2:length(df[1,]))
  {
    current_df <- df[,c(1,i)]
    current_df <- current_df[complete.cases(current_df),]
    t_values <- current_df[,2]
    
    current_weights <- weight_matrix
    
    names(current_df)[1] <- "ID"
    names(current_weights)[1] <- "ID"
    
    common_ids <- merge(current_df, current_weights, by = "ID")
    common_ids <- common_ids$ID
    common_ids <- as.character(common_ids)
    
    row.names(current_df) <- current_df$ID
    current_df <- as.data.frame(current_df[common_ids,-1])
    
    row.names(current_weights) <- current_weights$ID
    current_weights <- as.data.frame(current_weights[common_ids,-1])
    
    current_mat <- as.matrix(current_df)
    current_weights <- t(current_weights)
    
    scores <- as.data.frame(current_weights %*% current_mat)
    
    null_dist_t <- replicate(k, sample(t_values,length(current_mat[,1]), replace = F))
    
    null_dist_scores <- current_weights %*% null_dist_t
    
    if(get_nulldist)
    {
      nullDist_list[[i-1]] <- null_dist_scores
    }
    
    if(z_scores)
    {
      scores$mean <- apply(null_dist_scores,1,mean)
      scores$sd <- apply(null_dist_scores,1,sd)
      resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
    else
    {
      for(j in 1:length(weight_matrix[,-1]))
      {
        ecdf_function <- ecdf(null_dist_scores[j,])
        scores[j,1] <- ecdf_function(scores[j,1])
      }
      score_probas <- scores*2-1
      
      resListCurrent <- score_probas[,1]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
  }
  names(resList) <- names(df[,-1])
  resDf <- as.data.frame(resList)
  if(get_nulldist)
  {
    names(nullDist_list) <- names(df[,-1])
    return(list(resDf, nullDist_list))
  }
  else
  {
    return(resDf)
  }
  
}
