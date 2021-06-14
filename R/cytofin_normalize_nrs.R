#' Normalize CyTOF data panel
#'
#' This function normalize CyTOF data panel using internal anchors.
#'
#' @param sample_metadata_filename 
#' @param panel_filename 
#' @param input_file_dir 
#' @param val_file_dir 
#' @param output_file_dir 
#' @param nchannels 
#'
#' @return normalized files with specified panel
#' @export

cytofin_normalize_nrs <- 
  function(
    sample_metadata_filename, 
    panel_filename, 
    input_file_dir, 
    val_file_dir = "none", 
    output_file_dir, 
    nchannels
  ) {
    dir.create(output_file_dir)
    mode <- "meanshift_bulk"
    #read metadata table
    getExtension <- function(file){ 
      ex <- strsplit(basename(file), split="\\.")[[1]]
      return(ex[-1])
    } 
    
    
    if (getExtension(sample_metadata_filename)=="xlsx"){
      md <- readxl::read_excel(sample_metadata_filename, col_names = TRUE)
    }else if (getExtension(sample_metadata_filename)=="csv"){
      md <- read.csv(sample_metadata_filename)
    }
    
    #parse standard panel
    
    #read excel file
    if (getExtension(panel_filename)=="xlsx"){
      ref_panel <- readxl::read_excel(panel_filename)
    }else if (getExtension(panel_filename)=="csv"){
      ref_panel <- read.csv(panel_filename)
    }
    
    #ref_panel <- readxl::read_excel(panel_filename)
    
    (lineage_markers <- as.character(ref_panel$desc[ref_panel$Lineage == 1]))
    (functional_markers <- as.character(ref_panel$desc[ref_panel$Functional == 1]))
    all_markers <- c(lineage_markers, functional_markers)
    
    #transformation function
    norm_2 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] - mean(mean_ctr[selected_markers]) + mean(mean_uni[selected_markers])
      return(z)
    }#meanshift bulk
    
    
    # Define a function for arsinh transformation (based on definition from Nikolay).
    asinhNik <- function(value) {
      value <- value - 1
      for(i in 1:length(value)) {
        if((value[i] < 0) | is.na(value[i])) value[i] <- rnorm(1, mean = 0, sd = 0.01)
      }  # can also use max(value, 0)
      value <- value / 5  # Co-factor value
      value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))  
      return(value)
    }
    
    # Reverses arsinh transformation with a cofactor of 5.
    rev_asinh <- function(x) {
      e.x <- exp(x)
      rev.x <- (e.x^2 - 1) / (2 * e.x)
      x <- rev.x * 5
      return(x)
    }
    
    if (mode == "meanshift"){
      norm <- norm_1 
    }else if (mode == "meanshift_bulk"){
      norm <- norm_2
    }else if (mode == "variance"){
      norm <- norm_3
    }else if (mode == "z_score"){
      norm <- norm_4
    }else if (mode == "beadlike"){
      norm <- norm_5
    }
    
    ##compute stable channels
    NRS <- function(x, ncomp = 3){
      pr <- prcomp(x, center = TRUE, scale. = FALSE)
      score <- rowSums(outer(rep(1, ncol(x)),
                             pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
      return(score)
    }
    
    fcs <- flowCore::read.flowSet(paste0(input_file_dir,"homogenized_",md$filename), transformation = FALSE, truncate_max_range = FALSE)
    
    asinhNik <- flowCore::arcsinhTransform(a=0,b=0.2)
    colname <- flowCore::colnames(fcs)
    tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
    fcs_asinh <- flowCore::transform(fcs, tlist)
    expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
    mean_uni <- apply(expr, 2, mean)
    var_uni <- apply(expr, 2, var)
    nrs_sample <- flowCore::fsApply(fcs_asinh[, all_markers], NRS, use.exprs = TRUE)
    colnames(nrs_sample) <- as.character(ref_panel$range[match((colnames(nrs_sample)), ref_panel$desc)])
    nrs <- colMeans(nrs_sample, na.rm = TRUE)
    
    nrs_sample <- data.frame(nrs_sample)
    markers_ord <- names(sort(nrs, decreasing = TRUE))
    nrs_sample <- data.frame(nrs_sample)
    nrs_sample$sample_id <- rownames(nrs_sample)
    
    #colnames(nrs_sample) <- ref_panel$range[match((colnames(nrs_sample)), ref_panel$desc)]
    ggdf <- reshape2::melt(nrs_sample, id.var = "sample_id", value.name = "nrs", variable.name = "antigen")
    ggdf$antigen <- factor(ggdf$antigen, levels = markers_ord)
    p <- ggplot2::ggplot(ggdf, ggplot2::aes(x = antigen, y = nrs)) +
      ggplot2::geom_point(ggplot2::aes(color = sample_id), alpha = 0.9,
                          position = ggplot2::position_jitter(width = 0.3, height = 0)) +
      ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
      ggplot2::stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    print(p)
    
    #select top 5 markers for calibration
    selected_markers <- names(sort(nrs, decreasing = FALSE))[1:nchannels]
    selected_markers<- as.character(ref_panel$desc[match(selected_markers, ref_panel$range)])
    ##
    
    
    for (i in 1:length(md$filename)){
      #calculate adjustment parameters from control plate
      #i=1
      cat(md$filename[i], "\n")
      #filename_ctr <- md_control$filename[which((md_control$plate_number==md$plate_number[i])&(md_control$cohort==md$cohort[i]))]
      filename_ctr <- md$filename[i] 
      fcs <- flowCore::read.flowSet(paste0(input_file_dir,"homogenized_",filename_ctr), transformation = FALSE, truncate_max_range = FALSE)
      asinhNik <- flowCore::arcsinhTransform(a=0,b=0.2)
      colname <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_ctr <- flowCore::fsApply(fcs_asinh[,selected_markers], flowCore::exprs)
      mean_ctr <- apply(expr_ctr, 2, mean)
      mean_ctr_mean <- mean(mean_ctr)
      var_ctr <- apply(expr_ctr, 2, var)
      var_ctr_mean <- mean(var_ctr)
      markername <- flowCore::pData(flowCore::parameters(fcs[[1]]))$desc
      
      #normalize the control plate
      expr_ctr_norm <- t(flowCore::fsApply(fcs_asinh, function(x) apply(x, 1, norm), use.exprs=TRUE))
      mean_ctr_norm <- apply(expr_ctr_norm[,selected_markers], 2, mean)
      var_ctr_norm <- apply(expr_ctr_norm[,selected_markers], 2, var)
      mean_ctr_norm_mean <- mean(mean_ctr_norm)
      var_ctr_norm_mean <- mean(var_ctr_norm)
      
      #normalize the target plate
      ##before
      filename <- md$filename[i]
      fcs <- flowCore::read.flowSet(paste0(input_file_dir,"homogenized_",filename), transformation = FALSE, truncate_max_range = FALSE)
      asinhNik <- flowCore::arcsinhTransform(a=0,b=0.2)
      colname <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_b4norm <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
      mean_b4norm <- apply(expr_b4norm, 2, mean)
      var_b4norm <- apply(expr_b4norm, 2, var)
      mean_b4norm_mean <- mean(mean_b4norm)
      var_b4norm_mean <- mean(var_b4norm)
      ##after
      expr_norm <- t(flowCore::fsApply(fcs_asinh, function(x) apply(x, 1, norm), use.exprs=TRUE))
      mean_norm <- apply(expr_norm, 2, mean)
      var_norm <- apply(expr_norm, 2, var)
      mean_norm_mean <- mean(mean_norm)
      var_norm_mean <- mean(var_norm)
      fcs_norm <- flowCore::flowFrame(expr_norm)
      #normalization completed, reverse transformation
      tlist2 <- flowCore::transformList(from = colname, tfun = rev_asinh)
      fcs_asinh_rev <- flowCore::transform(fcs_norm, tlist2)
      flowCore::pData(flowCore::parameters(fcs_asinh_rev))$desc <- flowCore::pData(flowCore::parameters(fcs_asinh[[1]]))$desc
      fcs_name <- paste0(output_file_dir,"normalized_",filename)
      flowCore::write.FCS(fcs_asinh_rev, fcs_name)
      
      #compare to validation
      if (val_file_dir != "none"){
        filename_val <- md$validation[i]
        fcs <- flowCore::read.flowSet(paste0(val_file_dir,filename_val), transformation = FALSE, truncate_max_range = FALSE)
        asinhNik <- flowCore::arcsinhTransform(a=0,b=0.2)
        colname <- flowCore::colnames(fcs)
        tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
        fcs_asinh <- flowCore::transform(fcs, tlist)
        expr_val <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
        mean_val <- apply(expr_val, 2, mean)
        var_val <- apply(expr_val, 2, var)
        mean_val_mean <- mean(mean_val)
        var_val_mean <- mean(var_val)
        
      }
      
      
      #visualize antigen panel comparing universal vs plate
      par(mfrow=c(2,4))
      
      len<-length(mean_uni[all_markers])
      
      #expression (mean)
      #plot 1
      plot(mean_uni[all_markers], col='red', xlab='antigen', ylab='universal expression (mean)', xlim=c(0,len), ylim=c(-5,10), main='overall', cex.main=1)
      legend(1,10,legend=c("universal"),col=c("red"), lty=1:2, cex=0.8)
      
      #plot 2
      plot(mean_ctr[all_markers], col='cyan', xlab='antigen', ylab='control expression (mean)', xlim=c(0,len), ylim=c(-5,10), main=filename_ctr, cex.main=0.8)
      par(new=TRUE)
      plot(mean_ctr_norm[all_markers], col='blue', xlab='antigen', ylab='control expression (mean)', xlim=c(0, len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("blue", "cyan"), lty=1:2, cex=0.8)
      
      #plot 3
      plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='plate expression (mean)', xlim=c(0,len), ylim=c(-5,10), main=filename, cex.main=0.8)
      par(new=TRUE)
      plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='plate expression (mean)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("darkgreen", "green"), lty=1:2, cex=0.8)
      
      #plot 4
      if (val_file_dir != "none"){
        plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_val[all_markers], col='purple', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized", "validation"),col=c("green", "darkgreen", "purple"), lty=1:2, cex=0.8)
      }else{
        plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized"), col=c("green", "darkgreen"), lty=1:2, cex=0.8)
      }
      
      
      #expression (std)
      #plot 5
      plot(sqrt(var_uni[all_markers]), col='red', xlab='antigen', ylab='universal expression (std)', xlim=c(0,len), ylim=c(-5,10), main='overall', cex.main=1)
      legend(1,10,legend=c("universal"),col=c("red"), lty=1:2, cex=0.8)
      
      #plot 6
      plot(sqrt(var_ctr[all_markers]), col='cyan', xlab='antigen', ylab='control expression (std)', xlim=c(0,len), ylim=c(-5,10), main=filename_ctr, cex.main=0.8)
      par(new = TRUE)
      plot(sqrt(var_ctr_norm[all_markers]), col='blue', xlab='antigen', ylab='control expression (std)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1, 10, legend = c("normalized", "original"),col=c("blue", "cyan"), lty = 1:2, cex = 0.8)
      
      #plot 7
      plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='plate expression (std)', xlim=c(0,len), ylim=c(-5,10), main=filename, cex.main=0.8)
      par(new = TRUE)
      plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='plate expression (std)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("darkgreen", "green"), lty=1:2, cex=0.8)
      
      #plot 8
      if (val_file_dir != "none") {
        plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(var_val[all_markers], col='purple', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized", "validation"), col=c("green", "darkgreen", "purple"), lty=1:2, cex=0.8)
      } else {
        plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized"), col=c("green", "darkgreen"), lty=1:2, cex=0.8)
      }
      
      
    }
  }
