#' Generate A Histograms for A Seurat Object Feature
#'
#' Creates a histogram for a given feature in a Seurat object for QC purpose with interactive Plotly support.
#' #FEATURE_HIST COURTESY OF MATT CANYON via https://github.com/kidcancerlab/rrrSingleCellUtils
#. modified to change text angle to 45
#'
#' @param sobject A Seurat object.
#' @param feature A character of feature to plot.
#' @param cutoff Cutoff for x-axis in plot.
#' @param nbins Number of bins to use for the histogram. Default to 200.
#' @param interactive Logical, if TRUE, returns an interactive Plotly plot. Default is TRUE.
#' 
#' @return A ggplot2 object (if plotly = FALSE) or a Plotly interactive plot (if plotly = TRUE).
#' @examples
#' \dontrun{
#' single_histogram(sobject, features = c("nFeature_RNA", "nCount_RNA"))
#' }
#' @export

single_histogram <- function(sobject,
                             feature,
                             # a single feature
                             cutoff = 0.999,
                             # percentile cutoff
                             nbins = 200,
                             interactive = TRUE) {
  
  if (!feature %in% colnames(sobject@meta.data)) {
    stop(paste0("No match found in the metadata for ", feature))
  }
  
  
  temp_data <- Seurat::FetchData(sobject, vars = feature) %>%
    dplyr::rename(value = 1)
  
  x_cutoff <- quantile(temp_data$value, probs = cutoff)
  
  # calculate data in ggplot to be passed to plotly
  p <- ggplot2::ggplot(temp_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = nbins)
  
  
  if (interactive) {
    # construct a 0 gap bar plot in plotly instead of using ggplot2
    # so that we can have cdf on hover
    plot_data <- ggplot2::layer_data(p) %>%
      dplyr::mutate(cumprob = cumsum(count) / sum(count)) %>%
      dplyr::mutate(width = xmax - xmin)
    # bin width should be a single number
    bin_width <- plot_data$width[1]
    
    
    p <- plotly::plot_ly(
      data = plot_data,
      type = "bar",
      x = ~x,
      y = ~count,
      width = bin_width,
      marker = list(color = "#333333"), # dark grey
      customdata = ~cumprob,
      hovertemplate = paste0(
        "<b>Bin:</b> %{x:,d}<br>",
        "<b>Count:</b> %{y:,d}<br>",
        "<b>Cumulative Probability:</b> %{customdata:.2%}",
        "<extra></extra>"
      )
    ) %>%
      plotly::layout(
        margin = list(t = 100),
        bargap = 0,
        hoverlabel = list(yanchor = "top"),
        xaxis = list(
          range = c(0, x_cutoff),
          title = 'Intervals',
          showspikes = TRUE, 
          spikedistance = -1, 
          spikemode = "across",
          spikethickness = 1,
          spikedash = "dash",
          spikesnap = "data",
          spikecolor = "black"
        ),
        hovermode = "x",
        yaxis = list(title = "Count"),
        title = feature
      )
  } else {
    # if not interactive, return a ggplot with reformated xaxis
    p <- p +
      ggplot2::ggtitle(label = feature) +
      ggplot2::ylab('Counts') +
      ggplot2::xlab("") +
      ggplot2::coord_cartesian(xlim = c(0, x_cutoff)) +
      ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(p)
}


plot_qc <- function(obj,
                    wkdir,
                    stage  = "pre_filter",
                    title  = NULL,
                    feats  = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                    width  = 12,
                    height = 5,
                    dpi    = 200) {
  
  dir.create(wkdir, showWarnings = FALSE, recursive = TRUE)
  sample_name <- paste(unique(as.character(obj$orig.ident)), collapse = ", ")
  plot_title  <- if (is.null(title)) paste0(sample_name, " (", stage, ")") else title
  
  # static histograms
  p_static_list <- lapply(feats, function(f)
    single_histogram(obj, feature = f, interactive = FALSE))
  p_static <- patchwork::wrap_plots(p_static_list) +
    patchwork::plot_annotation(title = plot_title)
  ggplot2::ggsave(file.path(wkdir, paste0(stage, "_histograms.png")),
                  plot = p_static, width = width, height = height, dpi = dpi)
  
  # interactive histograms
  p_int_list <- lapply(feats, function(f)
    single_histogram(obj, feature = f, interactive = TRUE))
  
  axes <- list()
  for (i in seq_along(feats)) {
    suffix <- if (i == 1) "" else i
    xa <- p_int_list[[i]]$x$layoutAttrs[[1]][["xaxis"]]
    if (is.null(xa)) xa <- list()
    axes[[paste0("xaxis", suffix)]] <- modifyList(xa, list(title = feats[i]))
    axes[[paste0("yaxis", suffix)]] <- list(title = "Barcodes")
  }
  
  p_interactive <- do.call(
    plotly::layout,
    c(list(p          = plotly::subplot(p_int_list, nrows = 1),
           title      = list(text = plot_title),
           showlegend = FALSE,
           bargap     = 0,
           hovermode  = "x",
           hoverlabel = list(yanchor = "top"),
           margin     = list(t = 100)),
      axes)
  )
  htmlwidgets::saveWidget(p_interactive,
                          file.path(wkdir, paste0(stage, "_histograms.html")),
                          selfcontained = TRUE)
  
  invisible(list(histogram_static      = p_static,
                 histogram_interactive = p_interactive))
}
