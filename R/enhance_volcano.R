#' enhance_volcano
#'
#' an enhance_colcano function for DEG volcano plot
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param y which column is P value
#' @param label symbol word for groups,default is `c("Down","Stable","Up")`
#' @param label_ns which group is the stable group
#' @param palette color for every label
#' @param cut_FC a numeric vector in threshold value of log FC
#' @param cut_P threshold value of P value, can set for every cut_FC number in numeric vector format
#' @param top a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param size a single number of font size
#' @param expand if labels in plot display wrong, try increase it.
#' @param genes_list a gene names character vectors;
#' @param highlight a gene names character vectors; default NULL
#'
#' @importFrom stats sd
#'
#' @return a ggplot ob
#' @export
#'
#' @examples
#' enhance_volcano(DEG_df,x = 'log2FoldChange', y = 'pvalue',
#'                 label = c("Down","Stable","Up"), label_ns = "Stable",
#'                 palette = c('#66c2a5', 'grey50', '#fc8d62'),
#'                 cut_FC = "auto",cut_P = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
#'                 genes_list = "top", highlight = NULL)
enhance_volcano <- function(deg_data,x, y,
                            label = c("Down","Stable","Up"), label_ns = "Stable",
                            palette =  c("#2874C5", "grey", "#f87669"),
                            cut_FC = "auto", cut_P = 0.05, top = 10,size = 2.0,expand=c(0.25,0.25),
                            genes_list = "top", highlight = NULL) {

  if (!x %in% colnames(deg_data)) stop("x must be a character belongs to column names of deg_data")

  if (!y %in% colnames(deg_data)) stop("y must be a character belongs to column names of deg_data")

  if (!label_ns %in% label) stop("label_ns must be one of label")

  if (cut_FC == "auto") {
    cut_FC = mean(abs(deg_data[,x])) + 2*sd(abs(deg_data[,x]))
  }

  if (genes_list == "top") {
    genes_list = top_deg(deg_data, x = x, y = y, top, cut_P)
  }

  deg_data_g <- cut_much(deg_data, x= x,y = y, cut_FC = cut_FC, cut_P = cut_P,
                         label = label, label_ns = label_ns)

  plot <- enhance_volcano_plot(deg_data_g,x = x, y = y, label = label, cut_FC = cut_FC, cut_P = cut_P)
  #
  if (any(!is.null(genes_list),!is.null(highlight))) {
    plot <- show_genes(deg_data_g, x = x, y = y, genes_list, highlight,plot,size,expand)
  }
  #
  return(plot)

}

#' enhance DEG volcano basic plot
#'
#' more beautiful and clear DEG volcano plot
#'
#' @param data a grouped DEG data frame
#' @param x column name of log Fold Change
#' @param y column name of p value
#' @param label symbol word for groups,default is `c("Down","Stable","Up")`
#' @param cut_FC a numeric vector in threshold value of log FC
#' @param cut_P threshold value of P value, can set for every cut_FC number in numeric vector format
#' @param palette color for every label
#'
#' @import ggplot2
#'
#' @return a ggplot object
#'
#' @noRd
#' @examples
#' enhance_volcano_plot(data,x = 'log2FoldChange', y = 'pvalue', cut_FC = 1, cut_P = 0.05)
enhance_volcano_plot <- function(data,x,y,label = c("Down","Stable","Up"),cut_FC,cut_P,
                                 palette =  c("#2874C5", "grey", "#f87669")) {

  names(palette) <- label

  if (length(cut_FC) == 1) {
    cut_FC <- c(-cut_FC, cut_FC)
  }

  FC_data <- data.frame(
    cut_FC = cut_FC
  )

  count_updown <- table(data$group)

  p <- ggplot(data=data, aes(x = get(x), y = -log10(get(y)), colour = group, alpha = group)) +
    geom_point(aes(color = group),size = 0.6) +
    geom_vline(data=FC_data, mapping=aes(xintercept=cut_FC),color = c(palette[label[1]],palette[label[3]]), linetype="longdash",size = 0.4) +
    geom_hline(aes(yintercept = -log10(cut_P)),color = palette[label[2]], linetype="longdash", size = 0.4) +
    scale_alpha_manual(values = c(0.6,0.5,0.6),
                       guide = "none")+
    scale_color_manual(
      values = palette,
      breaks = names(palette),
      labels = label,
      guide = guide_legend(override.aes = list(label = ""))
    ) +
    theme_nice() +
    labs(colour = "Group",x = "log2FoldChange", y = "-log10(PValue)",
         subtitle = paste(sprintf('%s:  %.3f;',y, cut_P),
                         sprintf('FC: %.3f;',cut_FC),
                         sprintf('Up: %1.0f; Down: %1.0f;',count_updown[3],count_updown[1]),
                         sprintf('Total: %1.0f',nrow(data))
         ))
  return(p)
}

#' Showing gene labels and highlighting gene points in a volcano plot in ggplot
#'
#' specific some genes to label
#'
#' @param data a grouped DEG data frame
#' @param genes_list a gene names character vectors;
#' @param highlight a gene names character vectors; default NULL
#' @param p volcano plot in enhance_volcano_plot
#'
#' @importFrom data.table data.table
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @return a ggplot ob
#'
#' @noRd
show_genes <- function(data, x = "log2FoldChange", y = "pvalue", genes_list, highlight = NULL,p,size = 1.8, expand = c(0.12, 0.12)) {

  max_lfc <- max(data[,x])
  min_lfc <- min(data[,x])

  max_fdr <- max_fdr <- max(na.omit(-log10(data[,y])))

  if (max_fdr) {
    max_fdr <- 300
  }

  data <- data.table(data, key=x, keep.rownames = TRUE)

  if (any(!is.null(genes_list),!is.null(highlight))) {

    if (!is.null(genes_list)) {
      genes_list <- data.frame(
        rn = genes_list
      )
      # look up piont for labels
      label_data <- merge(genes_list,data)
    }

    if (!is.null(highlight)) {
      highlight <- data.frame(
        rn = highlight
      )
      hlight_data <- merge(highlight,data)
    }

    res <- p +
      # scale_y_continuous(expand = expand) +
      scale_x_continuous(expand = expand) +
      geom_text_repel(
        size = size,
        data = label_data[which(label_data[,x] > 0),],
        aes(label = rn),
        nudge_y      = -1,
        direction    = "y",
        hjust        = 0,
        segment.size = 0.1,
        segment.linetype = 6,
        max.overlaps = 10,
        max.iter = 1000000,
        max.time = 10,
        nudge_x =  label_data[which(label_data[,x] > 0),][,x] + max_lfc+max_lfc/4,
        min.segment.length = 0,
        fontface = "bold",family = "Times"
        # bg.color = "#e0e0e0", # shadow color
        # bg.r = 0.15          # shadow radius
      ) +
      geom_text_repel(
        size = size,
        data = label_data[which(label_data[,x] < 0),],
        aes(label = rn),
        nudge_y      = -1,
        direction    = "y",
        hjust        = 1,
        segment.size = 0.1,
        segment.linetype = 6,
        max.overlaps = 10,
        max.iter = 1000000,
        max.time = 10,
        nudge_x =  label_data[which(label_data[,x] < 0),][,x] + min_lfc+min_lfc/4,
        min.segment.length = 0,
        fontface = "bold",family = "Times"
        # bg.color = "#e0e0e0", # shadow color
        # bg.r = 0.15          # shadow radius
      )

    if (!is.null(highlight)) {

      res <- res +
        geom_point(size = 1.5,data = hlight_data,shape = 1,stroke = 0.35,
                   color = "#1474ff",fill = NA) +
        geom_text_repel(
          size = size,
          data = hlight_data[which(hlight_data[,x] > 0),],
          aes(label = rn),
          nudge_y      = -1,
          direction    = "y",
          hjust        = 0,
          segment.size = 0.1,
          segment.linetype = 6,
          max.overlaps = 1,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,x] > 0),][,x] + max_lfc+max_lfc/4,
          min.segment.length = 0,
          fontface = "bold",family = "Times"
          # bg.color = "#26b1fb", # shadow color
          # bg.r = 0.15          # shadow radius
        ) +
        geom_text_repel(
          size = size,
          data = hlight_data[which(hlight_data[,x] < 0),],
          aes(label = rn),
          nudge_y      = -1,
          direction    = "y",
          hjust        = 1,
          segment.size = 0.1,
          segment.linetype = 6,
          max.overlaps = 1,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,x] > 0),][,x] + min_lfc+min_lfc/4,
          min.segment.length = 0,
          fontface = "bold",family = "Times"
          # bg.color = "#26b1fb", # shadow color
          # bg.r = 0.15          # shadow radius
        )

    }


    return(res)

  }

}

#' theme_nice
#'
#' a nice theme for DEG volcano
#'
#' @param ... param passed from theme
#'
#' @importFrom ggplot2 theme element_line element_text element_rect
#'
#' @return a ggplot theme
#'
#' @noRd
theme_nice <- function(...) {theme(...,
                                  axis.line = element_line(size = 0.2, linetype = "solid"),
                                  axis.ticks = element_line(size = 0.2),
                                  panel.grid.major = element_line(colour = "gray30", size = 0.1,linetype = "dotted"),
                                  panel.grid.minor = element_line(linetype = "blank"),
                                  axis.title = element_text(family = "Times"),
                                  axis.text = element_text(family = "Times"),
                                  axis.text.x = element_text(family = "Times"),
                                  axis.text.y = element_text(family = "Times"),
                                  legend.text = element_text(family = "Times"),
                                  legend.title = element_text(family = "Times"),
                                  panel.background = element_rect(fill = NA),
                                  legend.key = element_rect(fill = NA),
                                  legend.background = element_rect(fill = NA),
                                  legend.position = "top", legend.direction = "horizontal",
                                  plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 6, face = "italic", colour = "black"))}

