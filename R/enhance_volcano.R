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
                            palette = c('#66c2a5', 'grey50', '#fc8d62'),
                            cut_FC = "auto", cut_P = 0.05, top = 10,size = 2.0,expand=c(0.25,0.25),
                            genes_list = "top", highlight = NULL) {

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
  if (any(is.null(genes_list),is.null(highlight))) {
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
enhance_volcano_plot <- function(data,x,y,label = c("Down","Stable","Up"),cut_FC,cut_P, palette = c('#d53e4f', 'grey50', '#3288bd')) {

  names(palette) <- label

  if (length(cut_FC) == 1) {
    cut_FC <- c(-cut_FC, cut_FC)
  }

  FC_data <- data.frame(
    cut_FC = cut_FC
  )

  p <- ggplot(data=data, aes(x = get(x), y = -log10(get(y)), colour = group, alpha = group)) +
    geom_point(aes(color = group),size = 0.5) +
    geom_vline(data=FC_data, mapping=aes(xintercept=cut_FC),color = c(palette[label[1]],palette[label[3]]), linetype="longdash",size = 0.4) +
    geom_hline(aes(yintercept = -log10(cut_P)),color = palette[label[2]], linetype="longdash", size = 0.4) +
    scale_alpha_manual(values = c(1,0.5,1),
                       guide = "none")+
    scale_color_manual(
      values = palette,
      breaks = names(palette),
      labels = label,
      guide = guide_legend(override.aes = list(label = ""))
    ) +
    theme_nice() +
    labs(colour = "Group",x = "-log10(PValue)", y = "log2FoldChange",
         caption = paste(sprintf('FDR:  %.3f\n', cut_P),
                         sprintf('log 2 Fold Change: %.3f & %.3f \n',cut_FC,-cut_FC),
                         sprintf('Total: %1.0f variables',nrow(data))
         ))
  return(p)
}

#' Choose top genes from DEG data frame
#'
#' For example, if top = 50, head 50 and tail 50 will return from a
#' ordered DEG data frame after filtering by threshold value of P value
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param y which column is P value
#' @param top a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param cut_P a single number for threshold value of P value
#'
#' @importFrom utils head tail
#' @importFrom data.table data.table .SD
#'
#' @return a character vector of top genes
#' @export
#'
#' @examples
#' top_deg(DEG_df, x = "log2FoldChange", y = "pvalue", 50, 0.05)
top_deg <- function(deg_data, x, y, top, cut_P = 0.05) {

  data_f <- deg_data[which(deg_data[,y] <= cut_P),]

  d <- data.table(data_f, key=x, keep.rownames = TRUE)

  if (length(top) == 1) {
    td <- d[, head(.SD, top)]
    hd <-  d[, tail(.SD, top)]
  } else if (length(top) == 2) {
    td <- d[, head(.SD, top[2])]
    hd <-  d[, tail(.SD, top[1])] # the descending order
  } else {
    stop("top should be a single number or a length of 2 numeric vector")
  }

  hd_names <- hd$rn
  td_names <- td$rn

  choose_names <- c(hd_names, td_names)

  return(choose_names)
}

#' cut data
#'
#' group data for plot DEG volcano by FDR and logFC
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param y which column is P value
#' @param cut_FC a single number character or numeric vector in threshold value of log FC
#' @param cut_P threshold value of P value, can set for every cut_FC number in numeric vector format
#' @param label symbol word for groups
#' @param label_ns which group is the stable group
#'
#' @return add group column in data frame
#' @export
#'
#' @examples
#' cut_much(DEG_df, x = "log2FoldChange", y = "pvalue", cut_FC = 1,cut_P = 0.05,
#' label = c("Down","Stable","Up"), label_ns = "Stable")
cut_much <- function(deg_data,x, y,cut_FC = 1,cut_P = 0.05,label = c("Down","Stable","Up"), label_ns = "Stable") {

  if (length(cut_FC) == 1) {
    cut_FC <- c(-cut_FC, cut_FC)
  }
  if (length(cut_P) == 1) {
    cut_P <- rep(cut_P,length(cut_FC))
  }

  label_cg <- base::setdiff(label,label_ns)
  names(cut_P) <- label_cg
  deg_data$group <- cut(deg_data[,x],
                    breaks = c(-Inf,cut_FC,Inf),
                    labels = label)
  index = list()
  for (i in label_cg) {
    index[[i]] <- setdiff(which(deg_data$group == i), which(deg_data$group == i & deg_data[,y] < cut_P[i]))
    deg_data$group[index[[i]]] <- label_ns
  }
  message(paste(label, as.character(table(deg_data$group)),"\n"))
  return(deg_data)
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

  max_fdr <- max(-log10(data[,y]))

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
        bg.color = "#e0e0e0", # shadow color
        bg.r = 0.15          # shadow radius
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
        bg.color = "#e0e0e0", # shadow color
        bg.r = 0.15          # shadow radius
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
          bg.color = "#26b1fb", # shadow color
          bg.r = 0.15          # shadow radius
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
          bg.color = "#26b1fb", # shadow color
          bg.r = 0.15          # shadow radius
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
                                  plot.caption = element_text(family = "Times", size = 6, face = "italic", colour = "dodgerblue"))}

