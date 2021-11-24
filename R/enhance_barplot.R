# https://community.rstudio.com/t/how-to-start-text-label-from-extreme-left-in-bar-plot-using-geom-col-and-geom-text-in-ggplot2/77256
#' Enhanced barplot
#'
#' plot enrich result in bar plot
#'
#' @param data 'enrichResult' object, enrichGO result
#' @param fillstrip fill strip rect or not
#' @param split_color color for CC, MF, and BP
#' @param bar_color color for bar
#' @param print logic for print plot
#' @param showCategory Category numbers to show
#' @param by one of Count and GeneRatio
#' @param order logical
#' @param drop logical
#' @param split separate result by 'split' variable
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggfun element_roundrect
#' @importFrom grid grid.draw
#'
#' @return ggplot or grid plot
#' @export
#'
#' @examples
#' \dontrun{
#' test <- enrich_go(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue")
#' enhance_barplot(test$Down,showCategory=10,split = "ONTOLOGY")
#' enhance_barplot(test$Down,showCategory=30)
#' }
enhance_barplot <- function(data, showCategory = 10,by = "Count",split=NULL,
                            order = FALSE,
                            drop = FALSE,
                            fillstrip = TRUE,
                            print = FALSE,
                            split_color = RColorBrewer::brewer.pal(3,"Dark2"),
                            bar_color = rev(RColorBrewer::brewer.pal(5,"GnBu")[3:5])
                            ) {

  if (missing(split)){
    dat <- fortify.enrichResult(model = data, showCategory = showCategory, by = by, order = order, drop = drop) %>%
      sort_goTerms(by=by,split = split)

    p <- barplot_base2(dat, bar_color = bar_color, text_color = split_color[[1]])

    return(p)

  } else {
    dat <- fortify.enrichResult(model = data, showCategory = showCategory, by = by, order = order, drop = drop, split = split) %>%
      sort_goTerms(by=by,split = split)

    p <- barplot_base(dat, bar_color = bar_color, split_color = split_color)
    p <- p +
      facet_grid(ONTOLOGY~., scales="free", space="free_y", switch = "both") +
      theme(strip.background=element_roundrect(fill=NA, color=NA, r=0.31415,size = 0.5,linetype = "dotted")
            # ,legend.background = element_roundrect(fill=NA, color="grey80", r=0.31415,size = 0.5,linetype = "solid")
            )

    g <- ggplot_gtable(ggplot_build(p))
    strip_both <- which(grepl('strip-', g$layout$name))

    k <- 1

    if (fillstrip) {
      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lty <- "solid"
        g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- "white"
        k <- k+1
      }
    } else {
      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- NA
        g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- split_color[k]
        k <- k+1
      }
    }

    if (print) {
      grid.draw(g)
    }

    return(g)

    }

}

#' barplot_base
#'
#' plot basic barplot for enhance barplot
#'
#' @param data a data frame from enrichGO result
#' @param bar_color color for bar
#' @param split_color color for ont
#'
#' @import ggplot2
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
barplot_base <- function(data, bar_color, split_color) {
  p <- ggplot(data,aes(Count, myY, fill = p.adjust)) +
    geom_shadowtext(aes(color = ONTOLOGY, label = Description),
                    x = 0,
                    nudge_y = -0.5,
                    hjust = -0.01, show.legend = F,
                    # position=position_dodge2(width=0.9),
                    size = 4,
                    bg.colour='#ffffff',bg.r = NA) +
    scale_color_manual(values = split_color)+
    new_scale_color() +
    geom_segment(aes(x=0,y=myY-1,xend = Count,yend = myY-1,color = p.adjust),size = 1, lineend = 'round') +
    scale_y_continuous(name = "Description", breaks = data$myY, labels = data$Description, expand = c(0, 0.5)) +
    # geom_col(width = .15,alpha = 1,aes(color = p.adjust),
    #          position=position_dodge2(padding = 0.9)) +
    scale_color_gradientn(colours = bar_color) +
    # geom_vline(aes(xintercept = 0,colour = ONTOLOGY),size = 0.5,show.legend = F)+
    # geom_text(aes(label = Description,color = ONTOLOGY),
    #           x = 0,
    #           hjust = 0,
    #           size = 4) +
    scale_x_continuous(expand = c(0, 0))+
    theme_minimal()+
    theme(legend.position = "right",
          # plot.margin=margin(t= 100,b=2,l=-2,r= 2,unit="pt"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dotted"),
          panel.grid.minor = element_blank(),
          # axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

    return(p)
}

#' barplot_base2
#'
#' plot basic barplot for enhance barplot not split
#'
#' @param data a data frame from enrichGO result
#' @param bar_color color for bar
#' @param text_color color for text
#'
#' @import ggplot2
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
barplot_base2 <- function(data, bar_color, text_color) {
  p <- ggplot(data,aes(Count, myY, fill = p.adjust)) +
    geom_vline(aes(xintercept = 0,colour = text_color),size = 0.5,show.legend = F)+
    geom_shadowtext(aes(color = text_color, label = Description),
                    x = 0,
                    nudge_y = -0.5,
                    hjust = -0.01, show.legend = F,
                    # position=position_dodge2(width=0.9),
                    size = 4,
                    bg.colour='#ffffff',bg.r = NA) +
    scale_color_manual(values = text_color)+
    new_scale_color() +
    geom_segment(aes(x=0,y=myY-1,xend = Count,yend = myY-1,color = p.adjust),size = 1, lineend = 'round') +
    scale_y_continuous(name = "Description", breaks = data$myY, labels = data$Description, expand = c(0, 0.5)) +
    # geom_col(width = .15,alpha = 1,aes(color = p.adjust),
    #          position=position_dodge2(padding = 0.9)) +
    scale_color_gradientn(colours = bar_color) +
    # geom_text(aes(label = Description,color = ONTOLOGY),
    #           x = 0,
    #           hjust = 0,
    #           size = 4) +
    scale_x_continuous(expand = c(0, 0))+
    theme_minimal()+
    theme(legend.position = "right",
          legend.background = element_roundrect(fill=NA, color="grey80", r=0.31415,size = 0.5,linetype = "solid"),
          # plot.margin=margin(t= 100,b=2,l=-2,r= 2,unit="pt"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dotted"),
          panel.grid.minor = element_blank(),
          # axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

  return(p)
}

#' KEGG bar plot
#'
#' plot for enrich_kegg result
#'
#' @param data a enrich_kegg result, list format contains Up and Down
#' @param top top rows of up and down
#' @param down_label which one is Down
#'
#' @importFrom stringr str_wrap
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @importFrom dplyr top_n mutate
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return a ggplot OB
#' @export
#'
#' @examples
#' \dontrun{
#' test <- enrich_kegg(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue")
#' kegg_barplot(test)
#' }
kegg_barplot <- function(data,top = -10,down_label = "Down") {

  lis <- lapply(seq_along(data), function(x)

    if (names(data)[[x]] == down_label) {
      data[[x]]@result$upordown <- -1
      data[[x]]@result
    }
    else {
      data[[x]]@result$upordown <- 1
      data[[x]]@result
    }
  )

  names(lis) <- names(data)

  lis_f <- lapply(lis, function(x)
    # x[order(x[,"pvalue"])[1:top],]
    x %>% top_n(top,wt = pvalue)
  )

  dat <- do.call(rbind,lis_f)

  dat$fulog10pvalue <- (-log10(dat$pvalue))*dat$upordown

  dat <- dat %>%
    mutate(Description = fct_reorder(Description,fulog10pvalue))

  p <- ggplot(dat, aes(fulog10pvalue, Description, fill = as.character(upordown))) +
    # scale_color_manual(values = c("#e66101","#5e3c99"),name = "",
    #                    breaks = c("1","-1"),
    #                    labels= c("Up","Down")) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#3685af","#af4236"),guide = "none") +
    # scale_alpha_continuous(range = c(0.3, 1)) +
    geom_text(aes(label = Description),data = dat[which(dat["fulog10pvalue"] >= 0),],
              x = 0,
              hjust = 0,
              size = 3) +
    geom_text(aes(label = Description),data = dat[which(dat["fulog10pvalue"] < 0),],
              x = 0,
              hjust = 1,
              size = 3) +
    labs(y = "Description", x = "-log10(PValue)")+
    theme_minimal()+
    theme(legend.position = "none",
          # plot.margin=margin(b=200,l=-2,unit="pt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

  return(p)

}

#' gseKEGG bar plot
#'
#' plot for enrich_gsekegg result
#'
#' @param data a enrich_gsekegg result
#' @param pvalue_cut filter by pvalue
#' @param enrichmentScore_cut filter by enrichmentScore
#' @param top top rows of up and down
#'
#' @import ggplot2
#' @importFrom dplyr top_n mutate
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
#' @export
#'
#' @examples
#' \dontrun{
#' test <- enrich_gsekegg(DEG_df,x = "log2FoldChange")
#' geskegg_barplot(test)
#' }
geskegg_barplot <- function(data, pvalue_cut = 0.1, enrichmentScore_cut = 0.5, top = -10) {
  down_kegg<-data[data$pvalue< pvalue_cut & data$enrichmentScore < -enrichmentScore_cut,];down_kegg$group=-1
  up_kegg<-data[data$pvalue< pvalue_cut & data$enrichmentScore > enrichmentScore_cut,];up_kegg$group=1

  # down_kegg <- down_kegg[order(down_kegg[,"pvalue"])[1:top],]
  # up_kegg <- up_kegg[order(up_kegg[,"pvalue"])[1:top],]
  down_kegg <- down_kegg %>% top_n(top,wt = pvalue)
  up_kegg <- up_kegg %>% top_n(top,wt = pvalue)

  dat=rbind(up_kegg,down_kegg)
  dat$fulog10pvalue = -log10(dat$pvalue)
  dat$fulog10pvalue = dat$fulog10pvalue*dat$group

  dat <- dat %>%
    mutate(Description = fct_reorder(Description,fulog10pvalue))

  ggplot(dat, aes(fulog10pvalue, Description, fill = as.character(group))) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#3685af","#af4236"),guide = "none") +
    # scale_alpha_continuous(range = c(0.3, 1)) +
    geom_text(aes(label = Description),data = dat[which(dat["fulog10pvalue"] >= 0),],
              x = 0,
              hjust = 0,
              size = 3) +
    geom_text(aes(label = Description),data = dat[which(dat["fulog10pvalue"] < 0),],
              x = 0,
              hjust = 1,
              size = 3) +
    labs(y = "Description", x = "-log10(PValue)")+
    theme_minimal()+
    theme(legend.position = "none",
          # plot.margin=margin(b=200,l=-2,unit="pt"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

}

fortify.enrichResult <- enrichplot:::fortify.enrichResult

sort_goTerms <- function(data,split,by) {
  if (is.null(split)) {
    data <- data[
      order( data[,by], decreasing = FALSE),
    ]
  } else {
    data <- data[
      order( data[,split], data[,by], decreasing = FALSE),
    ]
  }
  data$Description <- factor(data$Description,levels = data$Description)
  data$myY <- as.numeric(data$Description)
  return(data)
}
