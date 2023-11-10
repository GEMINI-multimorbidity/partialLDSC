###### LDSCpartial - other functions ######


#' Get list of conditions from partial_ldsc results
#' @param x an object containing the output from a call to \code{\link{partial_ldsc}()}
#' @param ... further arguments
#' 
#' @return the names of all the conditions from the analysis
#' @export
get_conditions <- function(x,...) {
  if(length(x) != 10) stop("only works for objects containing the output from a call to partial_ldsc()")
  
  if(!all(names(x) == c("res_diff", "S", "V", "S_Stand", "V_Stand",      
                    "partial.S", "partial.V", "partial.S_Stand",
                    "partial.V_Stand", "I"))) stop("only works for objects containing the output from a call to partial_ldsc()")
  conditions = colnames(x$partial.S)
  return(conditions)
}


#' Get list of pairs from partial_ldsc results
#' @param x an object containing the output from a call to \code{\link{partial_ldsc}()}
#' @param ... further arguments
#' 
#' @return the names of all the pairs from the analysis
#' @export
get_pairs <- function(x,...) {
  if(length(x) != 10) stop("only works for objects containing the output from a call to partial_ldsc()")
  
  if(!all(names(x) == c("res_diff", "S", "V", "S_Stand", "V_Stand",      
                        "partial.S", "partial.V", "partial.S_Stand",
                        "partial.V_Stand", "I"))) stop("only works for objects containing the output from a call to partial_ldsc()")
  pairs = x$res_diff %>% dplyr::mutate(pair = paste0(.data$condition.1, "-", .data$condition.2)) %>% dplyr::pull(.data$pair)
  return(pairs)
}



#' Creates a forest plot to visualise partial genetic correlation results
#'
#'
#' @param obj an object containing the output from a call to \code{\link{partial_ldsc}()}
#' @param save_file Should the graphic should be saved?
#'        \code{default=FALSE}, graphic will be displayed on the on-screen device (logical)
#' @param file_name The name of the file to be saved (if \code{save_file} is \code{TRUE}),
#'        \code{default=NULL}, will use <forest_plot.jpeg>, extensions allowed are ".jpeg"/".pdf"/".tiff" (character)
#' @param file_width The width of the figure to be saved (if \code{save_file} is \code{TRUE}),
#'        to help with readability if the conditions' names are too long, \code{default=NULL} (numeric)
#' @param conditions The names of the conditions to be included, \code{default=NULL}, 
#'         will use all conditions and all corresponding pairs, a maximum of 8 conditions (character)
#' @param pairs The names of the pairs to be included, \code{default=NULL}, 
#'         if provided, will overrule \code{conditions}, a maximum of 30 pairs (character)
#' @param order How to order the conditions,  possibility to use the difference (\code{"difference"}),
#'        the unadjusted genetic correlation (\code{"unadjusted"}), or the partial genetic correlation (\code{"partial"}),
#'        \code{default="difference"} (character)
#' 
#' 
#' @return a Forest plot
#'
#' @export

forest_plot <- function(obj, save_file = F, file_name = NULL, file_width = NULL,
                       conditions = NULL, pairs = NULL, order = "difference") {
  
  ## check parameters
  if(length(obj) != 10) stop("only works for objects containing the output from a call to partial_ldsc()")
  
  if(!all(names(obj) == c("res_diff", "S", "V", "S_Stand", "V_Stand",      
                        "partial.S", "partial.V", "partial.S_Stand",
                        "partial.V_Stand", "I"))) stop("only works for objects containing the output from a call to partial_ldsc()")
  
  if(!is.logical(save_file)) stop("save_file : should be logical")
  # if no name, use the one from the analysis (in log file)
  if(save_file && is.null(file_name)){
    file_name = "forest_plot.jpeg"
  }
  if(save_file && !is.character(file_name)) stop("file_name : should be a character")
  if(save_file){
    extension = dplyr::last(as.vector(stringr::str_split(file_name, stringr::fixed("."), simplify = T)))
    if(!extension %in% c("jpeg", "tiff", "pdf")) stop("file extension should be either \".jpeg\", \".pdf\" or \".tiff\"")
    
  }
  if(save_file && !is.numeric(file_width)) stop("file_width : should be numeric")
  if(save_file && !is.null(file_width)) file_width=15
  
  if(!is.null(conditions) && !is.character(conditions)) stop("conditions : should be character")
  if(!is.null(conditions) && length(conditions)>8) stop("conditions : can not use more that 8 conditions")
  if(!is.null(conditions) && !any(conditions %in% partialLDSC::get_conditions(obj))) stop("conditions : some conditions are incorrect, use get_conditions() to get the list of conditions that can be used")
  if(!is.null(pairs) && !is.character(pairs)) stop("pairs : should be character")
  if(!is.null(pairs) && length(pairs)>30) stop("pairs : can not use more that 30 pairs")
  if(!is.null(pairs) && !any(pairs %in% partialLDSC::get_pairs(obj))) stop("pairs : some pairs are incorrect, use get_pairs() to get the list of pairs that can be used")
  
  if(!order %in% c("difference", "unadjusted", "partial")) stop("order : should be either \"difference\", \"unadjusted\" or \"partial\"")
  
  # subset data if needed
  res_rg = obj$res_diff
  
  if(!is.null(pairs)){
    res_rg %>%
      dplyr::mutate(pair = paste0(.data$condition.1, "-", .data$condition.2)) %>%
      dplyr::filter(.data$pair %in% pairs) %>% 
      dplyr::mutate(pair = NULL) -> res_rg
  } else if(!is.null(conditions)) {
    res_rg %>%
      dplyr::filter(.data$condition.1 %in% conditions, 
                    .data$condition.2 %in%  conditions) -> res_rg
  }
    
  # order
  if(order == "difference") {
    res_rg %>%
      dplyr::arrange(.data$diff.P) %>% 
      dplyr::select(.data$condition.1, .data$condition.2, .data$rg ,.data$rg.SE ,
                    .data$partial_rg, .data$partial_rg.SE, .data$diff.P) -> res_plot
  } else if(order == "unadjusted") {
    res_rg %>%
      dplyr::arrange(dplyr::desc(.data$rg)) %>% 
      dplyr::select(.data$condition.1, .data$condition.2, .data$rg ,.data$rg.SE ,
                    .data$partial_rg, .data$partial_rg.SE, .data$diff.P) -> res_plot
  } else if(order == "difference") {
    res_rg %>%
      dplyr::arrange(dplyr::desc(.data$partial.rg)) %>% 
      dplyr::select(.data$condition.1, .data$condition.2, .data$rg ,.data$rg.SE ,
                    .data$partial_rg, .data$partial_rg.SE, .data$diff.P) -> res_plot
    
  }
  
  res_plot %>%
    dplyr::mutate(pair = paste0(.data$condition.1, " ~ ", .data$condition.2),
                  label = paste0("(", format(.data$diff.P, scientific=T, digits=2), ")  ")) -> res_plot
  
  # pivot
  res_plot %>%
    dplyr::select(.data$pair, .data$label, .data$rg, .data$partial_rg ) %>%
    tidyr::pivot_longer(!c(.data$pair, .data$label), names_to = "adjustment", values_to = "rg") %>%
    dplyr::mutate(adjustment = dplyr::case_when(
      adjustment == "rg" ~ "unadjusted",
      adjustment == "partial_rg" ~ "adjusted for BMI genetics")) -> plot_rg
  res_plot %>%
    dplyr::transmute(.data$pair,.data$label, .data$rg.SE, .data$partial_rg.SE ) %>%
    tidyr::pivot_longer(!c(.data$pair,.data$label), names_to = "adjustment", values_to = "SE") %>%
    dplyr::mutate(adjustment = dplyr::case_when(
      adjustment == "rg.SE" ~ "unadjusted",
      adjustment == "partial_rg.SE" ~ "adjusted for BMI genetics"))-> plot_SE
  # they need to have the same adjustment name
  
  plot_all = dplyr::full_join(plot_rg,plot_SE, by=c("pair", "label", "adjustment"))
  plot_all$adjustment = forcats::as_factor(plot_all$adjustment)
  plot_all$pair = forcats::fct_rev(forcats::as_factor(plot_all$pair))
  
  
  ggplot2::theme_set(ggplot2::theme_bw() + 
                       ggplot2::theme(axis.text = ggplot2::element_text(color="black"),
                                      axis.title = ggplot2::element_text(color="black"),
                                      axis.ticks = ggplot2::element_line(color="black"),
                                      axis.line = ggplot2::element_line(color="black"),
                                      legend.text = ggplot2::element_text(color="black")))
  
  ggplot2::ggplot(plot_all, ggplot2::aes(x=.data$pair, y=.data$rg, ymin=.data$rg-1.96*.data$SE, ymax=.data$rg+1.96*.data$SE, color=.data$adjustment, shape=.data$adjustment, label=.data$label)) +
    ggplot2::geom_pointrange( position=ggplot2::position_dodge(width=0.2)) +
    ggplot2::geom_hline(yintercept = 0, col="black", linetype="longdash") +
    ggplot2::ylim(min(plot_all$rg-1.96*plot_all$SE), max(plot_all$rg+1.96*plot_all$SE)+0.1)+
    # q-value attenuation
    ggplot2::geom_text(vjust=-0.5, hjust=0.5, mapping = ggplot2::aes(y=(max(.data$rg+1.96*.data$SE)+0.1), 
                                                   label=.data$label), color="black", size=2)+
    ggplot2::labs(x="", 
         y="genetic correlation (95% CI)",
         color="", shape="") +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_color_manual(values=c("#5270D4", "#FBCF8F")) +
    ggplot2::coord_flip() -> figure
  
  if(save_file){
    ggplot2::ggsave(file_name, 
                    figure, height = (4+1.5*nrow(res_rg)), width = file_width, units = "cm", dpi = 320)
  }
  
  return(figure)
}
