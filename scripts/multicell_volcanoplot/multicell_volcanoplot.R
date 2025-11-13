#' @import ggplot2
#' @import ggpubr
#' @importFrom scales hue_pal
#'
NULL

#' Label Differential Expression Significance
#' 
#' Appends a column to a table of features tested for differential expression
#' that labels the features as significantly upregulated, significantly
#' downregulated, or insigificantly different.
#' 
#' @param df data frame of features tested for differential expression. Should include columns for measured effect sizes (ex. log2 fold-change) and p-values derived from statistical test of differentual expression.
#' @param lfc_column name of the column containing effect size estimates in df
#' @param pval_column name of the column containing p-values in df
#' @param log2FC_cutoff absolute threshold of effect size required to call a feature differentially expressed
#' @param pval_cutoff threshold of p-value required to call a feature differentially expressed
#' @return a character vector of significance labels 
#' @examples
#' test_df <- data.frame(
#'     "avg_log2_foldchange" = rnorm(n = 300),
#'     "p_val_adj" = runif(n = 300),
#'     "cell_type" = c(rep("A", 100), rep("B", 100), rep("C", 100)),
#'     "gene" = rep(paste0("gene", 1:100), 3)
#' )
#' test_df$significance <- label_significance(
#'      test_df, lfc_column = "avg_log2_foldchange", pval_column = "p_val_adj"
#' )
#' head(test_df)
#' 
label_significance <- function(
        df, 
        lfc_column, 
        pval_column,
        log2FC_cutoff = 0.5,
        pval_cutoff = 0.05
) {
    signif <- apply(
        df, 1, FUN = function(x) {
            lfc <- as.numeric(x[[lfc_column]])
            pval <- as.numeric(x[[pval_column]])
            if(is.na(lfc) | is.na(pval)){
                return(NA)
            }
            if(lfc > log2FC_cutoff & pval < pval_cutoff) {
                return("sig_UP")
            } else if(lfc < (-1 * log2FC_cutoff) & pval < pval_cutoff) {
                return("sig_DOWN")
            } else {
                return("insignificant")
            }
        }
    )
    return(signif)
}

#' Plot Differential Expression Volcano Plot for Multiple Cell Types
#' 
#' A plotting function for displaying differentially expressed features identified
#' in multiple cell types. Each tested feature is displayed as a unique point and
#' the cell type labels on the x axis obscure all features with effect sizes below
#' the stated threshold. Position on the y scale corresponds to estimated effect
#' size. Points are jittered on the x axis.
#' 
#' @param df data frame of features tested for differential expression. Should include columns for measured effect sizes (ex. log2 fold-change), p-values derived from statistical test of differentual expression, annotation of cell types, and annotation of differentially expressed features.
#' @param lfc_column name of the column containing effect size estimates in df.
#' @param pval_column name of the column containing p-values in df.
#' @param significance_column name of the column annotating statistically significant features.
#' @param cell_type_column name of the column annotating the cell type in which the feature was tested.
#' @param label_column name of the column containing feature names.
#' @param log2FC_cutoff absolute threshold of effect size required to call a feature differentially expressed.
#' @param pval_cutoff threshold of p-value required to call a feature differentially expressed.
#' @param significance_pal named vector of colors used to label points for features that are differentially upregulated, downregulated, or insignificantly different.
#' @param cell_type_pal named vector of colors used in blocks deliniating cell types. If NULL, [scales::hue_pal()] is used.
#' @return ggplot object displaying differentially expressed features as points grouped by cell type.
#' @examples
#' test_df <- data.frame(
#'     "avg_log2_foldchange" = rnorm(n = 300),
#'     "p_val_adj" = runif(n = 300),
#'     "cell_type" = c(rep("A", 100), rep("B", 100), rep("C", 100)),
#'     "gene" = rep(paste0("gene", 1:100), 3)
#' )
#' test_df$significance <- label_significance(
#'      test_df, lfc_column = "avg_log2_foldchange", pval_column = "p_val_adj"
#' )
#' multi_volcano(test_df)
#' 
#' multi_volcano(
#'     test_df, 
#'     label_column = "gene",
#'     pval_column = "p_val_adj", 
#'     lfc_column = "avg_log2_foldchange",
#'     significance_pal = c("sig_UP" = "#00FF00", "sig_DOWN" = "#000000", "insignificant" = "#BBBBBB"),
#'     cell_type_pal = c("A" = "#FF991C", "B" = "#7FFFD4", "C" = "#330066")
#' )
#' 
multi_volcano <- function(
        df, 
        lfc_column = "avg_log2_foldchange", 
        pval_column = "padj",
        significance_column = "significance", 
        cell_type_column = "cell_type",
        label_column = "gene",
        log2FC_cutoff = 0.5,
        pval_cutoff = 0.05,
        significance_pal = c("sig_UP" = "#FF0000", "sig_DOWN" = "#0000FF", "insignificant" = "#BBBBBB"),
        cell_type_pal = NULL
) {
    cts <- unique(df[[cell_type_column]])
    cts <- cts[order(cts)]
    if(is.null(cell_type_pal)) {
        cell_type_pal <- setNames(scales::hue_pal()(length(cts)), nm = cts)
    }
    
    # short dataframe to create celltype labels without redundant layers
    ls_ct <- list()
    ls_ct[[cell_type_column]] <- cts
    df_ct <- as.data.frame(ls_ct)
    
    # subset of data points to highlight with gene names
    df_sig <- df[
        df[[pval_column]] < pval_cutoff & 
            abs(df[[lfc_column]]) > log2FC_cutoff &
            !is.na(df[[pval_column]]),
    ]
    
    # gene points
    pl1 <- ggplot(df, aes(x = .data[[cell_type_column]], y = .data[[lfc_column]])) +
        geom_jitter(aes(color = .data[[significance_column]]), height = 0, width = 0.1) +
        scale_color_manual(values = significance_pal)
    
    # x-axis cell type labels
    pl2 <- pl1 +
        geom_tile(
            data = df_ct,
            aes(x = .data[[cell_type_column]], y = 0, fill = .data[[cell_type_column]]),
            color = "#000000", height = log2FC_cutoff * 2, show.legend = FALSE
        ) +
        geom_text(
            data = df_ct,
            aes(x = .data[[cell_type_column]], y = 0, label = .data[[cell_type_column]])
        ) +
        scale_fill_manual(values = cell_type_pal)
    
    # theme changes
    pl3 <- pl2 +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    
    # label significant points
    if(nrow(df_sig)){
        pl3 <- pl3 +
            ggrepel::geom_text_repel(
                data = df_sig,
                aes(
                    x = .data[[cell_type_column]],
                    y = .data[[lfc_column]], 
                    label = .data[[label_column]]
                ),
                max.overlaps = 50
            )
    }
    return(pl3)
}
