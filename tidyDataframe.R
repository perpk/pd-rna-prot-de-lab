tidyDataframe <- function(expression_df) {
  expression_df %>% remove_rownames %>% column_to_rownames(var = 'GeneID') %>% t -> expression_df_t
  expression_df_t %>% rownames %>% str_sub(0, -30) %>% str_replace("PPMI.Phase\\d{1}.IR\\d{1}.", "") -> rownames(expression_df_t)
  numvals <- as.numeric(sub(".*\\.", "", rownames(expression_df_t)))
  return (expression_df_t[order(numvals), ])
}