convertToNumeric <- function(expression_df) {
  dt <- as.data.table(expression_df)
  df <- as.data.frame(dt)
  df[] <- lapply(df, type.convert, as.is = TRUE)
  rownames(df) <- rownames(expression_df)
  colnames(df) <- colnames(expression_df)
  return (df)
}

convertToNumeric(bulk_counts_t)
