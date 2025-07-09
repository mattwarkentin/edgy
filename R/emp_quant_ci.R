emp_quant_ci <- function(model, y, b = 100) {
  # 1. Get residuals from the model fit
  y_res <- residuals(model)

  # 2. Get B samples of size N from runif(0, 1)
  bs <- purrr::map(seq_len(b), \(b) runif(length(y)))

  # 3. For each B, generate inverse CDF for residuals
  y_res_new <- purrr::map(bs, \(x) {

  })

  # 4. For each B, gene

  # 5. For each B,

  # 6. Construct empirical CDF CI
}
'inverted_ecdf <- quantile(dat, probs = my_ecdf(dat), type = 3)'
