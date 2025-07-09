library(tidyverse)
library(edgy)

## NCI Example ----
df <-
  read_tsv('inst/Sample.Data.txt', col_names = FALSE) |>
  filter(X1 == 0L)

res <- trend_reg(log(X3) ~ X2, df,
                 opts = knot_opts(max_knot = 5, min_obs_between = 2))

res

df |>
  ggplot() +
  geom_point(aes(X2, X3)) +
  geom_line(data = res$data, aes(X2, fit)) +
  theme_minimal(14)

## Different time scale ----

library(nycflights13)

df <-
  flights |>
  summarise(
    delay = mean(arr_delay, na.rm = TRUE),
    .by = c(month, day)
  ) |>
  arrange(month, day) |>
  mutate(
    id = consecutive_id(month, day)
  )

res <- trend_reg(delay ~ id, df, opts = knot_opts(max_knot = 1,
                                                          min_obs_end = 4,
                                                          min_obs_between = 6))

df |>
  ggplot() +
  geom_point(aes(id, delay)) +
  geom_line(data = res$data, aes(id, fit)) +
  theme_minimal(14)
