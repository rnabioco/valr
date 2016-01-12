library(dplyr)
library(Rbedtools)

# bed1_df <- data.frame(
#   start = c(100, 150, 600),
#   end   = c(200, 250, 800)
# )
# bed2_df <- data.frame(
#   start = c(175, 175, 799),
#   end   = c(200, 225, 900)
# )
bed1_df <- tibble(
  ~start,    ~end,
  100,       200,
  150,       250,
  600,       800
)
bed2_df <- tibble(
  ~start,    ~end,
  175,       200,
  175,       225,
  799,       900
)

intersect_cpp(bed1_df, bed2_df)
