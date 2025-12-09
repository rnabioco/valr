# Quick benchmark script for comparing branches
# Run this on both main and cpp11-full-migration branches

library(valr)
library(bench)

genome <- read_genome(valr_example('hg19.chrom.sizes.gz'))

# Number of intervals - start small to see the effect
n_sizes <- c(1e4, 1e5, 5e5)

cat("Quick benchmark for key interval tree operations\n")
cat("=================================================\n\n")

for (n in n_sizes) {
  cat(sprintf("\n--- %s intervals ---\n", format(n, big.mark = ",")))

  x <- bed_random(genome, n = n, seed = 1010486)
  y <- bed_random(genome, n = n, seed = 9283019)

  res <- mark(
    bed_intersect(x, y),
    bed_closest(x, y),
    bed_subtract(x, y),
    min_time = 2,
    iterations = 3,
    check = FALSE,
    time_unit = 's'
  )

  print(res[, c("expression", "min", "median", "mem_alloc", "n_itr")])
}
