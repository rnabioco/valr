x <- tibble::tribble(
  ~chrom ,
  ~start ,
  ~end   ,
  'chr1' ,
     100 ,
     250 ,
  'chr2' ,
     250 ,
     500
)

y <- tibble::tribble(
  ~chrom ,
  ~start ,
  ~end   ,
  ~value ,
  'chr1' ,
     100 ,
     250 ,
      10 ,
  'chr1' ,
     150 ,
     250 ,
      20 ,
  'chr2' ,
     250 ,
     500 ,
     500
)

bed_glyph(bed_map(x, y, value = sum(value)), label = 'value')

# summary examples
bed_map(x, y, .sum = sum(value))

bed_map(x, y, .min = min(value), .max = max(value))

# identify non-intersecting intervals to include in the result
res <- bed_map(x, y, .sum = sum(value))
x_not <- bed_intersect(x, y, invert = TRUE)
dplyr::bind_rows(res, x_not)

# create a list-column
bed_map(x, y, .values = list(value))

# use `nth` family from dplyr
bed_map(x, y, .first = dplyr::first(value))

bed_map(x, y, .absmax = abs(max(value)))

bed_map(x, y, .count = length(value))

bed_map(x, y, .vals = values(value))

# count defaults are NA not 0; differs from bedtools2 ...
bed_map(x, y, .counts = dplyr::n())

# ... but NA counts can be coverted to 0's
dplyr::mutate(
  bed_map(x, y, .counts = dplyr::n()),
  .counts = ifelse(is.na(.counts), 0, .counts)
)

# a bigWig or bigBed file path (or `http(s)://` URL) can be used in place of
# `y`; only the regions spanned by `x` are read from the file
peaks <- tibble::tribble(
  ~chrom  , ~start   , ~end     ,
  "chr22" , 16100000 , 16150000 ,
  "chr22" , 16500000 , 16550000
)

bed_map(peaks, valr_example("hg19.dnase1.bw"), .signal = sum(value))
