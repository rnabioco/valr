library(shinydashboard)

# options for saving
savingOptions <- list(
  "dom" = 'lBfrtip',
  buttons = list('copy', 'print', list(
    extend = 'collection',
    buttons = list(list(extend = 'csv',
                        fieldBoundary = ""
    ),
    list(extend = 'csv',
         fieldSeparator = "\t",
         fieldBoundary = "",
         text = "TSV",
         extension = ".tsv"),
    list(extend = 'pdf')
    ),
    text = 'Download'
  )),
  list(orderClasses = TRUE),
  lengthMenu = list(c(10, 25, -1), c('10', '25', 'All'))
)


bedfile <- system.file('extdata', 'genes.hg19.chr22.bed.gz', package = 'valr')
bgfile  <- system.file('extdata', 'hela.h3k4.chip.bg.gz', package = 'valr')
genomefile <- system.file('extdata', 'hg19.chrom.sizes.gz', package = 'valr')


genes <- read_bed(bedfile, n_fields = 6)

plot_metagene <- function(window_size, region){
  bedfile <- system.file('extdata', 'genes.hg19.chr22.bed.gz', package = 'valr')
  bgfile  <- system.file('extdata', 'hela.h3k4.chip.bg.gz', package = 'valr')
  genomefile <- system.file('extdata', 'hg19.chrom.sizes.gz', package = 'valr')
  
  genes <- read_bed(bedfile, n_fields = 6)
  genome <- read_genome(genomefile)
  
  # generate 1 bp intervals, + strand only for now
  tss <- genes %>%
    filter(strand == '+') %>%
    mutate(end = start + 1)
  
  region_size <- region
  win_size <- window_size
  
  x <- tss %>%
    bed_slop(genome, both = region_size) %>%
    bed_makewindows(genome, win_size) %>%
    group_by(win_id)
  
  y <- read_bedgraph(bgfile)
  
  res <- bed_map(x, y, sums = sum(value.y)) %>%
    summarize(means = mean(sums), sds = sd(sums))
  
  x_labels <- pretty(-region_size:region_size, n = 11)
  x_breaks <- seq(1, length(res$win_id.x), length.out = length(x_labels))
  sd_limits <- aes(ymax = means + sds, ymin = means - sds)
  
  ggplot(res, aes(x = win_id.x, y = means)) +
    geom_point() + geom_pointrange(sd_limits) + 
    scale_x_continuous(labels = x_labels, breaks = x_breaks) + 
    ggtitle('H3K4me3 ChIP signal near TSSs') +
    xlab('Position\n(bp from TSS)') + ylab('Signal') +
    theme_bw()
}



ui <- dashboardPage(
  dashboardHeader(title = "valr", titleWidth = 100),
  dashboardSidebar(
    width = 100, 
    menuItem("Table", tabName = "table", icon = icon("table")),
    menuItem("Plot", tabName =  "plot", icon = icon("signal", lib = "glyphicon"), selected = T)
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "table",
        fluidRow(
          box(DT::dataTableOutput("mytable1"))
        )),
      tabItem(tabName = "plot",
        fluidRow(
          column(width = 3,
              box(
                title = "Controls",
                sliderInput("window", "Window Size:", 1, 100, 50),
                sliderInput("region", "Region Size:", 100, 10000, 1000),
                width = 12
            )
          ),     
          column(width = 9,  
            box(
              title = "Coverage Plot",
              plotOutput("coverage"),
              width = 12
            )
          )
          )
        )     
        
      )
    )
  )






server <- function(input, output) {
  
  # a large table, reative to input$show_vars
  output$mytable1 = DT::renderDataTable({
    genes 
  }, extensions = c('Buttons'),
  options = savingOptions, 
  filter = 'top', 
  style = 'bootstrap',
  selection = list(mode = 'single', 
                   target = 'row'),
  rownames = FALSE)
  
  output$window <- renderUI({
    #validate(
    # need(input$mytable1_rows_selected != "",
    #       "")
    #)
    sliderInput("window",
                "Window size",
                min = 1,
                max = 100,
                value = 50)
  })
  
  output$region <- renderUI({
    #validate(
    # need(input$mytable1_rows_selected != "",
    #       "")
    #)
    sliderInput("region",
                "Region size",
                min = 100,
                max = 10000,
                value = 1000)
  })
  
  output$coverage = renderPlot({
    #validate(
    #  need(input$mytable1_rows_selected != "",
    #       "Select an entry to display coverage data")
    #)
    #plot_coverage(input$mytable1_rows_selected, 
    #              in_bed, 
    #              in_bedgraph, 
    #              in_bed12,
    #              zoom_factor = input$Zoom)
    plot_metagene(window = input$window, 
                  region = input$region)
  })
}



shiny::shinyApp(ui, server)



