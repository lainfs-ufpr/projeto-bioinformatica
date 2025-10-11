library(shiny)

ui <- fluidPage(
  titlePanel("Trimagem com Trimmomatic"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fastq", "Selecione arquivo FASTQ", accept = ".fastq"),
      textInput("adapters", "Arquivo de adaptadores",
                value = "Trimmomatic-0.39/adapters/TruSeq3-SE.fa"),
      numericInput("sliding_window", "Tamanho da janela (SLIDINGWINDOW)", 4),
      numericInput("quality_cutoff", "Qualidade mínima", 20),
      numericInput("minlen", "Comprimento mínimo (MINLEN)", 50),
      textInput("output", "Nome do arquivo de saída", "reads_trimmed.fastq"),
      actionButton("run", "Rodar Trimmomatic")
    ),
    
    mainPanel(
      verbatimTextOutput("log"),
      tableOutput("stats")
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$run, {
    req(input$fastq)
    
    trimmomatic <- "Trimmomatic-0.39/trimmomatic-0.39.jar"
    input_fastq <- input$fastq$datapath
    output_fastq <- input$output
    
    cmd <- paste(
      "java -jar", trimmomatic, "SE -phred33",
      input_fastq, output_fastq,
      paste0("ILLUMINACLIP:", input$adapters, ":2:30:10"),
      paste0("SLIDINGWINDOW:", input$sliding_window, ":", input$quality_cutoff),
      paste0("MINLEN:", input$minlen)
    )
    
    # Rodar Trimmomatic
    system(cmd, intern = TRUE)
    
    output$log <- renderText({
      paste("Trimagem concluída!\nArquivo salvo em:", output_fastq)
    })
    
    # Estatísticas antes/depois (usando ShortRead)
    library(ShortRead)
    raw <- readFastq(input_fastq)
    trimmed <- readFastq(output_fastq)
    
    stats <- data.frame(
      Arquivo = c("Original", "Trimado"),
      Reads = c(length(raw), length(trimmed)),
      Tam_médio = c(mean(width(sread(raw))),
                    mean(width(sread(trimmed))))
    )
    
    output$stats <- renderTable(stats)
  })
}

shinyApp(ui, server)
