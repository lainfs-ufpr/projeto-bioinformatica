# Função Server do Módulo de Trimagem
trimagemServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    observeEvent(input$run_trim, {
      req(input$fastq_trim)
      
      # Paths
      trimmomatic <- "Trimmomatic-0.39/trimmomatic-0.39.jar"
      input_fastq <- input$fastq_trim$datapath
      output_fastq <- file.path(dirname(input_fastq), input$output_trim)
      
      # Trimmomatic Command
      cmd <- paste(
        "java -jar", shQuote(trimmomatic), "SE -phred33",
        shQuote(input_fastq), shQuote(output_fastq),
        paste0("ILLUMINACLIP:", input$adapters, ":2:30:10"),
        paste0("SLIDINGWINDOW:", input$sliding_window, ":", input$quality_cutoff),
        paste0("MINLEN:", input$minlen)
      )
      
      shinyjs::html(session$ns("trim_log"), "<b style='color:black'>Rodando Trimmomatic...</b>")
      
      # Executar Trimmomatic e capturar o status de saída
      status_saida <- tryCatch({
        # NOTE: Using ShortRead library functions requires it to be loaded, 
        # either globally (in server.R) or explicitly here.
        resultado <- system(cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)
        # Adicionar o log do trimmomatic
        shinyjs::html(session$ns("trim_log"), paste0("Log do Trimmomatic:\n<pre>", 
                                                     paste(resultado, collapse = "\n"), 
                                                     "</pre>"))
        attr(resultado, "status") # Retorna o status de saída
      }, error = function(e) {
        shinyjs::html(session$ns("trim_log"), paste0("Erro ao executar Trimmomatic: ", e$message))
        return(1) # Retorna um status de erro
      })
      
      # Check Trimmomatic success
      if (!is.null(status_saida) && status_saida != 0) {
        shinyjs::html(session$ns("trim_log"), paste0("<b style='color:red'>Erro na Trimagem. O Trimmomatic retornou um status de erro (", status_saida, "). Verifique o log acima e os caminhos.</b>"))
        return() 
      }
      
      # Check if output file exists
      if (!file.exists(output_fastq)) {
        shinyjs::html(session$ns("trim_log"), paste0("<b style='color:red'>Erro: Arquivo de saída esperado (", output_fastq, ") não foi encontrado após a execução do Trimmomatic.</b>"))
        return()
      }
      
      shinyjs::html(session$ns("trim_log"), paste("Trimagem concluída!\nArquivo salvo em:", output_fastq))
      
      # Calculate stats before and after
      raw <- ShortRead::readFastq(input_fastq)
      trimmed <- ShortRead::readFastq(output_fastq)
      
      stats <- data.frame(
        Arquivo = c("Original", "Trimado"),
        Reads = c(length(raw), length(trimmed)),
        Tam_médio = c(mean(width(ShortRead::sread(raw))),
                      mean(width(ShortRead::sread(trimmed))))
      )
      
      output$trim_stats <- renderTable(stats)
    })
  })
}