library(ShortRead)

source("scripts/funcoes_prontas/plotAdapterContamination.R")
source("scripts/funcoes_prontas/tableAdapterContamination.R")
source("scripts/funcoes_prontas/freqSequences.R")
source("scripts/funcoes_prontas/plotCycleQuality.R")
source("scripts/funcoes_prontas/plotNucleotideCount.R")
source("scripts/funcoes_prontas/readQualityScore.R")
source("scripts/funcoes_prontas/plotOcurrences.R")

# BACK
server <- function(input, output, session) {
  
  # Usa o conteudo html de um arquivo separado
  html_content <- HTML(paste(readLines("www/apresentacao.html",
                                       encoding = "UTF-8"), collapse = "\n"))
  
  # Cria uma tela inicial como se fosse uma caixa de dialogo
  showModal(modalDialog(title = NULL,
                        html_content,
                        size = "l",
                        easyClose = TRUE,
                        footer = modalButton("Iniciar"))
  )
  
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)
  resultado_qa <- reactiveVal(NULL)
  
  observeEvent(input$diretorio,
               {path <- parseDirPath(volumes, input$diretorio)
               dir_path(path)}
  )
  
  output$caminhoPasta <- renderText({req(dir_path())
    paste("Pasta selecionada:", dir_path())})
  
  # Rodar analise
  observeEvent(input$run_analysis,
               { shinyjs::hide("run_analysis")
                 shinyjs::show("loading_animation")
                 shinyjs::html("qa_output",'<b style="color: black;">Controle de Qualidade em andamento...</b>')
                 
                 fls <- if (!is.null(input$arquivos)) {
                   input$arquivos$datapath
                 } else if (!is.null(dir_path())) {
                   dir(dir_path(), pattern = "*fq$", full.names = TRUE)
                 } else {
                   showNotification("Nenhum arquivo fornecido!", type = "error")
                   shinyjs::show("run_analysis")
                   shinyjs::hide("loading_animation")
                   return()
                 }
                 
                 if (length(fls) == 0) {
                   showNotification("Nenhum arquivo FASTQ encontrado!", type = "error")
                   shinyjs::show("run_analysis")
                   shinyjs::hide("loading_animation")
                   return()
                 }
                 
                 tryCatch({
                   tempo_inicio <- Sys.time()
                   resultado <- qa(fls, type = "fastq")
                   resultado_qa(resultado)
                   
                   shinyjs::show("download_plot_ciclo")
                   shinyjs::show("download_plot_media")
                   shinyjs::show("download_plot_contagens")
                   shinyjs::show("download_plot_ocorrencias")
                   shinyjs::show("download_plot_adapters")
                   
                   tempo_execucao <- Sys.time() - tempo_inicio
                   shinyjs::html("qa_output",
                                 paste0('<b style="color: #31231a">Controle de Qualidade finalizado! Tempo de execução:</b> ',
                                        round(tempo_execucao, 2), 
                                        ' segundos'))
                   shinyjs::hide("loading_animation")
                   shinyjs::show("run_analysis")
                   
                 }, 
                 error = function(e) {
                   shinyjs::html("qa_output",
                                 paste0('<b style="color: #940e01">Erro: ', e$message, '</b>'))
                   shinyjs::hide("loading_animation")
                   shinyjs::show("run_analysis")
                 })
               })
  
  observeEvent(input$run_trim, {
    req(input$fastq_trim)
    
    # Caminhos
    trimmomatic <- "Trimmomatic-0.39/trimmomatic-0.39.jar"
    input_fastq <- input$fastq_trim$datapath
    output_fastq <- file.path(dirname(input_fastq), input$output_trim)
    
    # Comando do Trimmomatic
    cmd <- paste(
      "java -jar", shQuote(trimmomatic), "SE -phred33",
      shQuote(input_fastq), shQuote(output_fastq),
      paste0("ILLUMINACLIP:", input$adapters, ":2:30:10"),
      paste0("SLIDINGWINDOW:", input$sliding_window, ":", input$quality_cutoff),
      paste0("MINLEN:", input$minlen)
    )
    
    shinyjs::html("trim_log", "<b style='color:black'>Rodando Trimmomatic...</b>")
    
    # Executar Trimmomatic
    resultado <- tryCatch({
      system(cmd, intern = TRUE)
    }, error = function(e) {
      paste("Erro:", e$message)
    })
    
    output$trim_log <- renderText({
      paste("Trimagem concluída!\nArquivo salvo em:", output_fastq)
    })
    
    # Calcular estatísticas antes e depois
    library(ShortRead)
    raw <- readFastq(input_fastq)
    trimmed <- readFastq(output_fastq)
    
    stats <- data.frame(
      Arquivo = c("Original", "Trimado"),
      Reads = c(length(raw), length(trimmed)),
      Tam_médio = c(mean(width(sread(raw))),
                    mean(width(sread(trimmed))))
    )
    
    output$trim_stats <- renderTable(stats)
  })
  
  # Paleta de cores escolhida
  paleta_cores <- reactive({req(input$palette_choice)
    tolower(input$palette_choice)
  })
  
  # Nomes originais dos arquivos para ajustar
  nomes_arquivos <- reactive({req(input$arquivos)
    input$arquivos$name
  })
  
  # Plot qualidade ciclo
  output$plot_qualidade_ciclo_est <- renderPlot({req(resultado_qa())
    plotCycleQuality(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_qualidade_ciclo_int <- renderPlot({req(resultado_qa())
    plotCycleQuality(resultado_qa(), paleta_cores(), nomes_arquivos())$p_interativo
  })
  
  output$download_plot_ciclo <- downloadHandler(
    filename = function() {
      paste0("qualidade_ciclo.", input$formato_download_ciclo)
    },
    content = function(file) {
      req(resultado_qa(), paleta_cores(), input$formato_download_ciclo)
      
      p <- plotCycleQuality(
        resultado_qa(),
        paleta_cores(),
        nomes_arquivos()
      )$p_estatico
      
      ggsave(
        filename = file,
        plot = p,
        width = 8, height = 6,
        device = input$formato_download_ciclo
      )
    }
  )
  
  # Plot qualidade média 
  output$plot_qualidade_media_est <- renderPlot({req(resultado_qa())
    readQualityScore(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_qualidade_media_int <- renderPlot({req(resultado_qa())
    readQualityScore(resultado_qa(), paleta_cores(), nomes_arquivos())$p_interativo
  })
  
  output$download_plot_media <- downloadHandler(
    filename = function() {
      paste0("media.", input$formato_download_media)
    },
    content = function(file) {
      req(resultado_qa(), paleta_cores(), input$formato_download_media)
      
      p <- readQualityScore(
        resultado_qa(),
        paleta_cores(),
        nomes_arquivos()
      )$p_estatico
      
      ggsave(
        filename = file,
        plot = p,
        width = 8, height = 6,
        device = input$formato_download_media
      )
    }
  )
  
  # Plot contagem bases
  output$plot_contagens_est <- renderPlot({fls <- if (!is.null(input$arquivos)) {
    input$arquivos$datapath
  } else if (!is.null(dir_path())) {
    dir(dir_path(), pattern = "*fq$", full.names = TRUE)
  } else {NULL}
  req(fls)
  plotNucleotideCount(fls, paleta_cores(), 
                      nomes_arquivos())$p_estatico
  })
  
  output$plot_contagens_int <- renderPlot({fls <- if (!is.null(input$arquivos)) {
    input$arquivos$datapath
  } else if (!is.null(dir_path())) {
    dir(dir_path(), pattern = "*fq$", full.names = TRUE)
  } else {NULL}
  req(fls)
  plotNucleotideCount(fls, paleta_cores(), 
                      nomes_arquivos())$p_interativo
  })
  
  output$download_plot_contagens <- downloadHandler(
    filename = function() {
      paste0("contagens.", input$formato_download_contagens)
    },
    content = function(file) {
      req(resultado_qa(), paleta_cores(), input$formato_download_contagens)
      
      p <- plotNucleotideCount(
        resultado_qa(),
        paleta_cores(),
        nomes_arquivos()
      )$p_estatico
      
      ggsave(
        filename = file,
        plot = p,
        width = 8, height = 6,
        device = input$formato_download_contagens
      )
    }
  )
  
  # Plot adapters
  output$plot_adapters_est <- renderPlot({fls <- if (!is.null(input$arquivos)) {
    input$arquivos$datapath
  } else if (!is.null(dir_path())) {
    dir(dir_path(), pattern = "*fq$", full.names = TRUE)
  } else {NULL}
  req(fls)
  plotAdapterContamination(fls, paleta_cores())$p_estatico
  })
  
  output$plot_adapters_int <- renderPlot({fls <- if (!is.null(input$arquivos)) {
    input$arquivos$datapath
  } else if (!is.null(dir_path())) {
    dir(dir_path(), pattern = "*fq$", full.names = TRUE)
  } else {NULL}
  req(fls)
  plotAdapterContamination(fls, paleta_cores())$p_interativo
  })
  
  output$download_plot_adapters <- downloadHandler(
    filename = function() {
      paste0("adapters_Plot.", input$formato_download_adaptersP)
    },
    content = function(file) {
      req(resultado_qa(), paleta_cores(), input$formato_download_adaptersP)
      
      p <- plotAdapterContamination(
        resultado_qa(),
        paleta_cores(),
        nomes_arquivos()
      )$p_estatico
      
      ggsave(
        filename = file,
        plot = p,
        width = 8, height = 6,
        device = input$formato_download_adaptersP
      )
    }
  )
  
  # Plot ocorrencias
  output$plot_ocorrencias_est <- renderPlot({req(resultado_qa())
    plotOcurrences(resultado_qa(), 
                   paleta_cores(), 
                   nomes_arquivos())$p_estatico
  })
  
  output$plot_ocorrencias_int <- renderPlot({req(resultado_qa())
    plotOcurrences(resultado_qa(), 
                   paleta_cores(), 
                   nomes_arquivos())$p_interativo
  })
  
  output$download_plot_ocorrencias <- downloadHandler(
    filename = function() {
      paste0("ocorrencias.", input$formato_download_ocorrencias)
    },
    content = function(file) {
      req(resultado_qa(), paleta_cores(), input$formato_download_ocorrencias)
      
      p <- plotOcurrences(
        resultado_qa(),
        paleta_cores(),
        nomes_arquivos()
      )$p_estatico
      
      ggsave(
        filename = file,
        plot = p,
        width = 8, height = 6,
        device = input$formato_download_ocorrencias
      )
    }
  )
  
  # Tabela sequencias frequentes
  output$tabela_frequencias <- renderTable({req(resultado_qa())
    t <- freqSequences(resultado_qa())
    t$Arquivo <- input$arquivos$name[match(t$Arquivo,
                                           unique(t$Arquivo))]
    t$Arquivo <- ifelse(duplicated(t$Arquivo),"",t$Arquivo)
    t
  })
  
  # Tabela adapters
  output$tabela_adapters <- renderTable({fls <- if (!is.null(input$arquivos)) {
    input$arquivos$datapath
  }else if (!is.null(dir_path())) {
    dir(dir_path(), pattern = "*fq$", full.names = TRUE)
  }else {NULL}
  
  req(fls)
  t <- tableAdapterContamination(fls)
  t$Arquivo <- input$arquivos$name[
    match(
      t$Arquivo,unique(t$Arquivo))]
  t
  })
  
}