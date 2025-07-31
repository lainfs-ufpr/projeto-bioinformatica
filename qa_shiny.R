source("funcoes_prontas/plotAdapterContamination.R")
source("funcoes_prontas/tableAdapterContamination.R")
source("funcoes_prontas/freqSequences.R")
source("funcoes_prontas/plotCycleQuality.R")
source("funcoes_prontas/plotNucleotideCount.R")
source("funcoes_prontas/readQualityScore.R")
source("funcoes_prontas/plotOcurrences.R")

# Instale o BiocManager se ainda não estiver instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instale os pacotes
BiocManager::install(c("shiny", "shinyFiles", "ShortRead"), force=TRUE)

library(stats)
library(shiny)
library(shinyFiles)
library(ShortRead)
library(shinyjs)

# Aumenta o limite máximo de upload para 100MB
options(shiny.maxRequestSize = 100 * 1024^3) 

# Interface do usuário
ui <- fluidPage(
  shinyjs::useShinyjs(),  
  tags$style(HTML("
   #qa_output {
          position: fixed;
          bottom: 0;
          left: 0;
          width: 100%;
          background-color: #f2f2f2;
          padding: 10px;
          text-align: center;
          border-top: 1px solid #ccc;
          z-index: 9999;  
    }
  ")),
  tags$style(HTML("
    table {
      font-size: 10px;
    }
  ")),

  titlePanel("Seleção de Pasta, Upload e Controle de Qualidade"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("diretorio", "Escolher Pasta", "Selecionar"),
      verbatimTextOutput("caminhoPasta"),
      fileInput("arquivos", "Escolha os arquivos", 
                multiple = TRUE, 
                accept = c(".fasta", ".fa", ".fastq", ".fq", "text/plain")),
      actionButton("run_analysis", "Rodar Controle de Qualidade")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Qualidade por Ciclo",
                 plotOutput("plot_qualidade_ciclo")),
        tabPanel("Qualidade Média",
                 plotOutput("plot_qualidade_media")),
        tabPanel("Contagem de Bases",
                 plotOutput("plot_contagens")),
        tabPanel("Sequências Frequentes",
                 tableOutput("tabela_frequencias")),
        tabPanel("Contaminação por Adaptadores",
                 div(
                   plotOutput("plot_adapters"),
                   br(),  
                   tableOutput("tabela_adapters")
                 )
                )
      )
    )
  ),
  div(id = "qa_output")
)

# Lógica do servidor
server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)

  dir_path <- reactiveVal(NULL)
  plot_qualidade_ciclo_reactive <- reactiveVal(NULL)
  plot_qualidade_media_reactive <- reactiveVal(NULL)
  plot_contagens_reactive <- reactiveVal(NULL)
  tabela_frequencias_reactive <- reactiveVal(NULL)
  plot_adapters_reactive <- reactiveVal(NULL)
  tabela_adapters_reactive <- reactiveVal(NULL)

  observeEvent(input$diretorio, {
    path <- parseDirPath(volumes, input$diretorio)
    dir_path(path)

  })
  
  output$caminhoPasta <- renderText({
    req(dir_path())
    paste("Pasta selecionada:", dir_path())
  })
  
  output$conteudoArquivos <- renderTable({
    req(input$arquivos)
    data.frame(Nome_Arquivo = input$arquivos$name, Caminho = input$arquivos$datapath)
  })
  
  observeEvent(input$run_analysis, {
    # Atualiza imediatamente com "QC em andamento..."
    shinyjs::html("qa_output", '<b style="color:black">QC em andamento...</b>')
    
    # Determina arquivos
    if (!is.null(input$arquivos)) {
      fls <- input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      fls <- dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      showNotification("Nenhum arquivo fornecido!", type = "error")
      return()
    }
    
    if (length(fls) == 0) {
      showNotification("Nenhum arquivo FASTQ encontrado!", type = "error")
      return()
    }
    
    # print(fls)
    
    tryCatch({
    
      # QC
      tempo_inicio <- Sys.time()
      
      # Função original do ShortRead
      resultado_qa <- qa(fls, type = "fastq")

      # Plots personalizados
      print("Plot ciclo")
      qualidade_ciclo <- plotCycleQuality(resultado_qa)
      plot_qualidade_ciclo <- qualidade_ciclo$p_estatico
      plot_qualidade_ciclo$data$Arquivo <- input$arquivos$name[
        match(plot_qualidade_ciclo$data$Arquivo,
              unique(plot_qualidade_ciclo$data$Arquivo))
      ]
      plot_qualidade_ciclo_reactive(plot_qualidade_ciclo)
      
      print("Qualidade média")
      plot_qualidade_media <- readQualityScore(resultado_qa)
      plot_qualidade_media$data$Arquivo <- input$arquivos$name[
        match(plot_qualidade_media$data$Arquivo,
              unique(plot_qualidade_media$data$Arquivo))
      ]
      plot_qualidade_media_reactive(plot_qualidade_media)
      
      print("Contagens")
      plot_contagens <- plotNucleotideCount(fls)
      plot_contagens$data$Arquivo <- input$arquivos$name[
        match(plot_contagens$data$Arquivo,
              unique(plot_contagens$data$Arquivo))
      ]
      plot_contagens_reactive(plot_contagens)
      
      print("Freq seq")
      tabela_frequencias <- freqSequences(resultado_qa)
      tabela_frequencias$Arquivo <- input$arquivos$name[
        match(tabela_frequencias$Arquivo,
              unique(tabela_frequencias$Arquivo))
      ]
      tabela_frequencias$Arquivo <- ifelse(duplicated(tabela_frequencias$Arquivo),
                                           "",
                                           tabela_frequencias$Arquivo)
      tabela_frequencias_reactive(tabela_frequencias)
      
      print("Adapters plot")
      plot_adapters <- plotAdapterContamination(fls)
      plot_adapters_reactive(plot_adapters)
      
      print("Adapters table")
      tabela_adapters <- tableAdapterContamination(fls)
      tabela_adapters_reactive(tabela_adapters)
      
      tempo_fim <- Sys.time()
      tempo_execucao <- tempo_fim - tempo_inicio
      
      # Atualiza resultado final
      shinyjs::html("qa_output",
                    paste0('<b style="color:black">QC finalizado! Tempo de execução:</b> ',
                           round(tempo_execucao, 2), ' segundos'))
    }, error = function(e) {
      print(paste("Erro detectado:", e$message))
    })
    
  })
  
  output$plot_qualidade_ciclo <- renderPlot({
    req(plot_qualidade_ciclo_reactive())
    plot_qualidade_ciclo_reactive()
  })
  
  output$plot_contagens <- renderPlot({
    req(plot_contagens_reactive())
    plot_contagens_reactive()
  })
  
  output$plot_qualidade_media <- renderPlot({
    req(plot_qualidade_media_reactive())
    plot_qualidade_media_reactive()
  })
  
  output$tabela_frequencias <- renderTable({
    req(tabela_frequencias_reactive())
    tabela_frequencias_reactive()
  })
  
  output$plot_adapters <- renderPlot({
    req(plot_adapters_reactive())
    plot_adapters_reactive()
  })
  
  output$tabela_adapters <- renderTable({
    req(tabela_adapters_reactive())
    tabela_adapters_reactive()
  })
}

# Inicializa o app Shiny
shinyApp(ui = ui, server = server)
