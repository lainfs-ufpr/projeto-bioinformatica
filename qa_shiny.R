library(plotly)
library(viridis)
library(bslib)
library(shiny)
library(shinyjs)
library(shinyFiles)
library(ShortRead)
library(ggplot2)
library(fs)
library(scales)  


# Carregar pacotes

source("funcoes_prontas/plotAdapterContamination.R")
source("funcoes_prontas/tableAdapterContamination.R")
source("funcoes_prontas/freqSequences.R")
source("funcoes_prontas/plotCycleQuality.R")
source("funcoes_prontas/plotNucleotideCount.R")
source("funcoes_prontas/readQualityScore.R")
source("funcoes_prontas/plotOcurrences.R")

# Tamanho máximo de input
options(shiny.maxRequestSize = 100 * 1024^3)

# Front end
ui <- fluidPage(
  
  theme = bs_theme(version = 5, bg = "#FFFFFF", fg = "#31231a", 
                   bootswatch = "flatly", primary = "#1a754f", 
                   base_font = font_google("Lexend Deca")),
  useShinyjs(),
  
  # Procura o arquivo css na pasta www do diretório onde está rodando
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  
  #tags$h1("Controle de Qualidade de Sequências FASTQ", class = "titulo-app"),
  
  tags$div(
    tags$img(src = "logo-principal.png", id = "logo-fixo"),
    style = "position: fixed; bottom: 20px; left: 10px; z-index: 999;"
  ),
  
  
  
  page_navbar(
    input_dark_mode(id = "mode"), 
    textOutput("mode"),

    title = 
      tags$span("Controle de Qualidade de Sequências FASTQ", class = "titulo-app"),
    id = "page",
 
  nav_panel(
    "Análise",
    
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-custom",
      
      tags$h4("Para inicializar a análise, selecione o(s) arquivo(s) ou uma pasta", 
              class = "titulo-sidebar"),
      
      div(
        style = "display:flex; gap: 10px; margin-bottom: 20px; align-items: flex-start;", #ta com flex-start mas adicionar uma pasta nao esta alinhado com botao de input, esta alinhando com o texto
        
     #   div(class = "botao-diretorio",
     #   shinyDirButton("diretorio", "Adicionar uma pasta", "Selecionar")
      #  ),
        
        
        shinyDirButton("diretorio", "Adicionar uma pasta", "Selecionar"),
        
        #mais uma div para empilhar os botoes de input com o de rodar qualidade
        div(style = "display:flex; flex-direction:column; gap:10px;",
            
        tags$span("Adicionar um arquivo",class = "message-input"),
        div(class = "file-input-custom", 
            fileInput("arquivos",label = NULL , multiple = TRUE,
                      accept = c(".fasta", ".fa", ".fastq", ".fq", "text/plain"))),
        
      
      
      div(id = "loading_animation", class = "loading-spinner",
         style = "display: none;"),
      
      # Botão iniciar QA
      
        actionButton("run_analysis", "Rodar controle de qualidade",
                     class = "btnQA-custom")
        )
      ),
      # Aumentar espaço
      tags$br(),
      tags$br(),
      
      selectInput("palette_choice", "Mude a paleta de cores aqui:",
                  choices = c("viridis", "magma", "plasma", "rocket",
                              "cividis", "inferno", "turbo", "mako"),
                  selected = "viridis")
    ),
    
    
    # Paineis de gráficos
    mainPanel(
      navset_pill(
        nav_panel("Qualidade por Ciclo",
                 class = "titulo-plots",
                 plotOutput("plot_qualidade_ciclo_est"),
                 shinyjs::hidden(
                 downloadButton("download_plot_ciclo", "Baixar gráfico"))),
        nav_panel("Qualidade Média",
                 class = "titulo-plots",
                 plotOutput("plot_qualidade_media_est"),
                 shinyjs::hidden(
                 downloadButton("download_plot_media", "Baixar gráfico"))),
        nav_panel("Contagem de Bases",
                 class = "titulo-plots",
                 plotOutput("plot_contagens_est"),
                 shinyjs::hidden(
                 downloadButton("download_plot_contagens", "Baixar gráfico"))),
        nav_panel("Distribuição Cumulativa de Leituras",
                 class = "titulo-plots",
                 plotOutput("plot_ocorrencias_est"),
                 shinyjs::hidden(
                 downloadButton("download_plot_ocorrencias", "Baixar gráfico"))),
        nav_panel("Sequências Frequentes",
                 class = "titulo-plots",
                 tableOutput("tabela_frequencias")),
        nav_panel("Contaminação por Adaptadores",
                 class = "titulo-plots",
                 plotOutput("plot_adapters_est"),
                 shinyjs::hidden(
                 downloadButton("download_plot_adapters", "Baixar gráfico")),
                 tableOutput("tabela_adapters")),
      ),
      br(), hr(),
      div(id = "qa_output"),
  
  div(style = "margin-left: 20px;",
      textOutput("caminhoPasta"))
      )
    )
  )
  )
)
  


# Backend - servidor
server <- function(input, output, session) {
  
  # Usa o conteudo html de um arquivo separado
  html_content <- HTML(paste(readLines("www/apresentacao.html",
                                       encoding = "UTF-8"), collapse = "\n"))
  
  # Cria uma tela inicial como se fosse uma caixa de dialogo
  showModal(modalDialog(
    title = NULL,
    html_content,
    size = "l",
    easyClose = TRUE,
    footer = modalButton("Iniciar")
  ))
  
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)
  resultado_qa <- reactiveVal(NULL)
  
  observeEvent(input$diretorio, {
    path <- parseDirPath(volumes, input$diretorio)
    dir_path(path)
  })
  
  output$caminhoPasta <- renderText({
    req(dir_path())
    paste("Pasta selecionada:", dir_path())
  })
  
  # Rodar analise
  observeEvent(input$run_analysis, {
    shinyjs::hide("run_analysis")
    shinyjs::show("loading_animation")
    shinyjs::html("qa_output",
                  '<b style="color: #31231a">Controle de Qualidade em andamento...</b>')
    
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
      shinyjs::show("download_plot_adapters")
      shinyjs::show("download_plot_ocorrencias")
      
      
      tempo_execucao <- Sys.time() - tempo_inicio
      shinyjs::html("qa_output",
                    paste0('<b style="color: #31231a">Controle de Qualidade finalizado! Tempo de execução:</b> ',
                           round(tempo_execucao, 2), ' segundos'))
      shinyjs::hide("loading_animation")
      shinyjs::show("run_analysis")
      
    }, error = function(e) {
      shinyjs::html("qa_output",
                    paste0('<b style="color: #940e01">Erro: ', e$message, '</b>'))
      shinyjs::hide("loading_animation")
      shinyjs::show("run_analysis")
    })
  })
  
  # Paleta de cores escolhida
  paleta_cores <- reactive({
    req(input$palette_choice)
    tolower(input$palette_choice)
  })
  
  # Nomes originais dos arquivos para ajustar
  nomes_arquivos <- reactive({
    req(input$arquivos)
    input$arquivos$name
  })
  
  # Plot qualidade ciclo
  output$plot_qualidade_ciclo_est <- renderPlot({
    req(resultado_qa())
    plotCycleQuality(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_qualidade_ciclo_int <- renderPlot({
    req(resultado_qa())
    plotCycleQuality(resultado_qa(), paleta_cores(), nomes_arquivos())$p_interativo
  })
  
  output$download_plot_ciclo <- downloadHandler(
    filename = function() { "qualidade_ciclo.png" },
    content = function(file) {
      req(resultado_qa(), paleta_cores())
      p <- plotCycleQuality(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
      ggsave(file, p,
             width = 8, height = 6)
    }
  )
  
  # Plot qualidade média 
  output$plot_qualidade_media_est <- renderPlot({
    req(resultado_qa())
    readQualityScore(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_qualidade_media_int <- renderPlot({
    req(resultado_qa())
    readQualityScore(resultado_qa(), paleta_cores(), nomes_arquivos())$p_interativo
  })
  
  output$download_plot_media <- downloadHandler(
    filename = function() { "qualidade_media.png" },
    content = function(file) {
      req(resultado_qa())
      p <- readQualityScore(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
      ggsave(file, p,
             width = 8, height = 6)
    }
  )
  
  # Plot contagem bases
  output$plot_contagens_est <- renderPlot({
    fls <- if (!is.null(input$arquivos)) {
      input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      NULL
    }
    req(fls)
    plotNucleotideCount(fls, paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_contagens_int <- renderPlot({
    fls <- if (!is.null(input$arquivos)) {
      input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      NULL
    }
    req(fls)
    plotNucleotideCount(fls, paleta_cores(), nomes_arquivos())$p_interativo
  })

  output$download_plot_contagens <- downloadHandler(
    filename = function() { "contagem_bases.png" },
    content = function(file) {
      fls <- if (!is.null(input$arquivos)) {
        input$arquivos$datapath
      } else if (!is.null(dir_path())) {
        dir(dir_path(), pattern = "*fq$", full.names = TRUE)
      } else {
        NULL
      }
      req(fls)
      p <- plotNucleotideCount(fls, paleta_cores(), nomes_arquivos())$p_estatico
      ggsave(file, p,
             width = 8, height = 6)
    }
  )
  
  # Plot adapters
  output$plot_adapters_est <- renderPlot({
    fls <- if (!is.null(input$arquivos)) {
      input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      NULL
    }
    req(fls)
    plotAdapterContamination(fls, paleta_cores())$p_estatico
  })
  
  output$plot_adapters_int <- renderPlot({
    fls <- if (!is.null(input$arquivos)) {
      input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      NULL
    }
    req(fls)
    plotAdapterContamination(fls, paleta_cores())$p_interativo
  })
  
  output$download_plot_adapters <- downloadHandler(
    filename = function() { "contaminacao_adaptadores.png" },
    content = function(file) {
      fls <- if (!is.null(input$arquivos)) {
        input$arquivos$datapath
      } else if (!is.null(dir_path())) {
        dir(dir_path(), pattern = "*fq$", full.names = TRUE)
      } else {
        NULL
      }
      req(fls)
      p <- plotAdapterContamination(fls, paleta_cores())$p_estatico
      ggsave(file, p, width = 8, height = 6)
    }
  )
  
  # Plot ocorrencias
  output$plot_ocorrencias_est <- renderPlot({
    req(resultado_qa())
    plotOcurrences(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
  })
  
  output$plot_ocorrencias_int <- renderPlot({
    req(resultado_qa())
    plotOcurrences(resultado_qa(), paleta_cores(), nomes_arquivos())$p_interativo
  })
  
  output$download_plot_ocorrencias <- downloadHandler(
    filename = function() { "distribuicao_ocorrencias.png" },
    content = function(file) {
      req(resultado_qa())
      plotOcurrences(resultado_qa(), paleta_cores(), nomes_arquivos())$p_estatico
      ggsave(file, p,
             width = 8, height = 6)
    }
  )

  # Tabela sequencias frequentes
  output$tabela_frequencias <- renderTable({
    req(resultado_qa())
    t <- freqSequences(resultado_qa())
    t$Arquivo <- input$arquivos$name[
      match(t$Arquivo,
            unique(t$Arquivo))
    ]
    t$Arquivo <- ifelse(duplicated(t$Arquivo),
                                    "",
                                    t$Arquivo)
    t
  })
  
  # Tabela adapters
  output$tabela_adapters <- renderTable({
    fls <- if (!is.null(input$arquivos)) {
      input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    } else {
      NULL
    }
    req(fls)
    t <- tableAdapterContamination(fls)
    t$Arquivo <- input$arquivos$name[
      match(t$Arquivo,
            unique(t$Arquivo))
    ]
    t
  })
  
}

shinyApp(ui = ui, server = server)
