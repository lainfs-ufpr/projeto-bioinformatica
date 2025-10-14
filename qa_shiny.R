library(plotly)
library(viridis)
library(bslib)
library(shiny)
library(shinyjs)
library(shinyFiles)
library(shinyWidgets) #botão de download, precisa fazer o download dessa biblioteca e do svglite (para baixar svg)
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

#Tamanho máximo de input
options(shiny.maxRequestSize = 100 * 1024^3)

#Front
ui <- fluidPage(
  
  theme = bs_theme(version = 5, bg = "#FFFFFF", fg = "#31231a", 
                   bootswatch = "flatly", primary = "#1a754f", 
                   base_font = font_google("Lexend Deca")),
  useShinyjs(),
  
  #Link css
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  
  
  page_navbar(
      #Header
      title = tags$span(tags$img(src = "logo-principal.png", id = "logo-fixo",
                      style = "height:80px; margin-right:10px; vertical-align:middle;"), 
                      "Controle de Qualidade de Sequências FASTQ", 
                      class = "titulo-app"),
      #Painel de 'seleção' de função
      nav_panel("Qualidade"),
    
      nav_spacer(),
      
      #Controle de modo 
      nav_item(input_dark_mode(id = "mode", value = FALSE)), 
    
      id = "page",
    
      sidebarLayout(
          #sidebar vermelha
          sidebarPanel(width = 3,
                      class = "sidebar-custom",
                      tags$h4("Para inicializar a análise, selecione o(s) arquivo(s) ou uma pasta", 
                              class = "titulo-sidebar"),
      
                      div(style = "display:flex; flex-direction:column; gap:15px; margin-bottom:20px;",
                          div(style = "display:flex; gap:15px; align-items:flex-center;",
                              shinyDirButton(class = "botao-diretorio",
                                            "diretorio", "Adicionar uma pasta", "Selecionar"), 
                              div(style = "display:flex; flex-direction:column;",
                                  tags$span("Adicionar um arquivo", class = "message-input"),
                                  div(class = "botao-arquivo",
                                      fileInput("arquivos", label = NULL, multiple = TRUE,
                                                accept = c(".fasta", ".fa", ".fastq", ".fq", "text/plain"))
                                  ),
                              )
                          ),
                          div(id = "loading_animation", class = "loading-spinner", style = "display: none;"),
                          actionButton("run_analysis", "Rodar controle de qualidade", class = "btnQA-custom")
                      ),

                      tags$br(),
      
                      selectInput("palette_choice", "Mude a paleta de cores aqui:",
                                  choices = c("viridis", "magma", "plasma", "rocket",
                                              "cividis", "inferno", "turbo", "mako"),
                                  selected = "viridis")
          ),
    
    
          #Paineis de gráficos
          mainPanel(width = 9,
                    navset_pill(nav_panel("Qualidade por Ciclo",
                                          class = "titulo-plots",
                                          plotOutput("plot_qualidade_ciclo_est"),
                                          
                                          shinyjs::hidden(
                                            div(id = "download_plot_ciclo",
                                                dropdownButton(
                                                  inputId = "download_opts_ciclo",
                                                  label = "Baixar gráfico", 
                                                  icon = icon("download"),
                                                  status = "primary",
                                                  circle = FALSE,
                                                  inline = TRUE,
                                            
                                                selectInput("formato_download_ciclo", "Formato:",
                                                          choices = c("png", "pdf", "jpeg", "svg"),
                                                          selected = "pdf"),
                                                downloadButton("download_plot_ciclo", "Download")
                                              )
                                            )
                                          )
                                ),
                                nav_panel("Qualidade Média",
                                          class = "titulo-plots",
                                          plotOutput("plot_qualidade_media_est"),
                                          
                                          
                                          shinyjs::hidden(
                                            div(id = "download_plot_media",
                                                dropdownButton(
                                                  inputId = "download_opts_media",
                                                  label = "Baixar gráfico", 
                                                  icon = icon("download"),
                                                  status = "primary",
                                                  circle = FALSE,
                                                  inline = TRUE,
                                                  
                                                  selectInput("formato_download_media", "Formato:",
                                                              choices = c("png", "pdf", "jpeg", "svg"),
                                                              selected = "pdf"),
                                                  downloadButton("download_plot_media", "Download")
                                                )
                                            )
                                          )
                                ),
                                nav_panel("Contagem de Bases",
                                          class = "titulo-plots",
                                          plotOutput("plot_contagens_est"),
                                          
                                          shinyjs::hidden(
                                            div(id = "download_plot_contagens",
                                                dropdownButton(
                                                  inputId = "download_opts_contagens",
                                                  label = "Baixar gráfico", 
                                                  icon = icon("download"),
                                                  status = "primary",
                                                  circle = FALSE,
                                                  inline = TRUE,
                                                  
                                                  selectInput("formato_download_contagens", "Formato:",
                                                              choices = c("png", "pdf", "jpeg", "svg"),
                                                              selected = "pdf"),
                                                  downloadButton("download_plot_contagens", "Download")
                                                )
                                            )
                                          )
                                ),
                                nav_panel("Distribuição Cumulativa de Leituras",
                                          class = "titulo-plots",
                                          plotOutput("plot_ocorrencias_est"),
                                          
                                          shinyjs::hidden(
                                            div(id = "download_plot_ocorrencias",
                                                dropdownButton(
                                                  inputId = "download_opts_ocorrencias",
                                                  label = "Baixar gráfico", 
                                                  icon = icon("download"),
                                                  status = "primary",
                                                  circle = FALSE,
                                                  inline = TRUE,
                                                  
                                                  selectInput("formato_download_ocorrencias", "Formato:",
                                                              choices = c("png", "pdf", "jpeg", "svg"),
                                                              selected = "pdf"),
                                                  downloadButton("download_plot_ocorrencias", "Download")
                                                )
                                            )
                                          )
                                ),
                                nav_panel("Sequências Frequentes",
                                          class = "titulo-plots",
                                          tableOutput("tabela_frequencias")),
                                nav_panel("Contaminação por Adaptadores",
                                          class = "titulo-plots",
                                          plotOutput("plot_adapters_est"),
                                          
                                          shinyjs::hidden(
                                            div(id = "download_plot_adapters",
                                                dropdownButton(
                                                  inputId = "download_opts_adapters",
                                                  label = "Baixar gráfico", 
                                                  icon = icon("download"),
                                                  status = "primary",
                                                  circle = FALSE,
                                                  inline = TRUE,
                                                  
                                                  selectInput("formato_download_adaptersP", "Formato:",
                                                              choices = c("png", "pdf", "jpeg", "svg"),
                                                              selected = "pdf"),
                                                  downloadButton("download_plot_adapters", "Download")
                                                )
                                            )
                                          ),
                                          tableOutput("tabela_adapters")
                                ),
                    ),
                    br(),
                    hr(),
                    div(id = "qa_output"),
  
                    div(style = "margin-left: 20px;", textOutput("caminhoPasta"))
          )
      )
  )
)

#Backend - servidor
server <- function(input, output, session) {
  
  #Usa o conteudo html de um arquivo separado
  html_content <- HTML(paste(readLines("www/apresentacao.html",
                                       encoding = "UTF-8"), collapse = "\n"))
  
  #Cria uma tela inicial como se fosse uma caixa de dialogo
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
  
  #Rodar analise
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
  
  #Paleta de cores escolhida
  paleta_cores <- reactive({req(input$palette_choice)
                            tolower(input$palette_choice)
  })
  
  #Nomes originais dos arquivos para ajustar
  nomes_arquivos <- reactive({req(input$arquivos)
                              input$arquivos$name
  })
  
  #Plot qualidade ciclo
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
  
  #Plot qualidade média 
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
  
  #Plot contagem bases
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
  
  #Plot adapters
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
  
  #Plot ocorrencias
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

  #Tabela sequencias frequentes
  output$tabela_frequencias <- renderTable({req(resultado_qa())
                                            t <- freqSequences(resultado_qa())
                                            t$Arquivo <- input$arquivos$name[match(t$Arquivo,
                                                                            unique(t$Arquivo))]
                                            t$Arquivo <- ifelse(duplicated(t$Arquivo),"",t$Arquivo)
                                            t
  })
  
  #Tabela adapters
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

shinyApp(ui = ui, server = server)