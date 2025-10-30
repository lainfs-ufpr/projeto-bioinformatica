library(bslib)
library(shinyjs)
library(shinyFiles)
library(shinyWidgets)

# FRONT

ui <- fluidPage(

  theme = bs_theme(version = 5, bg = "#FFFFFF", fg = "#31231a", 
                   bootswatch = "flatly", primary = "#1a754f", 
                   base_font = font_google("Lexend Deca")),
  useShinyjs(),
  
  # CSS
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  
  
  page_navbar(
    # HEADER
    title = tags$span(tags$img(src = "logo-principal.png", id = "logo-fixo",
                               style = "height:80px; margin-right:10px; vertical-align:middle;"), 
                      "Controle de Qualidade de Sequências FASTQ", 
                      class = "titulo-app"),
    
    # QUALIDADE 
    nav_panel("Qualidade",
              sidebarLayout(
                #sidebar
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
                
                # Paineis de gráficos
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
    ), 
    
    # TRIMAGEM
    nav_panel("Trimagem",
              sidebarLayout(
                sidebarPanel(
                  class = "sidebar-trim",
                  fileInput("fastq_trim", "Selecione arquivo FASTQ", accept = c(".fastq", ".fq")),
                  textInput("adapters", "Arquivo de adaptadores",
                            value = "Trimmomatic-0.39/adapters/TruSeq3-SE.fa"),
                  numericInput("sliding_window", "Tamanho da janela (SLIDINGWINDOW)", 4, min = 1),
                  numericInput("quality_cutoff", "Qualidade mínima", 20, min = 1),
                  numericInput("minlen", "Comprimento mínimo (MINLEN)", 50, min = 1),
                  textInput("output_trim", "Nome do arquivo de saída", "reads_trimmed.fastq"),
                  actionButton("run_trim", "Rodar Trimmomatic", class = "btnQA-custom")
                ),
                mainPanel(
                  div(class = "trim-main-panel",
                      verbatimTextOutput("trim_log"),
                      div(style = "width: fit-content; margin-left: auto; margin-right: auto;",
                          tableOutput("trim_stats")
                      )
                  )
                )
              )
    ), 
    
    nav_spacer(),

    nav_item(input_dark_mode(id = "mode", value = FALSE)), 
    
    id = "page"
  )
)