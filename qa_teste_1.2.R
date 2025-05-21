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
    color: black !important;
    font-weight: bold;
    font-size: 16px;
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
      tableOutput("conteudoArquivos"),
      htmlOutput("qa_output")  # precisa ser htmlOutput pra funcionar com shinyjs
    )
  )
)


# Lógica do servidor
server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)
  
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
    
    # QC
    tempo_inicio <- Sys.time()
    qaSummary <- qa(fls, type = "fastq")
    report_path <- report(qaSummary)
    browseURL(report_path)
    tempo_fim <- Sys.time()
    tempo_execucao <- tempo_fim - tempo_inicio
    
    # Atualiza resultado final
    shinyjs::html("qa_output",
                  paste0('<b style="color:black">Tempo de execução:</b> ',
                         round(tempo_execucao, 2), ' segundos'))
    
  })
}

# Inicializa o app Shiny
shinyApp(ui = ui, server = server)
