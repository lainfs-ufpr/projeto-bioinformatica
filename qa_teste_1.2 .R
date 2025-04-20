# Instale o BiocManager se ainda não estiver instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instale os pacotes
BiocManager::install(c("shiny", "shinyFiles", "ShortRead"), force=TRUE)

library(stats)
library(shiny)
library(shinyFiles)
library(ShortRead)

# Aumenta o limite máximo de upload para 100MB
options(shiny.maxRequestSize = 100 * 1024^2) 

# Interface do usuário
ui <- fluidPage(
  titlePanel("Seleção de Pasta, Upload e Controle de Qualidade"),
  sidebarLayout(
    sidebarPanel(
      # Botão para escolher diretório
      shinyDirButton("diretorio", "Escolher Pasta", "Selecionar"),
      # Mostra o caminho da pasta escolhida
      verbatimTextOutput("caminhoPasta"),
      # Upload de arquivos
      fileInput("arquivos", "Escolha os arquivos", 
                multiple = TRUE, 
                accept = c(".fasta", ".fa", ".fastq", ".fq", "text/plain")),
      # Botão para iniciar análise
      actionButton("run_analysis", "Rodar Controle de Qualidade")
    ),
    mainPanel(
      # Tabela com os arquivos selecionados
      tableOutput("conteudoArquivos"),
      # Saída com tempo de execução
      verbatimTextOutput("qa_output")
    )
  )
)

# Lógica do servidor
server <- function(input, output, session) {
  # Define volumes acessíveis
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)  # Armazena o caminho da pasta selecionada
  
  # Atualiza caminho da pasta quando selecionada
  observeEvent(input$diretorio, {
    path <- parseDirPath(volumes, input$diretorio)
    dir_path(path)
  })
  
  # Mostra caminho da pasta
  output$caminhoPasta <- renderText({
    req(dir_path())
    paste("Pasta selecionada:", dir_path())
  })
  
  # Mostra tabela com arquivos selecionados
  output$conteudoArquivos <- renderTable({
    req(input$arquivos)
    data.frame(Nome_Arquivo = input$arquivos$name, Caminho = input$arquivos$datapath)
  })
  
  # Quando botão for clicado, roda controle de qualidade
  observeEvent(input$run_analysis, {
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
    } else {
      # Início da contagem de tempo
      tempo_inicio <- Sys.time()
      
      # Executa controle de qualidade
      qaSummary <- qa(fls, type = "fastq")
      report_path <- report(qaSummary)
      browseURL(report_path)
      
      # Fim da contagem de tempo
      tempo_fim <- Sys.time()
      tempo_execucao <- tempo_fim - tempo_inicio
      
      # Exibe tempo de execução
      output$qa_output <- renderText({
        paste("Tempo de execução do controle de qualidade:", round(tempo_execucao, 2), "segundos")
      })
    }
  })
}

# Inicializa o app Shiny
shinyApp(ui = ui, server = server)
