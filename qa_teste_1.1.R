library(shiny)
library(shinyFiles)
library(ShortRead)

ui <- fluidPage(
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
      verbatimTextOutput("qa_output")
    )
  )
)

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
    req(dir_path())
    fls <- dir(dir_path(), pattern = "*fq$", full.names = TRUE)
    
    if (length(fls) == 0) {
      showNotification("Nenhum arquivo FASTQ encontrado na pasta selecionada!", type = "error")
    } else {
      qaSummary <- qa(fls, type = "fastq")
      report_path <- report(qaSummary)
      browseURL(report_path)
    }
  })
}

shinyApp(ui = ui, server = server)
