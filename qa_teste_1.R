library(shiny)
library(shinyFiles)
library(ShortRead)

ui <- fluidPage(
  titlePanel("Controle de Qualidade de Arquivos de Sequenciamento"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("directory", "Selecionar Pasta", "Escolha a pasta de input"),
      verbatimTextOutput("dir_path"),
      actionButton("run_analysis", "Rodar QA")
    ),
    mainPanel(
      verbatimTextOutput("qa_output")
    )
  )
)

server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")  # Defina os volumes disponíveis
  shinyDirChoose(input, "directory", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)
  
  observeEvent(input$directory, {
    path <- parseDirPath(volumes, input$directory)
    dir_path(path)
  })
  
  output$dir_path <- renderText({
    req(dir_path())
    paste("Pasta selecionada:", dir_path())
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

shinyApp(ui, server)
