library(ShortRead)
library(shiny)
library(fs)
library(shinyFiles)
library(shinyjs)
library(ggplot2) 

# Carregar funções prontas e módulos de servidor
source("scripts/carrega_funcoes.R")
source("server_qualidade.R")
source("server_trimagem.R")

# BACK
server <- function(input, output, session) {
  
  html_content <- HTML(paste(readLines("www/apresentacao.html",
                                       encoding = "UTF-8"), collapse = "\n"))
  
  # Cria uma tela inicial como se fosse uma caixa de dialogo
  showModal(modalDialog(title = NULL,
                        html_content,
                        size = "l",
                        easyClose = TRUE,
                        footer = modalButton("Iniciar"))
  )
  
  dir_path <- reactiveVal(NULL)
  resultado_qa <- reactiveVal(NULL)
  
  # O módulo de qualidade recebe o ID e as referências às variáveis de estado
  qualidadeServer("qualidade_id", 
                  r_dir_path = dir_path, 
                  r_resultado_qa = resultado_qa)
  
  # O módulo de trimagem é independente
  trimagemServer("trimagem_id")
  
  observeEvent(input$mode, {
    toggleCssClass(id = "page", class = "dark-mode")
  })
}