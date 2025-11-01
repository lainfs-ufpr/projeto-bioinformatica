library(bslib)
library(shinyjs)
library(shinyFiles)
library(shinyWidgets)

# Carregando UI
source("ui_qualidade.R")
source("ui_trimagem.R")

# FRONT

ui <- fluidPage(
  
  theme = bs_theme(version = 5,
                   bootswatch = "flatly", 
                   primary = "#1a754f", 
                   base_font = font_google("Lexend Deca")),

  useShinyjs(),
  
  # CSS
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "qualidade.css")),
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "trimagem.css")),
  
  
  page_navbar(
    # HEADER
    title = tags$span(tags$img(src = "logo-principal.png", id = "logo-fixo"), 
                      "Controle de Qualidade de Sequências FASTQ", 
                      class = "titulo-app"),
    
    # QUALIDADE 
    nav_panel("Qualidade",
              qualidadeUI("qualidade_id") 
    ), 
    
    # TRIMAGEM
    nav_panel("Trimagem",
              trimagemUI("trimagem_id")
    ), 
    
    nav_spacer(),
    
    nav_item(input_dark_mode(id = "mode", value = FALSE)), 
    
    id = "page"
  )
)