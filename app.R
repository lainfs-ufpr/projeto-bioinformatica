library(shiny)
library(bslib)
library(shinyjs)
library(shinyFiles)
# botão de download, precisa fazer o download 
# dessa biblioteca e do svglite (para baixar svg)
library(shinyWidgets)
library(ShortRead)
library(ggplot2)
library(fs)
library(scales)
library(plotly)
library(viridis)

# Carregar funções prontas
source("scripts/carrega_funcoes.R")

# Tamanho máximo de input
options(shiny.maxRequestSize = 100 * 1024^3)

source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)