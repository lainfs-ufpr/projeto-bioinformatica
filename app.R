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

source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)